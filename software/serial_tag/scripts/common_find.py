#!/usr/bin/env python3
"""Step 3: retain complete pairs after independent filtering."""
from __future__ import annotations
import argparse
import json
import sqlite3
import sys
import tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from serial2brad16.common import (
    FastqRecord, discover_fastq_pairs, infer_fastq_identity, iter_fastq,
    open_text, safe_filename, write_fastq_record,
)


def output_paths(outdir: Path, prefix: str, compressed: bool):
    suffix = ".fastq.gz" if compressed else ".fastq"
    return (
        outdir / f"{prefix}_R1.GoodQuality.filter{suffix}",
        outdir / f"{prefix}_R2.GoodQuality.filter{suffix}",
    )


def memory_pair(forward: Path, reverse: Path, out1: Path, out2: Path):
    reverse_records = {}
    reverse_total = 0
    for rec in iter_fastq(reverse):
        reverse_total += 1
        if rec.read_id in reverse_records:
            raise ValueError(f"Duplicate read ID in reverse input: {rec.read_id}")
        reverse_records[rec.read_id] = rec
    forward_total = paired = 0
    seen_forward = set()
    with open_text(out1, "wt") as h1, open_text(out2, "wt") as h2:
        for rec1 in iter_fastq(forward):
            forward_total += 1
            if rec1.read_id in seen_forward:
                raise ValueError(f"Duplicate read ID in forward input: {rec1.read_id}")
            seen_forward.add(rec1.read_id)
            rec2 = reverse_records.get(rec1.read_id)
            if rec2 is not None:
                write_fastq_record(h1, rec1)
                write_fastq_record(h2, rec2)
                paired += 1
    return forward_total, reverse_total, paired


def sqlite_pair(forward: Path, reverse: Path, out1: Path, out2: Path, temp_dir: Path | None):
    temp_parent = str(temp_dir) if temp_dir else None
    with tempfile.TemporaryDirectory(prefix="serial2brad_pairs_", dir=temp_parent) as td:
        db_path = Path(td) / "pairs.sqlite3"
        con = sqlite3.connect(db_path)
        try:
            con.execute("PRAGMA journal_mode=OFF")
            con.execute("PRAGMA synchronous=OFF")
            con.execute("CREATE TABLE reverse_reads (id TEXT PRIMARY KEY, payload TEXT NOT NULL)")
            con.execute("CREATE TABLE forward_seen (id TEXT PRIMARY KEY)")
            reverse_total = 0
            try:
                batch = []
                for rec in iter_fastq(reverse):
                    reverse_total += 1
                    batch.append((rec.read_id, json.dumps([rec.header, rec.seq, rec.plus, rec.qual], separators=(",", ":"))))
                    if len(batch) >= 10000:
                        con.executemany("INSERT INTO reverse_reads VALUES (?,?)", batch)
                        batch.clear()
                if batch:
                    con.executemany("INSERT INTO reverse_reads VALUES (?,?)", batch)
                con.commit()
            except sqlite3.IntegrityError as exc:
                raise ValueError("Duplicate read ID in reverse input") from exc

            forward_total = paired = 0
            with open_text(out1, "wt") as h1, open_text(out2, "wt") as h2:
                for rec1 in iter_fastq(forward):
                    forward_total += 1
                    try:
                        con.execute("INSERT INTO forward_seen VALUES (?)", (rec1.read_id,))
                    except sqlite3.IntegrityError as exc:
                        raise ValueError(f"Duplicate read ID in forward input: {rec1.read_id}") from exc
                    row = con.execute("SELECT payload FROM reverse_reads WHERE id=?", (rec1.read_id,)).fetchone()
                    if row:
                        rec2 = FastqRecord(*json.loads(row[0]))
                        write_fastq_record(h1, rec1)
                        write_fastq_record(h2, rec2)
                        paired += 1
            con.commit()
            return forward_total, reverse_total, paired
        finally:
            con.close()


def infer_single_prefix(forward: Path, reverse: Path) -> str:
    left = infer_fastq_identity(forward, "goodquality")
    right = infer_fastq_identity(reverse, "goodquality")
    if left is None or right is None or left[1] != 1 or right[1] != 2 or left[0].casefold() != right[0].casefold():
        raise ValueError("-f/-b filenames must form one identifiable GoodQuality R1/R2 pair; otherwise supply --prefix")
    return left[0]


def process_pair(prefix: str, forward: Path, reverse: Path, outdir: Path, backend: str,
                 temp_dir: Path | None, compressed: bool):
    out1, out2 = output_paths(outdir, prefix, compressed)
    if backend == "memory":
        totals = memory_pair(forward, reverse, out1, out2)
    else:
        totals = sqlite_pair(forward, reverse, out1, out2, temp_dir)
    f_total, r_total, paired = totals
    qc = outdir / f"{prefix}.pair_qc.tsv"
    with open(qc, "w", encoding="utf-8") as handle:
        handle.write("library_id\tforward_input\treverse_input\tforward_reads\treverse_reads\tpaired_reads\tforward_retention_fraction\treverse_retention_fraction\tbackend\n")
        handle.write(f"{prefix}\t{forward}\t{reverse}\t{f_total}\t{r_total}\t{paired}\t{(paired/f_total if f_total else 0):.6f}\t{(paired/r_total if r_total else 0):.6f}\t{backend}\n")
    print(out1)
    print(out2)
    return prefix, forward, reverse, out1, out2, f_total, r_total, paired


def main() -> int:
    parser = argparse.ArgumentParser(description="Serial2BRAD16 step 3: batch-discover GoodQuality R1/R2 files and retain only complete pairs.")
    parser.add_argument("-d", "--input-dir", help="Directory mode: scan all GoodQuality R1/R2 files (recommended: goodqualitydata).")
    parser.add_argument("-f", "--forward", help="Expert compatibility mode: one forward GoodQuality FASTQ.")
    parser.add_argument("-b", "--reverse", help="Expert compatibility mode: one reverse GoodQuality FASTQ.")
    parser.add_argument("-o", "--outdir", default="goodquality-filter")
    parser.add_argument("--prefix", help="Expert mode library prefix; normally inferred from filenames.")
    parser.add_argument("--backend", choices=("sqlite", "memory"), default="sqlite")
    parser.add_argument("--temp-dir", help="Parent directory for the temporary SQLite database.")
    parser.add_argument("--gzip-output", choices=("yes", "no"), default="yes")
    args = parser.parse_args()

    if args.input_dir and (args.forward or args.reverse):
        parser.error("Use either directory mode (-d) or one-pair mode (-f/-b), not both")
    if bool(args.forward) != bool(args.reverse):
        parser.error("-f and -b must be supplied together")
    try:
        if args.forward:
            forward, reverse = Path(args.forward), Path(args.reverse)
            if not forward.is_file() or not reverse.is_file():
                parser.error("Both forward and reverse inputs must exist")
            prefix = safe_filename(args.prefix) if args.prefix else infer_single_prefix(forward, reverse)
            pairs = [(prefix, forward, reverse)]
        else:
            pairs = discover_fastq_pairs(Path(args.input_dir or "goodqualitydata"), "goodquality")
    except ValueError as exc:
        parser.error(str(exc))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(args.temp_dir) if args.temp_dir else None
    if temp_dir:
        temp_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for prefix, forward, reverse in pairs:
        rows.append(process_pair(prefix, forward, reverse, outdir, args.backend, temp_dir, args.gzip_output == "yes"))

    with open(outdir / "common_find_batch_qc.tsv", "w", encoding="utf-8") as handle:
        handle.write("library_id\tforward_input\treverse_input\tforward_output\treverse_output\tforward_reads\treverse_reads\tpaired_reads\n")
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")
    print(f"Completed {len(rows)} librar{'y' if len(rows) == 1 else 'ies'}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
