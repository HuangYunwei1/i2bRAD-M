#!/usr/bin/env python3
"""Step 2: reproduce the training-document quality rule."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from serial2brad16.common import iter_fastq, open_text, quality_passes, strip_fastq_suffix, write_fastq_record


def collect_inputs(values, input_dir: Path):
    if values:
        return [Path(v) for v in values]
    found = []
    for pattern in ("*.ATCG.fastq", "*.ATCG.fq", "*.ATCG.fastq.gz", "*.ATCG.fq.gz"):
        found.extend(input_dir.glob(pattern))
    return sorted(set(found))


def process(path: Path, outdir: Path, min_q: int, min_percent: float, phred_offset: int, gzip_output: str, write_rejected: bool):
    stem = strip_fastq_suffix(path.name)
    if stem.endswith(".ATCG"):
        stem = stem[:-5]
    compressed = path.name.lower().endswith(".gz") if gzip_output == "auto" else gzip_output == "yes"
    suffix = ".fastq.gz" if compressed else ".fastq"
    accepted = outdir / f"{stem}.GoodQuality{suffix}"
    rejected = outdir / f"{stem}.PoorQuality{suffix}"
    total = kept = removed = 0
    with open_text(accepted, "wt") as good:
        bad = open_text(rejected, "wt") if write_rejected else None
        try:
            for record in iter_fastq(path):
                total += 1
                if quality_passes(record.qual, min_q, min_percent, phred_offset):
                    kept += 1
                    write_fastq_record(good, record)
                else:
                    removed += 1
                    if bad:
                        write_fastq_record(bad, record)
        finally:
            if bad:
                bad.close()
    return accepted, rejected if write_rejected else None, total, kept, removed


def main() -> int:
    parser = argparse.ArgumentParser(description="Serial2BRAD16 step 2: retain reads whose required proportion of bases has quality strictly greater than the threshold.")
    parser.add_argument("-i", "--input", nargs="+", help="Input .ATCG FASTQ(.gz) files; if omitted, scan -d.")
    parser.add_argument("-d", "--input-dir", default="extractNdata")
    parser.add_argument("-o", "--outdir", default="goodqualitydata")
    parser.add_argument("--min-q", type=int, default=20, help="Base quality must be > this value (default: 20).")
    parser.add_argument("--min-percent", type=float, default=80.0, help="Required passing-base percentage (default: 80).")
    parser.add_argument("--phred-offset", type=int, choices=(33, 64), default=33)
    parser.add_argument("--gzip-output", choices=("auto", "yes", "no"), default="auto")
    parser.add_argument("--no-rejected", action="store_true")
    args = parser.parse_args()
    if not 0 <= args.min_percent <= 100:
        parser.error("--min-percent must be within 0-100")
    inputs = collect_inputs(args.input, Path(args.input_dir))
    if not inputs:
        parser.error("No .ATCG FASTQ input was found")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in inputs:
        if not path.is_file():
            parser.error(f"Input does not exist: {path}")
        rows.append((path, *process(path, outdir, args.min_q, args.min_percent, args.phred_offset, args.gzip_output, not args.no_rejected)))
    with open(outdir / "goodquality_qc.tsv", "w", encoding="utf-8", newline="") as handle:
        handle.write("input\tgood_output\tpoor_output\tinput_reads\tgood_reads\tpoor_reads\tretention_fraction\tmin_q_strictly_greater\tmin_percent\tphred_offset\n")
        for inp, accepted, rejected, total, kept, removed in rows:
            handle.write(f"{inp}\t{accepted}\t{rejected or ''}\t{total}\t{kept}\t{removed}\t{(kept/total if total else 0):.6f}\t{args.min_q}\t{args.min_percent}\t{args.phred_offset}\n")
    for _, accepted, *_ in rows:
        print(accepted)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
