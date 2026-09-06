#!/usr/bin/env python3
"""Step 1: separate reads containing ambiguous/non-ACGT bases."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from serial2brad16.common import iter_fastq, open_text, strip_fastq_suffix, write_fastq_record


def collect_inputs(values, input_dir: Path):
    if values:
        return [Path(v) for v in values]
    found = []
    for pattern in ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"):
        found.extend(input_dir.glob(pattern))
    return sorted(set(found))


def process(path: Path, outdir: Path, gzip_output: str, write_rejected: bool):
    stem = strip_fastq_suffix(path.name)
    compressed = path.name.lower().endswith(".gz") if gzip_output == "auto" else gzip_output == "yes"
    suffix = ".fastq.gz" if compressed else ".fastq"
    accepted = outdir / f"{stem}.ATCG{suffix}"
    rejected = outdir / f"{stem}.N{suffix}"
    total = kept = removed = 0
    with open_text(accepted, "wt") as good:
        bad = open_text(rejected, "wt") if write_rejected else None
        try:
            for record in iter_fastq(path):
                total += 1
                if record.seq and all(base in "ACGTacgt" for base in record.seq):
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
    parser = argparse.ArgumentParser(description="Serial2BRAD16 step 1: remove FASTQ reads containing N or another non-ACGT base.")
    parser.add_argument("-i", "--input", nargs="+", help="Input FASTQ(.gz) files. If omitted, scan -d.")
    parser.add_argument("-d", "--input-dir", default="rawdata", help="Input directory (historical default: rawdata).")
    parser.add_argument("-o", "--outdir", default="extractNdata", help="Output directory (default: extractNdata).")
    parser.add_argument("--gzip-output", choices=("auto", "yes", "no"), default="auto")
    parser.add_argument("--no-rejected", action="store_true", help="Do not retain rejected .N files.")
    args = parser.parse_args()
    inputs = collect_inputs(args.input, Path(args.input_dir))
    if not inputs:
        parser.error("No FASTQ input was found")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in inputs:
        if not path.is_file():
            parser.error(f"Input does not exist: {path}")
        rows.append((path, *process(path, outdir, args.gzip_output, not args.no_rejected)))
    with open(outdir / "extractN_qc.tsv", "w", encoding="utf-8", newline="") as handle:
        handle.write("input\tATCG_output\trejected_output\tinput_reads\tATCG_reads\trejected_reads\tretention_fraction\n")
        for inp, accepted, rejected, total, kept, removed in rows:
            handle.write(f"{inp}\t{accepted}\t{rejected or ''}\t{total}\t{kept}\t{removed}\t{(kept/total if total else 0):.6f}\n")
    for _, accepted, *_ in rows:
        print(accepted)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
