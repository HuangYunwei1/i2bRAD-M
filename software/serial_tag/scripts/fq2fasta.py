#!/usr/bin/env python3
"""Step 5: FASTQ to FASTA conversion."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from serial2brad16.common import iter_fastq, list_assembled_fastq, open_text, strip_fastq_suffix, write_fasta_record


def convert(inp: Path, out: Path):
    out.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open_text(out, "wt") as handle:
        for record in iter_fastq(inp):
            count += 1
            write_fasta_record(handle, record.header, record.seq.upper())
    qc = out.parent / f"{strip_fastq_suffix(out.name)}.fq2fasta_qc.tsv"
    with open(qc, "w", encoding="utf-8") as handle:
        handle.write("input\toutput\trecords\n")
        handle.write(f"{inp}\t{out}\t{count}\n")
    print(out)
    return inp, out, count


def main() -> int:
    parser = argparse.ArgumentParser(description="Serial2BRAD16 step 5: batch-convert PEAR .assembled.fastq files to FASTA.")
    parser.add_argument("-d", "--input-dir", help="Directory mode (recommended: 5tagData).")
    parser.add_argument("-i", "--input", help="Expert compatibility mode: one assembled FASTQ(.gz).")
    parser.add_argument("-o", "--output", help="Output FASTA(.gz); valid only with one -i input.")
    parser.add_argument("--outdir", default="5tagfasta")
    parser.add_argument("--gzip-output", choices=("yes", "no"), default="no")
    args = parser.parse_args()

    if args.input_dir and args.input:
        parser.error("Use either directory mode (-d) or single-file mode (-i), not both")
    if args.output and not args.input:
        parser.error("-o/--output is valid only with -i/--input")
    if args.input:
        inputs = [Path(args.input)]
        if not inputs[0].is_file():
            parser.error(f"Input does not exist: {inputs[0]}")
    else:
        try:
            inputs = list_assembled_fastq(Path(args.input_dir or "5tagData"))
        except ValueError as exc:
            parser.error(str(exc))
        if not inputs:
            parser.error("No .assembled.fastq(.gz) files were found")

    outdir = Path(args.outdir)
    suffix = ".fasta.gz" if args.gzip_output == "yes" else ".fasta"
    rows = []
    for inp in inputs:
        out = Path(args.output) if args.output else outdir / f"{strip_fastq_suffix(inp.name)}{suffix}"
        rows.append(convert(inp, out))
    batch_dir = rows[0][1].parent if args.output else outdir
    batch_dir.mkdir(parents=True, exist_ok=True)
    with open(batch_dir / "fq2fasta_batch_qc.tsv", "w", encoding="utf-8") as handle:
        handle.write("input\toutput\trecords\n")
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")
    print(f"Converted {len(rows)} assembled FASTQ file{'s' if len(rows) != 1 else ''}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
