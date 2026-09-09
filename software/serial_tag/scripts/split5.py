#!/usr/bin/env python3
"""Step 6: create five overlapping enzyme-specific candidate windows."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from serial2brad16.common import get_enzyme, iter_fasta, list_assembled_fasta, load_config, open_text, strip_fasta_suffix, write_fasta_record


def parse_windows(text: str):
    windows = []
    for item in text.split(","):
        start, end = item.split("-", 1)
        windows.append((int(start), int(end)))
    if len(windows) != 5 or any(start < 0 or end < start for start, end in windows):
        raise ValueError("exactly five valid inclusive 0-based start-end pairs are required")
    return windows


def split_one(inp: Path, outdir: Path, prefix: str, windows, compressed: bool):
    suffix = ".fasta.gz" if compressed else ".fasta"
    outputs = [outdir / f"{prefix}_T{i}{suffix}" for i in range(1, 6)]
    handles = [open_text(path, "wt") for path in outputs]
    total, truncated = 0, [0] * 5
    try:
        for header, seq in iter_fasta(inp):
            total += 1
            clean_header = header[1:] if header.startswith(">") else header
            base = f"{prefix}|{clean_header}"
            for index, (start, end_inclusive) in enumerate(windows):
                candidate = seq[start:end_inclusive + 1]
                if len(seq) <= end_inclusive:
                    truncated[index] += 1
                write_fasta_record(handles[index], f"{base}|T{index+1}|window0={start}-{end_inclusive}", candidate)
    finally:
        for handle in handles:
            handle.close()
    qc_path = outdir / f"{prefix}.split5_qc.tsv"
    with open(qc_path, "w", encoding="utf-8") as handle:
        handle.write("position\tstart_0based\tend_0based_inclusive\twindow_length\tinput_records\ttruncated_windows\n")
        for index, (start, end) in enumerate(windows):
            handle.write(f"T{index+1}\t{start}\t{end}\t{end-start+1}\t{total}\t{truncated[index]}\n")
    for output in outputs:
        print(output)
    return inp, prefix, total, sum(truncated), outputs


def main() -> int:
    parser = argparse.ArgumentParser(description="Serial2BRAD16 step 6: batch-split each assembled read into five overlapping candidate windows; this is not final fixed-length extraction.")
    parser.add_argument("-d", "--input-dir", help="Directory mode (recommended: 5tagfasta).")
    parser.add_argument("-i", "--input", help="One assembled FASTA(.gz) file.")
    parser.add_argument("-e", "--enzyme", required=True, type=int, choices=range(1, 17), help="One global enzyme preset ID (1-16 only).")
    parser.add_argument("-o", "--outdir", default="split5Data")
    parser.add_argument("--prefix", help="Optional prefix override for one -i file.")
    parser.add_argument("--windows", help="Optional five inclusive 0-based ranges, e.g. 0-49,32-81,65-114,98-147,131-190")
    parser.add_argument("--config")
    parser.add_argument("--gzip-output", choices=("yes", "no"), default="no")
    args = parser.parse_args()

    if args.input_dir and args.input:
        parser.error("Use either directory mode (-d) or single-file mode (-i), not both")
    if args.prefix and not args.input:
        parser.error("--prefix is valid only with one -i input")
    if args.input:
        inputs = [Path(args.input)]
        if not inputs[0].is_file():
            parser.error(f"Input does not exist: {inputs[0]}")
    else:
        try:
            inputs = list_assembled_fasta(Path(args.input_dir or "5tagfasta"))
        except ValueError as exc:
            parser.error(str(exc))
        if not inputs:
            parser.error("No .assembled.fasta(.gz) files were found")

    config = load_config(args.config)
    enzyme = get_enzyme(args.enzyme, config)
    try:
        windows = parse_windows(args.windows) if args.windows else [tuple(value) for value in enzyme["windows"]]
    except Exception as exc:
        parser.error(str(exc))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    seen = set()
    for inp in inputs:
        prefix = args.prefix or strip_fasta_suffix(inp.name).removesuffix(".assembled")
        key = prefix.casefold()
        if key in seen:
            parser.error(f"Duplicate output prefix detected: {prefix}")
        seen.add(key)
        rows.append(split_one(inp, outdir, prefix, windows, args.gzip_output == "yes"))
    with open(outdir / "split5_batch_qc.tsv", "w", encoding="utf-8") as handle:
        handle.write("input\tlibrary_id\tinput_records\ttruncated_windows_total\toutput_T1\toutput_T2\toutput_T3\toutput_T4\toutput_T5\n")
        for inp, prefix, total, truncated, outputs in rows:
            handle.write(f"{inp}\t{prefix}\t{total}\t{truncated}\t" + "\t".join(map(str, outputs)) + "\n")
    print(f"Split {len(rows)} assembled FASTA file{'s' if len(rows) != 1 else ''} into five positions each.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
