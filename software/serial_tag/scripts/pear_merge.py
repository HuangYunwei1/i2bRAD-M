#!/usr/bin/env python3
"""Step 4: merge filtered mates with PEAR using enzyme-specific broad limits."""
from __future__ import annotations
import argparse
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from serial2brad16.common import discover_fastq_pairs, get_enzyme, infer_fastq_identity, iter_fastq, load_config, safe_filename


def infer_single_prefix(forward: Path, reverse: Path) -> str:
    left = infer_fastq_identity(forward, "filtered")
    right = infer_fastq_identity(reverse, "filtered")
    if left is None or right is None or left[1] != 1 or right[1] != 2 or left[0].casefold() != right[0].casefold():
        raise ValueError("Cannot automatically derive one shared library prefix from -f/-r filenames; supply -op")
    return left[0]


def merge_one(prefix: str, forward: Path, reverse: Path, outdir: Path, enzyme: dict,
              minimum: int, maximum: int, pear: str, threads: int, dry_run: bool):
    output_prefix = outdir / prefix
    command = [pear, "-f", str(forward), "-r", str(reverse), "-m", str(maximum), "-n", str(minimum), "-o", str(output_prefix), "-j", str(threads)]
    print("COMMAND:", " ".join(command))
    if dry_run:
        return prefix, forward, reverse, Path(str(output_prefix) + ".assembled.fastq"), 0

    subprocess.run(command, check=True)
    assembled = Path(str(output_prefix) + ".assembled.fastq")
    if not assembled.is_file():
        raise RuntimeError(f"PEAR completed without expected output: {assembled}")
    lengths = Counter()
    assembled_count = 0
    for record in iter_fastq(assembled):
        assembled_count += 1
        lengths[len(record.seq)] += 1
    with open(outdir / f"{prefix}.assembled_length_hist.tsv", "w", encoding="utf-8") as handle:
        handle.write("length_bp\tassembled_reads\n")
        for length, count in sorted(lengths.items()):
            handle.write(f"{length}\t{count}\n")
    with open(outdir / f"{prefix}.pear_qc.tsv", "w", encoding="utf-8") as handle:
        handle.write("library_id\tenzyme_id\tenzyme\tpear_min\tpear_max\texpected_center\tassembled_reads\n")
        handle.write(f"{prefix}\t{enzyme['id']}\t{enzyme['name']}\t{minimum}\t{maximum}\t{enzyme['expected_assembled_length']}\t{assembled_count}\n")
    print(assembled)
    return prefix, forward, reverse, assembled, assembled_count


def main() -> int:
    parser = argparse.ArgumentParser(description="Serial2BRAD16 step 4: batch PEAR merge; library prefixes are normally collected from R1/R2 filenames.")
    parser.add_argument("-d", "--input-dir", help="Directory mode: scan filtered R1/R2 files (recommended: goodquality-filter).")
    parser.add_argument("-f", "--forward", help="One forward filtered FASTQ file.")
    parser.add_argument("-r", "--reverse", help="One reverse filtered FASTQ file.")
    parser.add_argument("-e", "--enzyme", required=True, type=int, choices=range(1, 17), help="One global enzyme preset ID (1-16 only).")
    parser.add_argument("-o", "--outdir", default="5tagData")
    parser.add_argument("-op", "--output-prefix", "--prefix", dest="prefix", help="Optional override for one input pair; normally omitted and inferred from filenames.")
    parser.add_argument("-p", "--threads", type=int, default=1)
    parser.add_argument("--pear", default="pear", help="PEAR executable or absolute path.")
    parser.add_argument("--pear-min", type=int, help="Custom minimum assembled-read length; supply together with --pear-max.")
    parser.add_argument("--pear-max", type=int, help="Custom maximum assembled-read length; supply together with --pear-min.")
    parser.add_argument("--config")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if args.input_dir and (args.forward or args.reverse):
        parser.error("Use either directory mode (-d) or one-pair mode (-f/-r), not both")
    if bool(args.forward) != bool(args.reverse):
        parser.error("-f and -r must be supplied together")
    if (args.pear_min is None) != (args.pear_max is None):
        parser.error("--pear-min and --pear-max must be supplied together")
    if args.threads < 1:
        parser.error("--threads must be at least 1")

    config = load_config(args.config)
    enzyme = get_enzyme(args.enzyme, config)
    minimum = args.pear_min if args.pear_min is not None else int(enzyme["pear_broad_min"])
    maximum = args.pear_max if args.pear_max is not None else int(enzyme["pear_broad_max"])
    if minimum <= 0 or maximum < minimum:
        parser.error("Invalid PEAR length range")

    try:
        if args.forward:
            forward, reverse = Path(args.forward), Path(args.reverse)
            if not forward.is_file() or not reverse.is_file():
                parser.error("Both paired FASTQ inputs must exist")
            prefix = safe_filename(args.prefix) if args.prefix else infer_single_prefix(forward, reverse)
            pairs = [(prefix, forward, reverse)]
        else:
            pairs = discover_fastq_pairs(Path(args.input_dir or "goodquality-filter"), "filtered")
            if args.prefix and len(pairs) != 1:
                parser.error("-op may only be used when the directory contains exactly one R1/R2 pair")
            if args.prefix:
                pairs = [(safe_filename(args.prefix), pairs[0][1], pairs[0][2])]
    except ValueError as exc:
        parser.error(str(exc))

    if not args.dry_run and not (Path(args.pear).is_file() or shutil.which(args.pear)):
        parser.error(f"PEAR executable not found: {args.pear}")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = [merge_one(prefix, forward, reverse, outdir, enzyme, minimum, maximum, args.pear, args.threads, args.dry_run)
            for prefix, forward, reverse in pairs]
    if not args.dry_run:
        with open(outdir / "pear_batch_qc.tsv", "w", encoding="utf-8") as handle:
            handle.write("library_id\tforward_input\treverse_input\tassembled_output\tassembled_reads\tenzyme_id\tenzyme\tpear_min\tpear_max\n")
            for prefix, forward, reverse, assembled, count in rows:
                handle.write(f"{prefix}\t{forward}\t{reverse}\t{assembled}\t{count}\t{enzyme['id']}\t{enzyme['name']}\t{minimum}\t{maximum}\n")
    print(f"Prepared {len(rows)} PEAR command{'s' if len(rows) != 1 else ''}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
