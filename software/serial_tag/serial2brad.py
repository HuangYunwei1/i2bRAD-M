#!/usr/bin/env python3
"""Controller for five-position serial 2b-RAD preprocessing."""
from __future__ import annotations
import argparse
import gzip
import shlex
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from serial2brad16.common import (get_enzyme, iter_fasta, load_config, open_text,
                                 read_manifest, safe_filename, strip_fastq_suffix)

SCRIPTS = ROOT / "scripts"


def command_text(command):
    return " ".join(shlex.quote(str(item)) for item in command)


def run_step(command, expected, log, resume=False, dry_run=False):
    command = [str(item) for item in command]
    log.write(f"[{datetime.now().isoformat(timespec='seconds')}] COMMAND {command_text(command)}\n")
    log.flush()
    if resume and all(Path(path).is_file() for path in expected):
        log.write("SKIP existing outputs\n")
        log.flush()
        return
    print(command_text(command))
    if not dry_run:
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as exc:
            log.write(f"FAILED exit_code={exc.returncode}\n")
            log.flush()
            raise
        missing = [str(path) for path in expected if not Path(path).is_file()]
        if missing:
            log.write("FAILED missing_outputs=" + ", ".join(missing) + "\n")
            log.flush()
            raise RuntimeError("Expected output missing: " + ", ".join(missing))
        log.write("OK\n")
        log.flush()


def count_fasta(path: Path) -> int:
    return sum(1 for _ in iter_fasta(path))


def print_enzymes(config):
    print("ID\tEnzyme\tTag_bp\tSerial_unit_x\tExpected_bp\tBroad_min\tBroad_max")
    for e in config["enzymes"]:
        print(f"{e['id']}\t{e['name']}\t{e['tag_length']}\t{e['serial_unit_length']}\t{e['expected_assembled_length']}\t{e['pear_broad_min']}\t{e['pear_broad_max']}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Serial2BRAD16: paired FASTQ to five position-resolved Type IIB analytical tags and a two-column downstream sample manifest."
    )
    parser.add_argument("-i", "--manifest", help="TSV: library_id,r1,r2,t1,t2,t3,t4,t5 (no enzyme column).")
    parser.add_argument("-e", "--enzyme", type=int, choices=range(1,17), help="One global enzyme preset ID, 1-16.")
    parser.add_argument("-o", "--outdir", default="Serial2BRAD16_results")
    parser.add_argument("-p", "--threads", type=int, default=1)
    parser.add_argument("--pear", default="pear")
    parser.add_argument("--pear-min", type=int, help="Custom minimum assembled-read length; requires --pear-max.")
    parser.add_argument("--pear-max", type=int, help="Custom maximum assembled-read length; requires --pear-min.")
    parser.add_argument("--min-q", type=int, default=20)
    parser.add_argument("--min-percent", type=float, default=80.0)
    parser.add_argument("--phred-offset", type=int, choices=(33,64), default=33)
    parser.add_argument("--pair-backend", choices=("sqlite", "memory"), default="sqlite")
    parser.add_argument("--hit-policy", choices=("first", "unique", "all"), default="first")
    parser.add_argument("--no-rejected", action="store_true", help="Do not retain N/PoorQuality rejected files.")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--config")
    parser.add_argument("--list-enzymes", action="store_true")
    args = parser.parse_args()

    try:
        config = load_config(args.config)
    except (OSError, ValueError) as exc:
        parser.error(f"Invalid enzyme configuration: {exc}")
    if args.list_enzymes:
        print_enzymes(config)
        return 0
    if not args.manifest or args.enzyme is None:
        parser.error("-i/--manifest and global -e/--enzyme are required")
    if (args.pear_min is None) != (args.pear_max is None):
        parser.error("--pear-min and --pear-max must be supplied together")
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if args.pear_min is not None and (args.pear_min <= 0 or args.pear_max < args.pear_min):
        parser.error("Invalid custom PEAR range")
    if not 0 <= args.min_percent <= 100:
        parser.error("--min-percent must be within 0-100")

    manifest_path = Path(args.manifest).resolve()
    try:
        rows = read_manifest(manifest_path)
        enzyme = get_enzyme(args.enzyme, config)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    config_path = Path(args.config).resolve() if args.config else ROOT / "config" / "enzymes.json"
    safe_libraries = {}
    safe_samples = {}
    for row in rows:
        if not row["r1_path"].is_file() or not row["r2_path"].is_file():
            parser.error(f"Missing FASTQ for library {row['library_id']}: {row['r1_path']} / {row['r2_path']}")
        library_key = safe_filename(row["library_id"])
        previous = safe_libraries.setdefault(library_key, row["library_id"])
        if previous != row["library_id"]:
            parser.error(f"library_id filename collision after sanitizing: {previous!r} and {row['library_id']!r}")
        for pos in ("t1", "t2", "t3", "t4", "t5"):
            sample_key = safe_filename(row[pos])
            previous = safe_samples.setdefault(sample_key, row[pos])
            if previous != row[pos]:
                parser.error(f"sample ID filename collision after sanitizing: {previous!r} and {row[pos]!r}")

    outroot = Path(args.outdir).resolve()
    work = outroot / "work"
    final_root = outroot / "final_tags" / enzyme["name"]
    qc_root = outroot / "qc"
    manifest_root = outroot / "manifests"
    for directory in (work, final_root, qc_root, manifest_root):
        directory.mkdir(parents=True, exist_ok=True)

    temp_tags = defaultdict(list)
    resolved = ["library_id\tr1\tr2\tenzyme_id\tenzyme\tt1\tt2\tt3\tt4\tt5"]
    log_path = outroot / "pipeline.log"
    with open(log_path, "a", encoding="utf-8") as log:
        log.write(f"\n[{datetime.now().isoformat(timespec='seconds')}] START enzyme={enzyme['id']}:{enzyme['name']} broad={enzyme['pear_broad_min']}-{enzyme['pear_broad_max']}\n")
        for row in rows:
            library = safe_filename(row["library_id"])
            resolved.append("\t".join([row["library_id"], str(row["r1_path"]), str(row["r2_path"]), str(enzyme["id"]), enzyme["name"], row["t1"], row["t2"], row["t3"], row["t4"], row["t5"]]))

            d1 = work / "extractN" / library
            d1.mkdir(parents=True, exist_ok=True)
            r1_stem, r2_stem = strip_fastq_suffix(row["r1_path"].name), strip_fastq_suffix(row["r2_path"].name)
            r1_atcg, r2_atcg = d1 / f"{r1_stem}.ATCG.fastq.gz", d1 / f"{r2_stem}.ATCG.fastq.gz"
            cmd = [sys.executable, SCRIPTS/"extractN.py", "-i", row["r1_path"], row["r2_path"], "-o", d1, "--gzip-output", "yes"]
            if args.no_rejected: cmd.append("--no-rejected")
            run_step(cmd, [r1_atcg, r2_atcg], log, args.resume, args.dry_run)

            d2 = work / "goodquality" / library
            d2.mkdir(parents=True, exist_ok=True)
            r1_good, r2_good = d2 / f"{r1_stem}.GoodQuality.fastq.gz", d2 / f"{r2_stem}.GoodQuality.fastq.gz"
            cmd = [sys.executable, SCRIPTS/"goodquality.py", "-i", r1_atcg, r2_atcg, "-o", d2,
                   "--min-q", args.min_q, "--min-percent", args.min_percent, "--phred-offset", args.phred_offset, "--gzip-output", "yes"]
            if args.no_rejected: cmd.append("--no-rejected")
            run_step(cmd, [r1_good, r2_good], log, args.resume, args.dry_run)

            d3 = work / "paired" / library
            d3.mkdir(parents=True, exist_ok=True)
            r1_pair = d3 / f"{library}_R1.GoodQuality.filter.fastq.gz"
            r2_pair = d3 / f"{library}_R2.GoodQuality.filter.fastq.gz"
            cmd = [sys.executable, SCRIPTS/"common_find.py", "-f", r1_good, "-b", r2_good, "-o", d3,
                   "--prefix", library, "--backend", args.pair_backend, "--gzip-output", "yes"]
            run_step(cmd, [r1_pair, r2_pair], log, args.resume, args.dry_run)

            d4 = work / "pear" / library
            d4.mkdir(parents=True, exist_ok=True)
            assembled = d4 / f"{library}.assembled.fastq"
            cmd = [sys.executable, SCRIPTS/"pear_merge.py", "-f", r1_pair, "-r", r2_pair, "-e", enzyme["id"],
                   "-o", d4, "--prefix", library, "-p", args.threads, "--pear", args.pear, "--config", config_path]
            if args.pear_min is not None:
                cmd += ["--pear-min", args.pear_min, "--pear-max", args.pear_max]
            run_step(cmd, [assembled], log, args.resume, args.dry_run)

            d5 = work / "fasta" / library
            d5.mkdir(parents=True, exist_ok=True)
            fasta = d5 / f"{library}.assembled.fasta"
            cmd = [sys.executable, SCRIPTS/"fq2fasta.py", "-i", assembled, "-o", fasta]
            run_step(cmd, [fasta], log, args.resume, args.dry_run)

            d6 = work / "split5" / library
            d6.mkdir(parents=True, exist_ok=True)
            split_files = [d6 / f"{library}_T{i}.fasta" for i in range(1,6)]
            cmd = [sys.executable, SCRIPTS/"split5.py", "-i", fasta, "-e", enzyme["id"], "-o", d6, "--prefix", library, "--config", config_path]
            run_step(cmd, split_files, log, args.resume, args.dry_run)

            d7 = work / "tags" / library
            d7.mkdir(parents=True, exist_ok=True)
            samples = [row[f"t{i}"] for i in range(1,6)]
            for index, sample in enumerate(samples, 1):
                part = d7 / f"{library}_T{index}.{enzyme['name']}.tags.fasta.gz"
                cmd = [sys.executable, SCRIPTS/"recover_sample_tags.py", "-r", split_files[index-1], "-e", enzyme["id"],
                       "-o", part, "--hit-policy", args.hit_policy, "--config", config_path]
                run_step(cmd, [part], log, args.resume, args.dry_run)
                temp_tags[sample].append(part)
        log.write(f"[{datetime.now().isoformat(timespec='seconds')}] STAGES COMPLETE\n")

    if args.dry_run:
        print("Dry run completed; no final aggregation was attempted.")
        return 0

    summary_rows = []
    for sample, parts in sorted(temp_tags.items()):
        output = final_root / f"{safe_filename(sample)}.{enzyme['name']}.fa.gz"
        count = 0
        with gzip.open(output, "wt", encoding="utf-8", newline="") as out:
            for part in parts:
                with open_text(part, "rt") as inp:
                    for line in inp:
                        out.write(line)
                        if line.startswith(">"):
                            count += 1
        summary_rows.append((sample, count, output, len(parts)))

    (outroot / "manifest_resolved.tsv").write_text("\n".join(resolved) + "\n", encoding="utf-8")
    with open(qc_root / "final_tag_summary.tsv", "w", encoding="utf-8") as handle:
        handle.write("sample_id\tenzyme_id\tenzyme\ttag_length_bp\tmerged_position_files\tfinal_tags\tfinal_fasta\n")
        for sample, count, output, part_count in summary_rows:
            handle.write(f"{sample}\t{enzyme['id']}\t{enzyme['name']}\t{enzyme['tag_length']}\t{part_count}\t{count}\t{output}\n")

    map_manifest = manifest_root / f"{enzyme['id']:02d}_{enzyme['name']}.samples.tsv"
    with open(map_manifest, "w", encoding="utf-8") as handle:
        for sample, _count, output, _parts in summary_rows:
            handle.write(f"{sample}\t{output}\n")
    print(f"Completed: {final_root}")
    print(f"Sample manifest: {map_manifest}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())


