#!/usr/bin/env python3
"""Batch-run Kraken2 on paired RNA reads and summarize species abundances."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_VERSION = "0.1.0"
SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")


@dataclass(frozen=True)
class Sample:
    sample_id: str
    read1: Path
    read2: Path


@dataclass
class ReportSummary:
    species_counts: dict[str, int]
    species_taxids: dict[str, str]
    classified_fragments: int
    unclassified_fragments: int

    @property
    def total_fragments(self) -> int:
        return self.classified_fragments + self.unclassified_fragments

    @property
    def species_level_fragments(self) -> int:
        return sum(self.species_counts.values())


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run Kraken2 for paired metatranscriptomic reads and create "
            "species-level count and relative-abundance matrices."
        )
    )
    parser.add_argument("-i", required=True, type=Path, help="RNA sample sheet (TSV)")
    parser.add_argument("-d", required=True, type=Path, help="Kraken2 database directory")
    parser.add_argument(
        "-o", type=Path, default=Path("RNA_profile"),
        help="output directory [RNA_profile]",
    )
    parser.add_argument("-p", type=int, default=8, help="Kraken2 threads [8]")
    parser.add_argument(
        "-c", type=float, default=0.0,
        help="Kraken2 confidence threshold, from 0 to 1 [0.0]",
    )
    return parser.parse_args(argv)


def normalize_species_name(name: str) -> str:
    """Normalize whitespace only; taxonomic synonyms are not inferred."""
    normalized = re.sub(r"\s+", "_", name.strip())
    normalized = re.sub(r"_+", "_", normalized)
    if not normalized:
        raise ValueError("Encountered an empty species name in a Kraken2 report.")
    return normalized


def resolve_read_path(text: str, base_dir: Path, line_no: int, field: str) -> Path:
    value = text.strip()
    if not value:
        raise ValueError(f"Sample sheet line {line_no}: empty {field} value.")
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base_dir / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(
            f"Sample sheet line {line_no}: {field} file not found: {path}"
        )
    return path


def load_samples(path: Path) -> list[Sample]:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"RNA sample sheet not found: {path}")

    samples: list[Sample] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Sample sheet has no header: {path}")
        reader.fieldnames = [str(x).strip() for x in reader.fieldnames]
        required = {"sample_id", "read1", "read2"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(
                "Sample sheet is missing required column(s): " + ", ".join(missing)
            )

        for line_no, row in enumerate(reader, start=2):
            extra_fields = row.get(None)
            if extra_fields and any(str(value).strip() for value in extra_fields):
                raise ValueError(
                    f"Sample sheet line {line_no}: more columns than declared in the header."
                )
            known_values = [row.get(field) or "" for field in reader.fieldnames]
            if all(not value.strip() for value in known_values):
                continue
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                raise ValueError(f"Sample sheet line {line_no}: empty sample_id.")
            if not SAMPLE_ID_RE.fullmatch(sample_id):
                raise ValueError(
                    f"Sample sheet line {line_no}: invalid sample_id {sample_id!r}; "
                    "use letters, numbers, '.', '_' or '-'."
                )
            if sample_id in seen:
                raise ValueError(f"Duplicate sample_id: {sample_id}")
            seen.add(sample_id)

            read1 = resolve_read_path(row.get("read1") or "", path.parent, line_no, "read1")
            read2 = resolve_read_path(row.get("read2") or "", path.parent, line_no, "read2")
            if read1 == read2:
                raise ValueError(
                    f"Sample sheet line {line_no}: read1 and read2 are the same file."
                )
            if read1.name.lower().endswith(".gz") != read2.name.lower().endswith(".gz"):
                raise ValueError(
                    f"Sample {sample_id}: read1 and read2 must both be gzip-compressed "
                    "or both be uncompressed."
                )
            samples.append(Sample(sample_id, read1, read2))

    if not samples:
        raise ValueError("RNA sample sheet contains no samples.")
    return samples


def prepare_output_directory(path: Path) -> Path:
    path = path.expanduser().resolve()
    if path.exists():
        if not path.is_dir():
            raise ValueError(f"Output path is not a directory: {path}")
        if any(path.iterdir()):
            raise FileExistsError(
                f"Output directory is not empty: {path}. Use a new or empty -o directory."
            )
    else:
        path.mkdir(parents=True)
    for child in ("kraken_reports", "kraken_output", "logs"):
        (path / child).mkdir()
    return path


def parse_count(text: str, path: Path, line_no: int, label: str) -> int:
    try:
        value = int(text)
    except ValueError as exc:
        raise ValueError(
            f"{path}:{line_no}: invalid integer for {label}: {text!r}"
        ) from exc
    if value < 0:
        raise ValueError(f"{path}:{line_no}: negative {label}: {value}")
    return value


def parse_kraken_report(path: Path) -> ReportSummary:
    """Parse the standard six-column Kraken2 report.

    Species counts are column 2 clade counts for exact rank code ``S``. A
    species clade count includes fragments assigned to subordinate strain or
    subspecies nodes, while subordinate rows are not counted separately.
    """
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Kraken2 report not found: {path}")

    species_counts: dict[str, int] = {}
    species_taxids: dict[str, str] = {}
    taxid_to_species: dict[str, str] = {}
    classified_fragments: int | None = None
    unclassified_fragments = 0

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t", 5)
            if len(fields) != 6:
                raise ValueError(
                    f"{path}:{line_no}: expected 6 tab-separated Kraken2 report "
                    f"fields, found {len(fields)}."
                )
            _percentage, clade_text, direct_text, rank_code, taxid, name = fields
            clade_count = parse_count(clade_text.strip(), path, line_no, "clade count")
            parse_count(direct_text.strip(), path, line_no, "direct count")
            rank_code = rank_code.strip()
            taxid = taxid.strip()

            if rank_code == "U":
                unclassified_fragments += clade_count
                continue
            if rank_code == "R":
                if classified_fragments is not None:
                    raise ValueError(f"{path}:{line_no}: multiple root (R) rows found.")
                classified_fragments = clade_count
                continue
            if rank_code != "S":
                continue

            species = normalize_species_name(name)
            old_species = taxid_to_species.get(taxid)
            if old_species is not None and old_species != species:
                raise ValueError(
                    f"{path}:{line_no}: TaxID {taxid} maps to both "
                    f"{old_species!r} and {species!r}."
                )
            old_taxid = species_taxids.get(species)
            if old_taxid is not None and old_taxid != taxid:
                raise ValueError(
                    f"{path}:{line_no}: normalized name {species!r} maps to "
                    f"multiple TaxIDs ({old_taxid}, {taxid}); refusing to merge them."
                )
            if species in species_counts:
                raise ValueError(
                    f"{path}:{line_no}: duplicate species row for {species!r}."
                )
            taxid_to_species[taxid] = species
            species_taxids[species] = taxid
            species_counts[species] = clade_count

    if classified_fragments is None:
        classified_fragments = 0
    return ReportSummary(
        species_counts=species_counts,
        species_taxids=species_taxids,
        classified_fragments=classified_fragments,
        unclassified_fragments=unclassified_fragments,
    )


def get_kraken_version(executable: str) -> str:
    result = subprocess.run(
        [executable, "--version"],
        text=True,
        capture_output=True,
        check=False,
    )
    text = (result.stdout or result.stderr).strip()
    if result.returncode != 0:
        raise RuntimeError(
            f"Could not obtain Kraken2 version (exit {result.returncode}): {text}"
        )
    return text.splitlines()[0] if text else "unknown"


def shell_join(command: list[str]) -> str:
    return shlex.join(command)


def run_kraken_sample(
    executable: str,
    database: Path,
    sample: Sample,
    output_dir: Path,
    threads: int,
    confidence: float,
) -> tuple[ReportSummary, list[str]]:
    report_path = output_dir / "kraken_reports" / f"{sample.sample_id}.report.tsv"
    classification_path = output_dir / "kraken_output" / f"{sample.sample_id}.kraken.tsv"
    log_path = output_dir / "logs" / f"{sample.sample_id}.log"

    command = [
        executable,
        "--db", str(database),
        "--threads", str(threads),
        "--confidence", str(confidence),
        "--paired",
        "--report", str(report_path),
        "--output", str(classification_path),
    ]
    if sample.read1.name.lower().endswith(".gz"):
        command.append("--gzip-compressed")
    command.extend([str(sample.read1), str(sample.read2)])

    result = subprocess.run(command, text=True, capture_output=True, check=False)
    with log_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"command\t{shell_join(command)}\n")
        handle.write(f"exit_code\t{result.returncode}\n")
        handle.write("\n[stdout]\n")
        handle.write(result.stdout or "")
        handle.write("\n[stderr]\n")
        handle.write(result.stderr or "")

    if result.returncode != 0:
        raise RuntimeError(
            f"Kraken2 failed for sample {sample.sample_id} with exit code "
            f"{result.returncode}. See {log_path}"
        )
    if not report_path.is_file():
        raise RuntimeError(
            f"Kraken2 finished for {sample.sample_id}, but report was not created: {report_path}"
        )
    if not classification_path.is_file():
        raise RuntimeError(
            f"Kraken2 finished for {sample.sample_id}, but classification output "
            f"was not created: {classification_path}"
        )

    return parse_kraken_report(report_path), command


def format_number(value: float | int) -> str:
    if isinstance(value, int):
        return str(value)
    if not math.isfinite(value):
        raise ValueError("Attempted to write a non-finite abundance value.")
    return f"{value:.12g}"


def write_matrix(
    path: Path,
    species_order: list[str],
    sample_order: list[str],
    data: dict[str, dict[str, float | int]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Species", *sample_order])
        for species in species_order:
            writer.writerow(
                [species]
                + [format_number(data.get(sample, {}).get(species, 0)) for sample in sample_order]
            )


def safe_fraction(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator > 0 else 0.0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.p < 1:
            raise ValueError("-p must be at least 1.")
        if not 0.0 <= args.c <= 1.0:
            raise ValueError("-c must be between 0 and 1.")

        database = args.d.expanduser().resolve()
        if not database.is_dir():
            raise FileNotFoundError(f"Kraken2 database directory not found: {database}")
        executable = shutil.which("kraken2")
        if executable is None:
            raise FileNotFoundError(
                "Kraken2 executable was not found in PATH. Activate/install Kraken2 first."
            )

        samples = load_samples(args.i)
        output_dir = prepare_output_directory(args.o)
        kraken_version = get_kraken_version(executable)

        summaries: dict[str, ReportSummary] = {}
        commands: dict[str, list[str]] = {}
        global_species_taxids: dict[str, str] = {}
        global_taxid_species: dict[str, str] = {}
        for sample in samples:
            print(f"[MTKrakenProfiler] Running {sample.sample_id}...", file=sys.stderr)
            summary, command = run_kraken_sample(
                executable, database, sample, output_dir, args.p, args.c
            )
            for species, taxid in summary.species_taxids.items():
                old_taxid = global_species_taxids.get(species)
                if old_taxid is not None and old_taxid != taxid:
                    raise ValueError(
                        f"Across reports, normalized species {species!r} maps to "
                        f"multiple TaxIDs ({old_taxid}, {taxid})."
                    )
                old_species = global_taxid_species.get(taxid)
                if old_species is not None and old_species != species:
                    raise ValueError(
                        f"Across reports, TaxID {taxid} maps to both "
                        f"{old_species!r} and {species!r}."
                    )
                global_species_taxids[species] = taxid
                global_taxid_species[taxid] = species
            summaries[sample.sample_id] = summary
            commands[sample.sample_id] = command

        sample_order = [sample.sample_id for sample in samples]
        species_order = sorted(
            {species for summary in summaries.values() for species in summary.species_counts}
        )
        counts = {
            sample_id: dict(summary.species_counts)
            for sample_id, summary in summaries.items()
        }
        abundances: dict[str, dict[str, float]] = {}
        for sample_id, summary in summaries.items():
            denominator = summary.species_level_fragments
            abundances[sample_id] = {
                species: (count / denominator if denominator > 0 else 0.0)
                for species, count in summary.species_counts.items()
            }

        write_matrix(
            output_dir / "RNA_species_counts.tsv",
            species_order,
            sample_order,
            counts,
        )
        write_matrix(
            output_dir / "RNA_species_abundance.tsv",
            species_order,
            sample_order,
            abundances,
        )

        with (output_dir / "sample_qc.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "sample_id",
                    "total_fragments",
                    "classified_fragments",
                    "unclassified_fragments",
                    "species_level_fragments",
                    "detected_species",
                    "classified_fraction",
                    "species_level_fraction",
                ]
            )
            for sample_id in sample_order:
                summary = summaries[sample_id]
                writer.writerow(
                    [
                        sample_id,
                        summary.total_fragments,
                        summary.classified_fragments,
                        summary.unclassified_fragments,
                        summary.species_level_fragments,
                        sum(count > 0 for count in summary.species_counts.values()),
                        format_number(
                            safe_fraction(summary.classified_fragments, summary.total_fragments)
                        ),
                        format_number(
                            safe_fraction(summary.species_level_fragments, summary.total_fragments)
                        ),
                    ]
                )

        config = {
            "script": "MTKrakenProfiler.py",
            "script_version": SCRIPT_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "sample_sheet": str(args.i.expanduser().resolve()),
            "database": str(database),
            "output_directory": str(output_dir),
            "kraken2_executable": executable,
            "kraken2_version": kraken_version,
            "threads": args.p,
            "confidence": args.c,
            "paired_end": True,
            "species_rank_code": "S",
            "species_count_field": "Kraken2 report column 2 (clade fragments)",
            "abundance_denominator": (
                "sum of clade fragment counts across exact species-rank (S) rows "
                "within each sample"
            ),
            "name_normalization": "trim whitespace; replace internal whitespace with underscore",
            "commands": {sample: shell_join(command) for sample, command in commands.items()},
        }
        with (output_dir / "config.json").open("w", encoding="utf-8") as handle:
            json.dump(config, handle, ensure_ascii=False, indent=2)
            handle.write("\n")

        print(
            f"[MTKrakenProfiler] Completed {len(samples)} sample(s); "
            f"wrote {len(species_order)} species to {output_dir}",
            file=sys.stderr,
        )
        return 0
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

