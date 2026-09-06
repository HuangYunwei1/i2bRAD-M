#!/usr/bin/env python3
"""Match MAP2B DNA and Kraken2-derived RNA abundances and compute RNA/DNA ratios."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_VERSION = "0.1.0"


@dataclass
class AbundanceTable:
    samples: list[str]
    values: dict[str, dict[str, float]]
    original_names: dict[str, str]
    column_sums: dict[str, float]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Match MAP2B DNA and RNA species relative-abundance tables by sample "
            "and Species, filter by abundance, and compute RNA/DNA ratios."
        )
    )
    parser.add_argument("-d", required=True, type=Path, help="MAP2B Abundance.xls")
    parser.add_argument("-r", required=True, type=Path, help="RNA_species_abundance.tsv")
    parser.add_argument(
        "-o", type=Path, default=Path("RNA_DNA_ratio"),
        help="output directory [RNA_DNA_ratio]",
    )
    parser.add_argument(
        "-a", type=float, default=0.001,
        help="minimum DNA relative abundance; 0 disables the DNA threshold [0.001]",
    )
    parser.add_argument(
        "-b", type=float, default=0.001,
        help="minimum RNA relative abundance; 0 disables the RNA threshold [0.001]",
    )
    return parser.parse_args(argv)


def normalize_species_name(name: str) -> str:
    """Normalize whitespace only; taxonomic synonyms are not inferred."""
    normalized = re.sub(r"\s+", "_", name.strip())
    normalized = re.sub(r"_+", "_", normalized)
    if not normalized:
        raise ValueError("Encountered an empty species name.")
    return normalized


def clean_header(raw_header: list[str], path: Path) -> list[str]:
    header = [cell.strip() for cell in raw_header]
    if not header:
        raise ValueError(f"Table has an empty header: {path}")
    duplicates = sorted({name for name in header if header.count(name) > 1})
    if duplicates:
        raise ValueError(
            f"Duplicate column name(s) in {path}: {', '.join(duplicates)}"
        )
    return header


def parse_abundance(text: str, path: Path, line_no: int, sample: str) -> float:
    value_text = text.strip()
    if value_text == "":
        raise ValueError(f"{path}:{line_no}: empty abundance for sample {sample!r}.")
    try:
        value = float(value_text)
    except ValueError as exc:
        raise ValueError(
            f"{path}:{line_no}: non-numeric abundance {value_text!r} "
            f"for sample {sample!r}."
        ) from exc
    if not math.isfinite(value) or value < 0:
        raise ValueError(
            f"{path}:{line_no}: abundance must be finite and non-negative; "
            f"found {value_text!r} for sample {sample!r}."
        )
    if value > 1.0 + 1e-9:
        raise ValueError(
            f"{path}:{line_no}: abundance {value} for sample {sample!r} is greater "
            "than 1. Inputs must be proportions rather than percentages or raw counts."
        )
    return value


def read_abundance_table(path: Path, table_kind: str) -> AbundanceTable:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{table_kind} abundance table not found: {path}")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = clean_header(next(reader), path)
        except StopIteration as exc:
            raise ValueError(f"Empty abundance table: {path}") from exc

        if "Species" not in header:
            raise ValueError(f"{path} has no required 'Species' column.")
        species_index = header.index("Species")
        if table_kind == "DNA":
            samples = header[species_index + 1 :]
            if not samples:
                raise ValueError(
                    f"{path}: no DNA sample columns occur after the Species column."
                )
        elif table_kind == "RNA":
            if species_index != 0:
                raise ValueError(
                    f"{path}: the RNA abundance table must use Species as its first column."
                )
            samples = header[1:]
            if not samples:
                raise ValueError(f"{path}: no RNA sample columns found.")
        else:
            raise ValueError(f"Unknown table kind: {table_kind}")

        if any(not sample for sample in samples):
            raise ValueError(f"{path}: one or more sample column names are empty.")
        if len(samples) != len(set(samples)):
            raise ValueError(f"{path}: duplicate sample names were found.")

        values: dict[str, dict[str, float]] = {}
        original_names: dict[str, str] = {}
        column_sums = {sample: 0.0 for sample in samples}
        expected_columns = len(header)

        for line_no, row in enumerate(reader, start=2):
            if not row or all(not cell.strip() for cell in row):
                continue
            if len(row) != expected_columns:
                raise ValueError(
                    f"{path}:{line_no}: expected {expected_columns} columns, "
                    f"found {len(row)}."
                )
            original_species = row[species_index].strip()
            species = normalize_species_name(original_species)
            if species in values:
                previous = original_names[species]
                raise ValueError(
                    f"{path}:{line_no}: duplicate species after formatting normalization: "
                    f"{previous!r} and {original_species!r} both become {species!r}."
                )
            original_names[species] = original_species
            sample_values: dict[str, float] = {}
            for sample in samples:
                value = parse_abundance(row[header.index(sample)], path, line_no, sample)
                sample_values[sample] = value
                column_sums[sample] += value
            values[species] = sample_values

    if not values:
        raise ValueError(f"{path}: no species rows found.")
    return AbundanceTable(samples, values, original_names, column_sums)


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
    return path


def format_number(value: float | None) -> str:
    if value is None:
        return "NA"
    if not math.isfinite(value):
        raise ValueError("Attempted to write a non-finite value.")
    return f"{value:.12g}"


def detection_status(dna_value: float, rna_value: float) -> str:
    if dna_value > 0 and rna_value > 0:
        return "both"
    if dna_value > 0:
        return "DNA_only"
    if rna_value > 0:
        return "RNA_only"
    return "missing"


def passes_threshold(value: float, threshold: float) -> bool:
    return threshold == 0 or value >= threshold


def write_unmatched_samples(
    path: Path, dna_samples: list[str], rna_samples: list[str]
) -> None:
    dna_set = set(dna_samples)
    rna_set = set(rna_samples)
    unmatched = sorted(dna_set ^ rna_set)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Sample", "Present_in_DNA", "Present_in_RNA"])
        for sample in unmatched:
            writer.writerow(
                [sample, "yes" if sample in dna_set else "no", "yes" if sample in rna_set else "no"]
            )


def write_unmatched_species(
    path: Path, dna_species: set[str], rna_species: set[str]
) -> None:
    unmatched = sorted(dna_species ^ rna_species)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Species", "Present_in_DNA", "Present_in_RNA"])
        for species in unmatched:
            writer.writerow(
                [
                    species,
                    "yes" if species in dna_species else "no",
                    "yes" if species in rna_species else "no",
                ]
            )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        for flag, threshold in (("-a", args.a), ("-b", args.b)):
            if not math.isfinite(threshold) or not 0.0 <= threshold <= 1.0:
                raise ValueError(f"{flag} must be a finite number between 0 and 1.")

        dna_path = args.d.expanduser().resolve()
        rna_path = args.r.expanduser().resolve()
        dna = read_abundance_table(dna_path, "DNA")
        rna = read_abundance_table(rna_path, "RNA")

        rna_sample_set = set(rna.samples)
        shared_samples = [sample for sample in dna.samples if sample in rna_sample_set]
        if not shared_samples:
            raise ValueError("DNA and RNA abundance tables have no shared sample names.")

        output_dir = prepare_output_directory(args.o)
        write_unmatched_samples(
            output_dir / "unmatched_samples.tsv", dna.samples, rna.samples
        )

        dna_species = set(dna.values)
        rna_species = set(rna.values)
        all_species = sorted(dna_species | rna_species)
        write_unmatched_species(
            output_dir / "unmatched_species.tsv", dna_species, rna_species
        )

        summary_rows: list[list[str | int]] = []
        output_rows = 0
        ratio_path = output_dir / "RNA_DNA_ratio.tsv"
        with ratio_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "Sample",
                    "Species",
                    "DNA_relative_abundance",
                    "RNA_relative_abundance",
                    "RNA_DNA_ratio",
                    "log2_ratio",
                    "Status",
                ]
            )
            for sample in shared_samples:
                retained = 0
                status_counts = {"both": 0, "DNA_only": 0, "RNA_only": 0, "missing": 0}
                for species in all_species:
                    dna_value = dna.values.get(species, {}).get(sample, 0.0)
                    rna_value = rna.values.get(species, {}).get(sample, 0.0)
                    if not (
                        passes_threshold(dna_value, args.a)
                        and passes_threshold(rna_value, args.b)
                    ):
                        continue

                    status = detection_status(dna_value, rna_value)
                    status_counts[status] += 1
                    if dna_value > 0:
                        ratio = rna_value / dna_value
                        log2_ratio = math.log2(ratio) if ratio > 0 else None
                    else:
                        ratio = None
                        log2_ratio = None

                    writer.writerow(
                        [
                            sample,
                            species,
                            format_number(dna_value),
                            format_number(rna_value),
                            format_number(ratio),
                            format_number(log2_ratio),
                            status,
                        ]
                    )
                    retained += 1
                    output_rows += 1

                summary_rows.append(
                    [
                        sample,
                        len(all_species),
                        retained,
                        status_counts["both"],
                        status_counts["DNA_only"],
                        status_counts["RNA_only"],
                        status_counts["missing"],
                    ]
                )

        with (output_dir / "sample_summary.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "Sample",
                    "candidate_species",
                    "retained_species",
                    "both",
                    "DNA_only",
                    "RNA_only",
                    "missing",
                ]
            )
            writer.writerows(summary_rows)

        warnings: list[str] = []
        for label, table in (("DNA", dna), ("RNA", rna)):
            for sample, total in table.column_sums.items():
                if not math.isclose(total, 1.0, rel_tol=1e-6, abs_tol=1e-6):
                    warnings.append(
                        f"{label} sample {sample}: relative-abundance column sum is "
                        f"{total:.12g}, not approximately 1."
                    )

        with (output_dir / "run.log").open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(f"shared_samples\t{len(shared_samples)}\n")
            handle.write(f"candidate_species\t{len(all_species)}\n")
            handle.write(f"output_rows\t{output_rows}\n")
            handle.write(f"dna_threshold\t{args.a}\n")
            handle.write(f"rna_threshold\t{args.b}\n")
            for warning in warnings:
                handle.write(f"WARNING\t{warning}\n")

        config = {
            "script": "RNADNARatio.py",
            "script_version": SCRIPT_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "dna_table": str(dna_path),
            "rna_table": str(rna_path),
            "output_directory": str(output_dir),
            "dna_minimum_relative_abundance": args.a,
            "rna_minimum_relative_abundance": args.b,
            "shared_samples": shared_samples,
            "species_matching": (
                "exact after trimming whitespace and replacing internal whitespace "
                "with underscores; no taxonomic synonym conversion"
            ),
            "filter_rule": "DNA >= a AND RNA >= b; a or b equal to 0 disables that threshold",
            "ratio_rule": "RNA relative abundance / DNA relative abundance",
            "zero_rule": (
                "DNA=0: ratio and log2 are NA; DNA>0 and RNA=0: ratio=0 and log2=NA; "
                "no pseudocount is added"
            ),
        }
        with (output_dir / "config.json").open("w", encoding="utf-8") as handle:
            json.dump(config, handle, ensure_ascii=False, indent=2)
            handle.write("\n")

        for warning in warnings:
            print(f"WARNING: {warning}", file=sys.stderr)
        print(
            f"[RNADNARatio] Wrote {output_rows} row(s) for "
            f"{len(shared_samples)} shared sample(s) to {ratio_path}",
            file=sys.stderr,
        )
        return 0
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

