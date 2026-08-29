#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MAP2BAbsoluteQuantifier.py v0.2.0

Absolute concentration:
    Ci = CH * RAi / RAH

Optional library cell equivalent (-c):
    Ni = Ci * VE * VL / VD

Metadata without -c:
    sample_id  CH

Metadata with -c:
    sample_id  CH  VE  VD  VL

Units:
    CH: cells/mL
    VE: mL
    VD and VL: same unit, usually uL
"""

import argparse
import csv
import math
import sys
from pathlib import Path

VERSION = "0.2.0"


def error(message):
    raise SystemExit("ERROR: " + message)


def number_text(value):
    if value == 0:
        return "0"
    return f"{value:.10g}"


def read_abundance(filename):
    path = Path(filename)
    if not path.is_file():
        error(f"input file not found: {path}")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = [
            row for row in csv.reader(handle, delimiter="\t")
            if row and any(x.strip() for x in row)
        ]

    if len(rows) < 2:
        error("Abundance.xls has no data rows")

    header = rows[0]
    data = rows[1:]
    names = [x.lstrip("#").strip().lower() for x in header]

    if "species" not in names:
        error("cannot find the Species column")

    species_col = names.index("species")
    sample_names = [x.strip() for x in header[species_col + 1:]]

    if not sample_names:
        error("no sample columns were found after Species")
    if any(not x for x in sample_names):
        error("an empty sample name was found")
    if len(sample_names) != len(set(sample_names)):
        error("duplicate sample names were found")

    for line_no, row in enumerate(data, 2):
        if len(row) != len(header):
            error(
                f"line {line_no} has {len(row)} columns; "
                f"expected {len(header)}"
            )
        for col, sample in enumerate(sample_names, species_col + 1):
            try:
                value = float(row[col])
            except ValueError:
                error(
                    f"non-numeric abundance at line {line_no}, "
                    f"sample {sample}: {row[col]!r}"
                )
            if not math.isfinite(value) or value < 0:
                error(
                    f"invalid abundance at line {line_no}, "
                    f"sample {sample}: {row[col]!r}"
                )

    return header, data, species_col, sample_names


def read_metadata(filename, calculate_cells):
    path = Path(filename)
    if not path.is_file():
        error(f"metadata file not found: {path}")

    lines = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            lines.append((line_no, line.split()))

    if len(lines) < 2:
        error("metadata must contain a header and at least one sample")

    header = [x.lower() for x in lines[0][1]]
    required = ["sample_id", "ch"]
    if calculate_cells:
        required += ["ve", "vd", "vl"]

    for name in required:
        if name not in header:
            error(f"metadata is missing column: {name}")

    columns = {name: header.index(name) for name in required}
    metadata = {}

    for line_no, fields in lines[1:]:
        if len(fields) != len(header):
            error(
                f"metadata line {line_no} has {len(fields)} columns; "
                f"expected {len(header)}"
            )

        sample = fields[columns["sample_id"]]
        if sample in metadata:
            error(f"duplicate sample_id in metadata: {sample}")

        values = {}
        for name in required[1:]:
            raw = fields[columns[name]]
            try:
                value = float(raw)
            except ValueError:
                error(f"non-numeric {name} for sample {sample}: {raw!r}")
            if not math.isfinite(value) or value <= 0:
                error(f"{name} must be greater than 0 for sample {sample}")
            values[name] = value

        metadata[sample] = values

    return metadata


def find_host(data, species_col, host):
    matches = [
        index for index, row in enumerate(data)
        if row[species_col].strip() == host
    ]

    if not matches:
        error(f"host species not found: {host}")
    if len(matches) > 1:
        error(f"host species matched more than one row: {host}")

    return matches[0]


def write_table(filename, header, rows):
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Convert MAP2B relative abundance to absolute concentration."
    )
    parser.add_argument("-i", required=True, help="MAP2B Abundance.xls")
    parser.add_argument("-m", required=True, help="sample metadata file")
    parser.add_argument(
        "-s", required=True,
        help="host name exactly as written in the Species column"
    )
    parser.add_argument(
        "-o", default="AbsoluteQuantification",
        help="output directory [default: AbsoluteQuantification]"
    )
    parser.add_argument(
        "-c", action="store_true",
        help="also calculate library cell equivalents using VE, VD and VL"
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {VERSION}"
    )
    args = parser.parse_args()

    header, data, species_col, samples = read_abundance(args.i)
    metadata = read_metadata(args.m, args.c)

    missing = [sample for sample in samples if sample not in metadata]
    if missing:
        error("metadata is missing samples: " + ", ".join(missing))

    extra = [sample for sample in metadata if sample not in samples]
    if extra:
        print(
            "WARNING: metadata contains unused samples: " + ", ".join(extra),
            file=sys.stderr,
        )

    host_row = find_host(data, species_col, args.s)
    host_ra = {}

    for col, sample in enumerate(samples, species_col + 1):
        host_ra[sample] = float(data[host_row][col])
        if host_ra[sample] <= 0:
            error(f"host abundance is zero for sample: {sample}")

    concentration_rows = []
    cell_rows = []

    for row in data:
        concentration_row = row[:species_col + 1]
        cell_row = row[:species_col + 1]

        for col, sample in enumerate(samples, species_col + 1):
            rai = float(row[col])
            ch = metadata[sample]["ch"]
            ci = ch * rai / host_ra[sample]
            concentration_row.append(number_text(ci))

            if args.c:
                ve = metadata[sample]["ve"]
                vd = metadata[sample]["vd"]
                vl = metadata[sample]["vl"]
                ni = ci * ve * vl / vd
                cell_row.append(number_text(ni))

        concentration_rows.append(concentration_row)
        if args.c:
            cell_rows.append(cell_row)

    output_dir = Path(args.o)

    if output_dir.exists() and not output_dir.is_dir():
        error(
            f"output path exists but is not a directory: {output_dir}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    concentration_output = output_dir / "AbsoluteConcentration.xls"
    write_table(concentration_output, header, concentration_rows)
    print(f"Absolute concentration table: {concentration_output}")

    if args.c:
        cell_output = output_dir / "LibraryCellEquivalent.xls"
        write_table(cell_output, header, cell_rows)
        print(f"Library cell-equivalent table: {cell_output}")


if __name__ == "__main__":
    main()
