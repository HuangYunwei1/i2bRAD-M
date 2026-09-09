#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterator, List, Tuple

IUPAC: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}

COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)

PRESETS = {
    "MspJI": {
        "aliases": {"1", "mspji"},
        "contexts": {
            "CpG": {
                "motif": "YNCGNR",
                "motif_start_1based": 14,
                "left_flank": 13,
                "right_flank": 13,
                "max_shift": 0,
                "both_orientations": False,
            },
            "CHG": {
                "motif": "YCHGR",
                "motif_start_1based": 14,
                "left_flank": 13,
                "right_flank": 13,
                "max_shift": 0,
                "both_orientations": True,
            },
        },
    },
    "AspBHI": {
        "aliases": {"2", "aspbhi"},
        "contexts": {
            "CpG": {
                "motif": "YSCGSR",
                "motif_start_1based": 14,
                "left_flank": 13,
                "right_flank": 13,
                "max_shift": 0,
                "both_orientations": False,
            },
            "CHG": {
                "motif": "YNCHGNR",
                "motif_start_1based": 13,
                "left_flank": 12,
                "right_flank": 12,
                "max_shift": 1,
                "both_orientations": True,
            },
        },
    },
    "RlaI": {
        "aliases": {"3", "rlai"},
        "contexts": {
            "CHG": {
                "motif": "VCWGB",
                "motif_start_1based": 14,
                "left_flank": 13,
                "right_flank": 13,
                "max_shift": 1,
                "both_orientations": False,
            },
        },
    },
    "SgrTI": {
        "aliases": {"4", "sgrti"},
        "contexts": {
            "CpG": {
                "motif": "SCGS",
                "motif_start_1based": 15,
                "left_flank": 14,
                "right_flank": 14,
                "max_shift": 0,
                "both_orientations": False,
            },
            "CHG": {
                "motif": "BCWGV",
                "motif_start_1based": 14,
                "left_flank": 13,
                "right_flank": 13,
                "max_shift": 1,
                "both_orientations": False,
            },
        },
    },
    "FspEI": {
        "aliases": {"5", "fspei"},
        "contexts": {
            "CpG": {
                "motif": "CCGG",
                "motif_start_1based": 15,
                "left_flank": 14,
                "right_flank": 14,
                "max_shift": 0,
                "both_orientations": False,
            },
            "CHG": {
                "motif": "CCWGG",
                "motif_start_1based": 14,
                "left_flank": 13,
                "right_flank": 13,
                "max_shift": 1,
                "both_orientations": False,
            },
        },
    },
}

CONTEXT_ALIASES = {
    "CpG": {"1", "cpg", "cg"},
    "CHG": {"2", "chg"},
}


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def normalize_enzyme(value: str | None) -> str:
    if value is None:
        return "MspJI"
    text = value.strip()
    for name, cfg in PRESETS.items():
        if text == name or text.lower() in cfg["aliases"]:
            return name
    supported = ", ".join(f"{next(iter(sorted(cfg['aliases'], key=len)))}={name}" for name, cfg in PRESETS.items())
    raise ValueError(
        f"Unsupported enzyme preset '{value}'. Supported enzymes: {supported}."
    )


def normalize_context(value: str | None, enzyme: str) -> str:
    supported = PRESETS[enzyme]["contexts"]
    if value is None:
        if "CpG" in supported:
            return "CpG"
        return next(iter(supported))
    text = value.strip()
    canonical = None
    for name, aliases in CONTEXT_ALIASES.items():
        if text == name or text.lower() in aliases:
            canonical = name
            break
    if canonical is None or canonical not in supported:
        choices = ", ".join(supported)
        raise ValueError(
            f"Unsupported methylation context '{value}' for {enzyme}. Supported context(s): {choices}."
        )
    return canonical


def preset_name(enzyme: str, context: str) -> str:
    return f"{enzyme}_{context}"


def preset_config(enzyme: str, context: str) -> dict:
    return PRESETS[enzyme]["contexts"][context]


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open("rt")


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    name = None
    chunks: List[str] = []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}.")
                name = header.split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(
                        f"FASTA sequence encountered before the first header at line {line_number}."
                    )
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def motif_to_regex(motif: str) -> re.Pattern[str]:
    parts = []
    for char in motif.upper():
        if char not in IUPAC:
            raise ValueError(f"Unsupported IUPAC code '{char}' in motif '{motif}'.")
        bases = IUPAC[char]
        parts.append(bases if len(bases) == 1 else f"[{bases}]")
    # Look-ahead preserves overlapping motif occurrences.
    return re.compile("(?=(" + "".join(parts) + "))")


def encode_acgt(seq: str) -> int | None:
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    code = 0
    for base in seq:
        value = mapping.get(base)
        if value is None:
            return None
        code = (code << 2) | value
    return code


def find_sites_in_sequence(
    seq: str,
    motif: str,
    left_flank: int,
    right_flank: int,
    both_orientations: bool,
) -> Iterator[Tuple[int, str, str]]:
    """Yield (window_start_0based, strand, oriented_reference_read)."""
    motif = motif.upper()
    motif_len = len(motif)
    orientations = [("+", motif_to_regex(motif))]
    reverse_motif = reverse_complement(motif)
    if both_orientations and reverse_motif != motif:
        orientations.append(("-", motif_to_regex(reverse_motif)))

    # A non-palindromic degenerate motif and its reverse-complement can
    # occasionally match the same genomic position. Treat that as one
    # physical site and choose a deterministic canonical orientation,
    # rather than emitting the same site twice.
    hits_by_start: Dict[int, List[Tuple[int, str, str]]] = {}
    for strand, regex in orientations:
        for match in regex.finditer(seq):
            motif_start = match.start()
            window_start = motif_start - left_flank
            window_end = motif_start + motif_len + right_flank
            if window_start < 0 or window_end > len(seq):
                continue
            read_seq = seq[window_start:window_end]
            if strand == "-":
                read_seq = reverse_complement(read_seq)
            hits_by_start.setdefault(motif_start, []).append(
                (window_start, strand, read_seq)
            )

    for motif_start in sorted(hits_by_start):
        options = hits_by_start[motif_start]
        if len(options) == 1:
            yield options[0]
            continue
        # If both orientations match the same locus, normalize both to a
        # single canonical sequence. This keeps reference construction
        # consistent with extract_fixed_reads.py.
        window_start, strand, read_seq = min(
            options, key=lambda item: (item[2], item[1])
        )
        yield window_start, strand, read_seq


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a fixed-length reference-read library from a reference genome."
    )
    parser.add_argument("-r", "--reference", required=True, help="Reference FASTA or FASTA.GZ.")
    parser.add_argument(
        "-e", "--enzyme", default=None,
        help="Enzyme preset: 1=MspJI, 2=AspBHI, 3=RlaI, 4=SgrTI, 5=FspEI. Default: MspJI.",
    )
    parser.add_argument(
        "-t", "--context", default=None,
        help="Methylation-context preset: 1=CpG, 2=CHG. Default: CpG when supported.",
    )
    parser.add_argument(
        "-o", "--output-dir", default=None,
        help="Output directory. Default: reference_reads/<enzyme>_<context>.",
    )

    custom = parser.add_argument_group("custom structure")
    custom.add_argument(
        "-u", "--custom", action="store_true",
        help="Enable custom structure mode. Cannot be combined with -e.",
    )
    custom.add_argument("-c", "--motif", default=None, help="Custom IUPAC motif.")
    custom.add_argument("-L", "--left-flank", type=int, default=None,
                        help="Bases retained to the left of the motif.")
    custom.add_argument("-R", "--right-flank", type=int, default=None,
                        help="Bases retained to the right of the motif.")
    custom.add_argument(
        "-b", "--both-orientations", action="store_true",
        help="For a non-palindromic motif, also scan its reverse-complement orientation.",
    )
    return parser.parse_args()


def resolve_structure(args: argparse.Namespace) -> dict:
    if args.custom:
        if args.enzyme is not None or args.context is not None:
            raise ValueError("-u/--custom cannot be combined with -e/--enzyme or -t/--context.")
        missing = [
            flag for flag, value in [
                ("-c", args.motif), ("-L", args.left_flank), ("-R", args.right_flank)
            ] if value is None
        ]
        if missing:
            raise ValueError("Custom mode requires " + ", ".join(missing) + ".")
        motif = args.motif.upper()
        motif_to_regex(motif)
        if args.left_flank < 0 or args.right_flank < 0:
            raise ValueError("-L and -R must be >= 0.")
        return {
            "mode": "custom",
            "enzyme": None,
            "context": None,
            "preset_name": "custom",
            "motif": motif,
            "motif_start_1based": args.left_flank + 1,
            "left_flank": args.left_flank,
            "right_flank": args.right_flank,
            "read_length": args.left_flank + len(motif) + args.right_flank,
            "both_orientations": bool(args.both_orientations),
        }

    forbidden = []
    for flag, value in [("-c", args.motif), ("-L", args.left_flank), ("-R", args.right_flank)]:
        if value is not None:
            forbidden.append(flag)
    if args.both_orientations:
        forbidden.append("-b")
    if forbidden:
        raise ValueError("Custom structure parameters require -u/--custom: " + ", ".join(forbidden))

    enzyme = normalize_enzyme(args.enzyme)
    context = normalize_context(args.context, enzyme)
    cfg = preset_config(enzyme, context)
    motif = cfg["motif"]
    return {
        "mode": "preset",
        "enzyme": enzyme,
        "context": context,
        "preset_name": preset_name(enzyme, context),
        "motif": motif,
        "motif_start_1based": cfg["left_flank"] + 1,
        "left_flank": cfg["left_flank"],
        "right_flank": cfg["right_flank"],
        "read_length": cfg["left_flank"] + len(motif) + cfg["right_flank"],
        "both_orientations": cfg["both_orientations"],
    }


def first_pass(reference: Path, structure: dict) -> Tuple[Counter, int, int]:
    sequence_counts: Counter[int] = Counter()
    total_sites = 0
    ambiguous_sites = 0
    for _, seq in read_fasta(reference):
        for _, _, read_seq in find_sites_in_sequence(
            seq, structure["motif"], structure["left_flank"],
            structure["right_flank"], structure["both_orientations"]
        ):
            total_sites += 1
            encoded = encode_acgt(read_seq)
            if encoded is None:
                ambiguous_sites += 1
            else:
                sequence_counts[encoded] += 1
    return sequence_counts, total_sites, ambiguous_sites


def second_pass(reference: Path, output_dir: Path, structure: dict,
                sequence_counts: Counter) -> None:
    sites_path = output_dir / "sites.tsv.gz"
    reads_path = output_dir / "reference_reads.fa"
    site_number = 0

    with gzip.open(sites_path, "wt") as sites_out, reads_path.open("wt") as reads_out:
        sites_out.write(
            "site_id\tseqname\tstart_1based\tstrand\tsequence\tcopy_number\n"
        )
        for seqname, seq in read_fasta(reference):
            for window_start, strand, read_seq in find_sites_in_sequence(
                seq, structure["motif"], structure["left_flank"],
                structure["right_flank"], structure["both_orientations"]
            ):
                site_number += 1
                site_id = f"site{site_number:09d}"
                encoded = encode_acgt(read_seq)
                copy_number = "NA" if encoded is None else str(sequence_counts[encoded])
                sites_out.write(
                    f"{site_id}\t{seqname}\t{window_start + 1}\t{strand}\t{read_seq}\t{copy_number}\n"
                )
                if encoded is not None:
                    reads_out.write(f">{site_id}\n{read_seq}\n")


def main() -> None:
    args = parse_args()
    try:
        structure = resolve_structure(args)
        reference = Path(args.reference).expanduser().resolve()
        if not reference.is_file():
            raise FileNotFoundError(f"Reference FASTA not found: {reference}")

        default_name = structure["preset_name"] if structure["mode"] == "preset" else "custom"
        output_dir = Path(args.output_dir).expanduser() if args.output_dir else Path("reference_reads") / default_name
        output_dir.mkdir(parents=True, exist_ok=True)

        for filename in ["reference_reads.fa", "sites.tsv.gz", "summary.tsv", "config.json"]:
            path = output_dir / filename
            if path.exists():
                path.unlink()
        index_dir = output_dir / ".bowtie_index"
        if index_dir.exists():
            shutil.rmtree(index_dir)


        sequence_counts, total_sites, ambiguous_sites = first_pass(reference, structure)
        usable_sites = sum(sequence_counts.values())
        unique_sites = sum(1 for count in sequence_counts.values() if count == 1)
        repeated_sites = sum(count for count in sequence_counts.values() if count > 1)
        if usable_sites + ambiguous_sites != total_sites:
            raise RuntimeError("Internal site-count consistency check failed.")

        second_pass(reference, output_dir, structure, sequence_counts)

        with (output_dir / "summary.tsv").open("wt") as out:
            out.write("metric\tvalue\n")
            out.write(f"reference\t{reference}\n")
            out.write(f"enzyme\t{structure['enzyme'] or 'custom'}\n")
            out.write(f"context\t{structure.get('context') or 'custom'}\n")
            out.write(f"preset_name\t{structure.get('preset_name') or 'custom'}\n")
            out.write(f"motif\t{structure['motif']}\n")
            out.write(f"read_length\t{structure['read_length']}\n")
            out.write(f"total_theoretical_sites\t{total_sites}\n")
            out.write(f"unique_site_records\t{unique_sites}\n")
            out.write(f"repeated_site_records\t{repeated_sites}\n")
            out.write(f"ambiguous_site_records\t{ambiguous_sites}\n")

        config = dict(structure)
        config["reference"] = str(reference)
        config["output_dir"] = str(output_dir.resolve())
        with (output_dir / "config.json").open("wt") as out:
            json.dump(config, out, indent=2, ensure_ascii=False)
            out.write("\n")

        print(f"Output directory        : {output_dir}")
        print(f"Total theoretical sites: {total_sites}")
        print(f"Unique site records     : {unique_sites}")
        print(f"Repeated site records   : {repeated_sites}")
        print(f"Ambiguous site records  : {ambiguous_sites}")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
