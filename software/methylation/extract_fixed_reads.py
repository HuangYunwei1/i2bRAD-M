#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
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


def validate_motif(motif: str) -> str:
    motif = motif.upper()
    if not motif:
        raise ValueError("Motif cannot be empty.")
    for char in motif:
        if char not in IUPAC:
            raise ValueError(f"Unsupported IUPAC code '{char}' in motif '{motif}'.")
    return motif


def motif_matches(seq: str, motif: str) -> bool:
    if len(seq) != len(motif):
        return False
    for base, code in zip(seq.upper(), motif.upper()):
        allowed = IUPAC.get(code)
        if allowed is None or base not in allowed:
            return False
    return True


def open_fastq_text(path: Path, mode: str):
    return gzip.open(path, mode) if str(path).endswith(".gz") else path.open(mode)


def read_fastq(path: Path) -> Iterator[Tuple[str, str, str, str]]:
    with open_fastq_text(path, "rt") as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            record_number += 1
            if not (seq and plus and qual):
                raise ValueError(f"Incomplete FASTQ record {record_number} in {path}.")
            header = header.rstrip("\r\n")
            seq = seq.rstrip("\r\n").upper()
            plus = plus.rstrip("\r\n")
            qual = qual.rstrip("\r\n")
            if not header.startswith("@"):
                raise ValueError(f"FASTQ record {record_number} in {path} does not start with '@'.")
            if not plus.startswith("+"):
                raise ValueError(f"FASTQ record {record_number} in {path} has an invalid '+' line.")
            if len(seq) != len(qual):
                raise ValueError(f"Sequence/quality length mismatch in FASTQ record {record_number} of {path}.")
            yield header, seq, plus, qual


def load_samples(samples_path: Path) -> List[dict]:
    base_dir = samples_path.parent
    with samples_path.open("rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("samples.tsv has no header.")
        required = {"sample_id", "fastq"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError("samples.tsv is missing required column(s): " + ", ".join(sorted(missing)))

        samples = []
        seen = set()
        valid_id = re.compile(r"^[A-Za-z0-9._-]+$")
        for row_number, row in enumerate(reader, 2):
            sample_id = (row.get("sample_id") or "").strip()
            fastq_text = (row.get("fastq") or "").strip()
            if not sample_id:
                raise ValueError(f"Empty sample_id at line {row_number}.")
            if not valid_id.match(sample_id):
                raise ValueError(
                    f"Invalid sample_id '{sample_id}' at line {row_number}. Use letters, numbers, '.', '_' or '-'."
                )
            if sample_id in seen:
                raise ValueError(f"Duplicate sample_id '{sample_id}'.")
            seen.add(sample_id)
            if not fastq_text:
                raise ValueError(f"Empty FASTQ path for sample '{sample_id}'.")
            fastq = Path(fastq_text).expanduser()
            fastq = (base_dir / fastq).resolve() if not fastq.is_absolute() else fastq.resolve()
            if not fastq.is_file():
                raise FileNotFoundError(f"FASTQ not found for sample '{sample_id}': {fastq}")
            samples.append({"sample_id": sample_id, "fastq": str(fastq)})

    if not samples:
        raise ValueError("samples.tsv contains no samples.")
    return samples


def candidate_offset_groups(max_shift: int) -> List[List[int]]:
    groups = [[0]]
    for distance in range(1, max_shift + 1):
        groups.append([-distance, distance])
    return groups


def extract_candidate(seq: str, qual: str, structure: dict,
                      adapter: str | None) -> Tuple[str, str, str]:
    motif = structure["motif"]
    motif_len = len(motif)
    expected_start0 = structure["input_motif_start_1based"] - 1
    left_flank = structure["left_flank"]
    right_flank = structure["right_flank"]
    max_shift = structure["max_shift"]
    both = structure["both_orientations"]
    reverse_motif = reverse_complement(motif)
    check_reverse = both and reverse_motif != motif

    adapter_start = None
    if adapter:
        adapter_index = seq.find(adapter)
        if adapter_index >= 0:
            adapter_start = adapter_index

    # Priority: expected position -> ±1 -> ±2 -> ...
    for offsets in candidate_offset_groups(max_shift):
        hits = []
        for offset in offsets:
            motif_start = expected_start0 + offset
            motif_end = motif_start + motif_len
            if motif_start < 0 or motif_end > len(seq):
                continue
            observed = seq[motif_start:motif_end]
            orientations = []
            if motif_matches(observed, motif):
                orientations.append("+")
            if check_reverse and motif_matches(observed, reverse_motif):
                orientations.append("-")

            orientation_hits = []
            for orientation in orientations:
                window_start = motif_start - left_flank
                window_end = motif_end + right_flank
                if window_start < 0 or window_end > len(seq):
                    continue
                if adapter_start is not None and window_end > adapter_start:
                    continue

                extracted_seq = seq[window_start:window_end]
                extracted_qual = qual[window_start:window_end]
                if orientation == "-":
                    extracted_seq = reverse_complement(extracted_seq)
                    extracted_qual = extracted_qual[::-1]
                orientation_hits.append(
                    (motif_start, orientation, extracted_seq, extracted_qual)
                )

            if len(orientation_hits) == 1:
                hits.append(orientation_hits[0])
            elif len(orientation_hits) > 1:
                # Degenerate non-palindromic motif and its reverse-complement
                # can both match the same position. This is one physical
                # fragment, not an ambiguous positional hit. Normalize to the
                # same deterministic orientation used by the reference builder.
                hits.append(min(orientation_hits, key=lambda item: (item[2], item[1])))

        if len(hits) == 1:
            _, _, output_seq, output_qual = hits[0]
            return "matched", output_seq, output_qual
        if len(hits) > 1:
            return "ambiguous", "", ""

    return "no_match", "", ""


def process_sample(task: dict) -> dict:
    sample_id = task["sample_id"]
    fastq = Path(task["fastq"])
    output_path = Path(task["output_path"])
    structure = task["structure"]
    quality_threshold = task["quality_threshold"]
    max_low_quality_bases = task["max_low_quality_bases"]
    adapter = task["adapter"]

    counts = {
        "sample_id": sample_id,
        "input_reads": 0,
        "motif_matched": 0,
        "ambiguous_removed": 0,
        "low_quality_removed": 0,
        "output_reads": 0,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        output_path.unlink()
    temp_path = output_path.with_name(output_path.name + f".tmp.{os.getpid()}")

    try:
        with gzip.open(temp_path, "wt", compresslevel=6) as out:
            for header, seq, plus, qual in read_fastq(fastq):
                counts["input_reads"] += 1
                status, output_seq, output_qual = extract_candidate(
                    seq, qual, structure, adapter
                )
                if status == "no_match":
                    continue
                if status == "ambiguous":
                    counts["ambiguous_removed"] += 1
                    continue

                counts["motif_matched"] += 1
                if any(base not in "ACGT" for base in output_seq):
                    counts["ambiguous_removed"] += 1
                    continue

                low_quality = sum((ord(char) - 33) < quality_threshold for char in output_qual)
                if low_quality > max_low_quality_bases:
                    counts["low_quality_removed"] += 1
                    continue

                out.write(f"{header}\n{output_seq}\n{plus}\n{output_qual}\n")
                counts["output_reads"] += 1

        temp_path.replace(output_path)
    except Exception:
        try:
            temp_path.unlink()
        except FileNotFoundError:
            pass
        raise

    counts["retention_fraction"] = (
        counts["output_reads"] / counts["input_reads"] if counts["input_reads"] else 0.0
    )
    return counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch-extract fixed-length reads from single-end or R1 FASTQ files."
    )
    parser.add_argument("-i", "--samples", required=True,
                        help="Two-column samples.tsv with sample_id and fastq.")
    parser.add_argument("-e", "--enzyme", default=None,
                        help="Enzyme preset: 1=MspJI, 2=AspBHI, 3=RlaI, 4=SgrTI, 5=FspEI. Default: MspJI.")
    parser.add_argument("-t", "--context", default=None,
                        help="Methylation-context preset: 1=CpG, 2=CHG. Default: CpG when supported.")
    parser.add_argument("-o", "--output-dir", default=None,
                        help="Output directory. Default: clean_reads/<enzyme>_<context>.")
    parser.add_argument("-p", "--parallel-samples", type=int, default=1,
                        help="Number of samples processed in parallel. Default: 1.")
    parser.add_argument("-q", "--quality-threshold", type=int, default=10,
                        help="Phred threshold defining a low-quality base. Default: 10.")
    parser.add_argument("-n", "--max-low-quality-bases", type=int, default=5,
                        help="Maximum number of bases below -q allowed per extracted read. Default: 5.")

    custom = parser.add_argument_group("custom structure")
    custom.add_argument("-u", "--custom", action="store_true",
                        help="Enable custom structure mode. Cannot be combined with -e.")
    custom.add_argument("-c", "--motif", default=None, help="Custom IUPAC motif.")
    custom.add_argument("-s", "--motif-start", type=int, default=None,
                        help="Expected motif start in the input read, 1-based.")
    custom.add_argument("-L", "--left-flank", type=int, default=None,
                        help="Bases retained to the left of the motif.")
    custom.add_argument("-R", "--right-flank", type=int, default=None,
                        help="Bases retained to the right of the motif.")
    custom.add_argument("-m", "--max-shift", type=int, default=0,
                        help="Maximum motif shift around the expected position. Default: 0.")
    custom.add_argument("-a", "--adapter", default=None,
                        help="Optional exact adapter sequence. The extracted read may not extend into it.")
    custom.add_argument("-b", "--both-orientations", action="store_true",
                        help="Accept reverse-complement structure and normalize to reference-read orientation.")
    return parser.parse_args()


def resolve_structure(args: argparse.Namespace) -> dict:
    if args.custom:
        if args.enzyme is not None or args.context is not None:
            raise ValueError("-u/--custom cannot be combined with -e/--enzyme or -t/--context.")
        missing = [
            flag for flag, value in [
                ("-c", args.motif), ("-s", args.motif_start),
                ("-L", args.left_flank), ("-R", args.right_flank)
            ] if value is None
        ]
        if missing:
            raise ValueError("Custom mode requires " + ", ".join(missing) + ".")
        motif = validate_motif(args.motif)
        if args.motif_start < 1:
            raise ValueError("-s/--motif-start must be >= 1.")
        if args.left_flank < 0 or args.right_flank < 0:
            raise ValueError("-L and -R must be >= 0.")
        if args.max_shift < 0:
            raise ValueError("-m/--max-shift must be >= 0.")
        return {
            "mode": "custom",
            "enzyme": None,
            "context": None,
            "preset_name": "custom",
            "motif": motif,
            "input_motif_start_1based": args.motif_start,
            "motif_start_1based": args.left_flank + 1,
            "left_flank": args.left_flank,
            "right_flank": args.right_flank,
            "read_length": args.left_flank + len(motif) + args.right_flank,
            "max_shift": args.max_shift,
            "both_orientations": bool(args.both_orientations),
        }

    forbidden = []
    for flag, value in [
        ("-c", args.motif), ("-s", args.motif_start),
        ("-L", args.left_flank), ("-R", args.right_flank)
    ]:
        if value is not None:
            forbidden.append(flag)
    if args.max_shift != 0:
        forbidden.append("-m")
    if args.adapter is not None:
        forbidden.append("-a")
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
        "input_motif_start_1based": cfg["motif_start_1based"],
        "motif_start_1based": cfg["left_flank"] + 1,
        "left_flank": cfg["left_flank"],
        "right_flank": cfg["right_flank"],
        "read_length": cfg["left_flank"] + len(motif) + cfg["right_flank"],
        "max_shift": cfg["max_shift"],
        "both_orientations": cfg["both_orientations"],
    }


def main() -> None:
    args = parse_args()
    try:
        if args.parallel_samples < 1:
            raise ValueError("-p/--parallel-samples must be >= 1.")
        if args.quality_threshold < 0:
            raise ValueError("-q/--quality-threshold must be >= 0.")
        if args.max_low_quality_bases < 0:
            raise ValueError("-n/--max-low-quality-bases must be >= 0.")

        structure = resolve_structure(args)

        samples_path = Path(args.samples).expanduser().resolve()
        if not samples_path.is_file():
            raise FileNotFoundError(f"samples.tsv not found: {samples_path}")
        samples = load_samples(samples_path)

        default_name = structure["preset_name"] if structure["mode"] == "preset" else "custom"
        output_dir = Path(args.output_dir).expanduser() if args.output_dir else Path("clean_reads") / default_name
        output_dir.mkdir(parents=True, exist_ok=True)

        for path in output_dir.glob("*.fastq.gz"):
            path.unlink()
        for filename in ["extraction_qc.tsv", "config.json"]:
            path = output_dir / filename
            if path.exists():
                path.unlink()

        adapter = args.adapter.upper() if args.adapter else None
        if adapter and any(base not in "ACGT" for base in adapter):
            raise ValueError("-a/--adapter currently accepts an exact A/C/G/T sequence.")

        tasks = []
        for sample in samples:
            tasks.append({
                **sample,
                "output_path": str(output_dir / f"{sample['sample_id']}.fastq.gz"),
                "structure": structure,
                "quality_threshold": args.quality_threshold,
                "max_low_quality_bases": args.max_low_quality_bases,
                "adapter": adapter,
            })

        workers = min(args.parallel_samples, len(tasks))
        results_by_id = {}
        if workers == 1:
            for task in tasks:
                print(f"Processing {task['sample_id']}", flush=True)
                result = process_sample(task)
                results_by_id[result["sample_id"]] = result
        else:
            print(f"Processing {len(tasks)} samples with {workers} parallel workers.", flush=True)
            with ProcessPoolExecutor(max_workers=workers) as executor:
                future_to_id = {executor.submit(process_sample, task): task["sample_id"] for task in tasks}
                for future in as_completed(future_to_id):
                    sample_id = future_to_id[future]
                    result = future.result()
                    results_by_id[sample_id] = result
                    print(f"Finished {sample_id}", flush=True)

        ordered_results = [results_by_id[sample["sample_id"]] for sample in samples]
        fields = [
            "sample_id", "input_reads", "motif_matched", "ambiguous_removed",
            "low_quality_removed", "output_reads", "retention_fraction"
        ]
        with (output_dir / "extraction_qc.tsv").open("wt", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in ordered_results:
                row = dict(row)
                row["retention_fraction"] = f"{row['retention_fraction']:.6f}"
                writer.writerow({key: row[key] for key in fields})

        config = dict(structure)
        config.update({
            "samples": str(samples_path),
            "quality_threshold": args.quality_threshold,
            "max_low_quality_bases": args.max_low_quality_bases,
            "adapter": adapter,
            "parallel_samples": workers,
            "output_dir": str(output_dir.resolve()),
        })
        with (output_dir / "config.json").open("wt") as out:
            json.dump(config, out, indent=2, ensure_ascii=False)
            out.write("\n")

        print(f"Output directory : {output_dir}")
        print(f"Samples processed: {len(samples)}")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
