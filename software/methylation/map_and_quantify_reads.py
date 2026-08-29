#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import json
import math
import re
import shutil
import statistics
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Iterator, List, Tuple

PRESETS = {
    "MspJI": {"aliases": {"1", "mspji"}, "legacy": False, "contexts": {"CpG", "CHG"}},
    "AspBHI": {"aliases": {"2", "aspbhi"}, "legacy": False, "contexts": {"CpG", "CHG"}},
    "RlaI": {"aliases": {"3", "rlai"}, "legacy": False, "contexts": {"CHG"}},
    "SgrTI": {"aliases": {"4", "sgrti"}, "legacy": False, "contexts": {"CpG", "CHG"}},
    "FspEI": {"aliases": {"5", "fspei"}, "legacy": True, "contexts": {"CpG", "CHG"}},
}

CONTEXT_ALIASES = {
    "CpG": {"1", "cpg", "cg"},
    "CHG": {"2", "chg"},
}

STRUCTURE_KEYS = [
    "context",
    "preset_name",
    "motif",
    "motif_start_1based",
    "left_flank",
    "right_flank",
    "read_length",
    "both_orientations",
]


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


def load_samples(samples_path: Path) -> List[str]:
    with samples_path.open("rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("samples.tsv has no header.")
        if "sample_id" not in reader.fieldnames:
            raise ValueError("samples.tsv is missing required column: sample_id")

        sample_ids = []
        seen = set()
        valid_id = re.compile(r"^[A-Za-z0-9._-]+$")
        for row_number, row in enumerate(reader, 2):
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                raise ValueError(f"Empty sample_id at line {row_number}.")
            if not valid_id.match(sample_id):
                raise ValueError(
                    f"Invalid sample_id '{sample_id}' at line {row_number}. Use letters, numbers, '.', '_' or '-'."
                )
            if sample_id in seen:
                raise ValueError(f"Duplicate sample_id '{sample_id}'.")
            seen.add(sample_id)
            sample_ids.append(sample_id)

    if not sample_ids:
        raise ValueError("samples.tsv contains no samples.")
    return sample_ids


def load_json(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"Required configuration file not found: {path}")
    with path.open("rt") as handle:
        return json.load(handle)


def normalize_legacy_config(config: dict) -> dict:
    """Treat the original single-preset MspJI config as MspJI_CpG."""
    config = dict(config)
    if (
        config.get("mode") == "preset"
        and config.get("enzyme") == "MspJI"
        and config.get("motif") == "YNCGNR"
    ):
        config.setdefault("context", "CpG")
        config.setdefault("preset_name", "MspJI_CpG")
    return config


def compare_structure_configs(reference_cfg: dict, clean_cfg: dict) -> None:
    differences = []
    for key in STRUCTURE_KEYS:
        if reference_cfg.get(key) != clean_cfg.get(key):
            differences.append(
                f"{key}: reference={reference_cfg.get(key)!r}, clean={clean_cfg.get(key)!r}"
            )
    if differences:
        raise ValueError(
            "Reference-read and clean-read configurations do not match:\n  "
            + "\n  ".join(differences)
        )


def index_complete(prefix: Path) -> bool:
    small = [
        Path(str(prefix) + suffix)
        for suffix in [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"]
    ]
    large = [
        Path(str(prefix) + suffix)
        for suffix in [".1.ebwtl", ".2.ebwtl", ".3.ebwtl", ".4.ebwtl", ".rev.1.ebwtl", ".rev.2.ebwtl"]
    ]
    return all(path.is_file() for path in small) or all(path.is_file() for path in large)


def ensure_bowtie_index(reference_fasta: Path, reference_dir: Path, log_path: Path) -> Path:
    bowtie_build = shutil.which("bowtie-build")
    if not bowtie_build:
        raise RuntimeError("Required program not found in PATH: bowtie-build")

    # Internal index: created automatically on the first mapping run and reused later.
    index_dir = reference_dir / ".bowtie_index"
    index_dir.mkdir(parents=True, exist_ok=True)
    prefix = index_dir / "reference_reads"

    if index_complete(prefix):
        with log_path.open("wt") as log:
            log.write(f"Bowtie index already exists and was reused: {prefix}\n")
        return prefix

    for path in index_dir.glob("reference_reads*.ebwt*"):
        path.unlink()

    command = [bowtie_build, str(reference_fasta), str(prefix)]
    with log_path.open("wt") as log:
        log.write("COMMAND: " + " ".join(command) + "\n\n")
        log.flush()
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)

    if result.returncode != 0 or not index_complete(prefix):
        raise RuntimeError(f"bowtie-build failed; see {log_path}")
    return prefix


def parse_bowtie_log(log_path: Path) -> dict:
    text = log_path.read_text(errors="replace")
    patterns = {
        "processed": r"# reads processed:\s*([0-9,]+)",
        "at_least_one": r"# reads with at least one alignment:\s*([0-9,]+)",
        "failed": r"# reads that failed to align:\s*([0-9,]+)",
        "suppressed": r"# reads with alignments suppressed due to -m:\s*([0-9,]+)",
        "reported": r"Reported\s+([0-9,]+)\s+alignments",
    }
    values = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if match:
            values[key] = int(match.group(1).replace(",", ""))

    if "processed" not in values:
        raise RuntimeError(f"Could not parse Bowtie summary from {log_path}.")

    processed = values["processed"]
    failed = values.get("failed")
    at_least_one = values.get("at_least_one")
    suppressed = values.get("suppressed", 0)
    reported = values.get("reported")

    if at_least_one is None and failed is not None:
        at_least_one = processed - failed
    if failed is None and at_least_one is not None:
        failed = processed - at_least_one
    if reported is None and at_least_one is not None:
        reported = at_least_one - suppressed

    if at_least_one is None or failed is None or reported is None:
        raise RuntimeError(f"Incomplete Bowtie summary in {log_path}.")

    return {
        "processed": processed,
        "at_least_one": at_least_one,
        "failed": failed,
        "suppressed": suppressed,
        "reported": reported,
    }


def site_number(site_id: str) -> int:
    match = re.fullmatch(r"site(\d+)", site_id)
    if not match:
        raise ValueError(f"Unexpected reference site ID: {site_id}")
    return int(match.group(1))


def write_counts_file(counts: Counter[str], path: Path) -> None:
    with gzip.open(path, "wt") as out:
        out.write("site_id\traw_count\n")
        for site_id in sorted(counts, key=site_number):
            out.write(f"{site_id}\t{counts[site_id]}\n")


def calculate_depthmax(counts: Counter[str], top_fraction: float,
                       method: str) -> Tuple[int, float]:
    depths = list(counts.values())
    if not depths:
        return 0, float("nan")
    top_n = max(1, math.ceil(len(depths) * top_fraction))
    top_depths = heapq.nlargest(top_n, depths)
    if method == "mean":
        value = statistics.fmean(top_depths)
    else:
        value = float(statistics.median(top_depths))
    return top_n, float(value)


def map_sample(sample_id: str, clean_fastq: Path, index_prefix: Path,
               output_dir: Path, work_dir: Path, threads: int,
               mismatches: int, minimum_count: int, top_fraction: float,
               depthmax_method: str, keep_intermediate: bool) -> dict:
    bowtie = shutil.which("bowtie")
    if not bowtie:
        raise RuntimeError("Required program not found in PATH: bowtie")
    gzip_program = shutil.which("gzip")
    if clean_fastq.suffix == ".gz" and not gzip_program:
        raise RuntimeError("Required program not found in PATH: gzip")

    log_path = output_dir / "logs" / f"{sample_id}.bowtie.log"
    counts_path = work_dir / f"{sample_id}.counts.tsv.gz"
    unique_ids_path = output_dir / "intermediate" / f"{sample_id}.unique_site_ids.txt.gz"

    # Positional index syntax is compatible with Bowtie 1.x releases.
    command = [
        bowtie, "-q", "-v", str(mismatches), "-m", "1", "--best", "--strata",
        "--suppress", "1,2,4,5,6,7,8", "-p", str(threads), str(index_prefix), "-"
    ]

    counts: Counter[str] = Counter()
    with log_path.open("wb") as log:
        log.write(("COMMAND: " + " ".join(command) + "\n\n").encode())
        log.flush()
        decompressor = None
        input_handle = None
        try:
            if clean_fastq.suffix == ".gz":
                decompressor = subprocess.Popen(
                    [gzip_program, "-cd", str(clean_fastq)],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                bowtie_stdin = decompressor.stdout
            else:
                input_handle = clean_fastq.open("rb")
                bowtie_stdin = input_handle

            bowtie_process = subprocess.Popen(
                command,
                stdin=bowtie_stdin,
                stdout=subprocess.PIPE,
                stderr=log,
            )

            if decompressor is not None and decompressor.stdout is not None:
                decompressor.stdout.close()

            unique_out = gzip.open(unique_ids_path, "wt") if keep_intermediate else None
            try:
                assert bowtie_process.stdout is not None
                for raw_line in bowtie_process.stdout:
                    site_id = raw_line.decode("utf-8", errors="strict").strip()
                    if not site_id:
                        continue
                    site_number(site_id)
                    counts[site_id] += 1
                    if unique_out is not None:
                        unique_out.write(site_id + "\n")
            finally:
                if unique_out is not None:
                    unique_out.close()

            bowtie_return = bowtie_process.wait()
            decompressor_return = 0
            decompressor_error = b""
            if decompressor is not None:
                decompressor_return = decompressor.wait()
                if decompressor.stderr is not None:
                    decompressor_error = decompressor.stderr.read()

            # Check Bowtie first: a Bowtie failure can make gzip receive SIGPIPE.
            if bowtie_return != 0:
                raise RuntimeError(f"Bowtie failed for {sample_id}; see {log_path}")
            if decompressor_return != 0:
                message = decompressor_error.decode(errors="replace").strip()
                if message:
                    raise RuntimeError(f"FASTQ decompression failed for {sample_id}: {message}")
                raise RuntimeError(f"FASTQ decompression failed for {sample_id}.")
        finally:
            if input_handle is not None:
                input_handle.close()
            if decompressor is not None and decompressor.poll() is None:
                decompressor.kill()
                decompressor.wait()

    mapping = parse_bowtie_log(log_path)
    unique_reads = sum(counts.values())
    if mapping["reported"] != unique_reads:
        raise RuntimeError(
            f"Bowtie/count mismatch for {sample_id}: log reports {mapping['reported']} unique alignments, "
            f"but {unique_reads} were counted."
        )

    write_counts_file(counts, counts_path)
    detected_sites = len(counts)
    qualified_sites = sum(count >= minimum_count for count in counts.values())
    top_site_number, depthmax = calculate_depthmax(counts, top_fraction, depthmax_method)

    processed = mapping["processed"]
    if processed > 0:
        total_mapping_fraction = mapping["at_least_one"] / processed
        unique_mapping_fraction = unique_reads / processed
        multi_mapping_fraction = mapping["suppressed"] / processed
    else:
        total_mapping_fraction = 0.0
        unique_mapping_fraction = 0.0
        multi_mapping_fraction = 0.0

    return {
        "sample_id": sample_id,
        "counts_path": str(counts_path),
        "clean_reads": processed,
        "tags_with_at_least_one_alignment": mapping["at_least_one"],
        "failed_to_align": mapping["failed"],
        "multi_mapped_suppressed": mapping["suppressed"],
        "uniquely_mapped_reads": unique_reads,
        "total_mapping_fraction": total_mapping_fraction,
        "unique_mapping_fraction": unique_mapping_fraction,
        "multi_mapping_fraction": multi_mapping_fraction,
        "detected_sites": detected_sites,
        "qualified_sites": qualified_sites,
        "top_site_number": top_site_number,
        "depthmax": depthmax,
    }


def iter_counts(path: Path) -> Iterator[Tuple[int, str, int]]:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header != ["site_id", "raw_count"]:
            raise ValueError(f"Unexpected counts header in {path}: {header}")
        for line in handle:
            if not line.strip():
                continue
            site_id, count_text = line.rstrip("\n").split("\t")
            yield site_number(site_id), site_id, int(count_text)


def union_count_rows(count_paths: List[Path]) -> Iterator[Tuple[str, List[int]]]:
    iterators = [iter_counts(path) for path in count_paths]
    heap = []
    for sample_index, iterator in enumerate(iterators):
        try:
            number, site_id, count = next(iterator)
            heapq.heappush(heap, (number, sample_index, site_id, count))
        except StopIteration:
            pass

    while heap:
        number = heap[0][0]
        row_counts = [0] * len(iterators)
        site_id_for_row = None
        while heap and heap[0][0] == number:
            _, sample_index, site_id, count = heapq.heappop(heap)
            if site_id_for_row is None:
                site_id_for_row = site_id
            elif site_id_for_row != site_id:
                raise RuntimeError(
                    f"Inconsistent site IDs for numeric site {number}: {site_id_for_row} vs {site_id}"
                )
            row_counts[sample_index] = count
            try:
                next_number, next_site_id, next_count = next(iterators[sample_index])
                heapq.heappush(heap, (next_number, sample_index, next_site_id, next_count))
            except StopIteration:
                pass
        if site_id_for_row is None:
            raise RuntimeError("Internal matrix merge error.")
        yield site_id_for_row, row_counts


def iter_reference_sites(sites_path: Path) -> Iterator[Tuple[int, dict]]:
    with gzip.open(sites_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"site_id", "seqname", "start_1based", "strand", "sequence", "copy_number"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected sites.tsv.gz header in {sites_path}.")
        for row in reader:
            yield site_number(row["site_id"]), row


def format_metric(value: float) -> str:
    if isinstance(value, float) and math.isnan(value):
        return "NA"
    return f"{value:.8f}"


def rpm_value(count: int, unique_reads: int, minimum_count: int) -> str:
    if count < minimum_count or unique_reads <= 0:
        return "NA"
    return f"{count / unique_reads * 1_000_000:.8f}"


def m_index_value(count: int, depthmax: float, minimum_count: int) -> str:
    if count < minimum_count or not math.isfinite(depthmax) or depthmax <= 1:
        return "NA"
    return f"{math.log(count) / math.log(depthmax):.8f}"


def write_final_outputs(sample_ids: List[str], sample_stats: List[dict],
                        reference_sites_path: Path, output_dir: Path,
                        minimum_count: int) -> None:
    count_paths = [Path(stats["counts_path"]) for stats in sample_stats]
    unique_reads = [stats["uniquely_mapped_reads"] for stats in sample_stats]
    depthmax_values = [stats["depthmax"] for stats in sample_stats]

    samples_dir = output_dir / "samples"
    samples_dir.mkdir(parents=True, exist_ok=True)
    sample_handles = []
    for sample_id in sample_ids:
        handle = gzip.open(samples_dir / f"{sample_id}.sites.tsv.gz", "wt")
        handle.write("site_id\tseqname\tstart_1based\tsequence\traw_count\trpm\tm_index\n")
        sample_handles.append(handle)

    count_out = gzip.open(output_dir / "count_matrix.tsv.gz", "wt")
    rpm_out = gzip.open(output_dir / "rpm_matrix.tsv.gz", "wt")
    mindex_out = gzip.open(output_dir / "m_index_matrix.tsv.gz", "wt")

    count_out.write("site_id\t" + "\t".join(sample_ids) + "\n")
    matrix_prefix = "site_id\tseqname\tstart_1based\tsequence"
    rpm_out.write(matrix_prefix + "\t" + "\t".join(sample_ids) + "\n")
    mindex_out.write(matrix_prefix + "\t" + "\t".join(sample_ids) + "\n")

    reference_iterator = iter_reference_sites(reference_sites_path)
    try:
        reference_number, reference_row = next(reference_iterator)
    except StopIteration:
        reference_number, reference_row = None, None

    try:
        for site_id, counts in union_count_rows(count_paths):
            number = site_number(site_id)
            while reference_number is not None and reference_number < number:
                try:
                    reference_number, reference_row = next(reference_iterator)
                except StopIteration:
                    reference_number, reference_row = None, None
                    break

            if reference_number != number or reference_row is None:
                raise RuntimeError(f"Detected site {site_id} was not found in reference sites.tsv.gz.")
            if reference_row["site_id"] != site_id:
                raise RuntimeError(
                    f"Reference site ID mismatch: expected {site_id}, found {reference_row['site_id']}."
                )

            seqname = reference_row["seqname"]
            start = reference_row["start_1based"]
            sequence = reference_row["sequence"]

            count_out.write(site_id + "\t" + "\t".join(str(count) for count in counts) + "\n")

            rpms = [
                rpm_value(counts[i], unique_reads[i], minimum_count)
                for i in range(len(sample_ids))
            ]
            m_indexes = [
                m_index_value(counts[i], depthmax_values[i], minimum_count)
                for i in range(len(sample_ids))
            ]

            for i, count in enumerate(counts):
                if count <= 0:
                    continue
                sample_handles[i].write(
                    f"{site_id}\t{seqname}\t{start}\t{sequence}\t{count}\t{rpms[i]}\t{m_indexes[i]}\n"
                )

            if any(count >= minimum_count for count in counts):
                prefix = f"{site_id}\t{seqname}\t{start}\t{sequence}"
                rpm_out.write(prefix + "\t" + "\t".join(rpms) + "\n")
                mindex_out.write(prefix + "\t" + "\t".join(m_indexes) + "\n")
    finally:
        count_out.close()
        rpm_out.close()
        mindex_out.close()
        for handle in sample_handles:
            handle.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map clean reads to reference reads and generate raw-count, RPM and M-index matrices."
    )
    parser.add_argument("-i", "--samples", required=True,
                        help="The same samples.tsv used for read extraction.")
    parser.add_argument("-e", "--enzyme", default=None,
                        help="Enzyme preset: 1=MspJI, 2=AspBHI, 3=RlaI, 4=SgrTI, 5=FspEI (legacy). Default: MspJI.")
    parser.add_argument("-t", "--context", default=None,
                        help="Methylation-context preset: 1=CpG, 2=CHG. Default: CpG when supported.")
    parser.add_argument("-o", "--output-dir", default=None,
                        help="Final output directory. Default: results/<enzyme>_<context>.")
    parser.add_argument("-p", "--threads", type=int, default=4,
                        help="Bowtie threads per sample. Default: 4.")
    parser.add_argument("-v", "--mismatches", type=int, default=2,
                        help="Maximum mismatches per read. Default: 2.")
    parser.add_argument("-x", "--minimum-count", type=int, default=3,
                        help="Minimum raw count required to calculate RPM/M-index. Default: 3.")
    parser.add_argument("-f", "--top-fraction", type=float, default=0.02,
                        help="Top fraction of detected sites used to calculate depthmax. Default: 0.02.")
    parser.add_argument("-d", "--depthmax-method", choices=["mean", "median"], default="mean",
                        help="How to summarize top-depth sites. Default: mean.")
    parser.add_argument("-r", "--reference-dir", default=None,
                        help="Reference-read directory. Default: reference_reads/<enzyme>_<context>.")
    parser.add_argument("-c", "--clean-dir", default=None,
                        help="Clean-read directory. Default: clean_reads/<enzyme>_<context>.")
    parser.add_argument("-k", "--keep-intermediate", action="store_true",
                        help="Keep per-read uniquely mapped site IDs.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        if args.threads < 1:
            raise ValueError("-p/--threads must be >= 1.")
        if args.mismatches < 0:
            raise ValueError("-v/--mismatches must be >= 0.")
        if args.minimum_count < 1:
            raise ValueError("-x/--minimum-count must be >= 1.")
        if not (0 < args.top_fraction <= 1):
            raise ValueError("-f/--top-fraction must be in (0, 1].")
        if shutil.which("bowtie") is None:
            raise RuntimeError("Required program not found in PATH: bowtie")
        if shutil.which("bowtie-build") is None:
            raise RuntimeError("Required program not found in PATH: bowtie-build")

        samples_path = Path(args.samples).expanduser().resolve()
        if not samples_path.is_file():
            raise FileNotFoundError(f"samples.tsv not found: {samples_path}")
        sample_ids = load_samples(samples_path)

        explicit_preset = args.enzyme is not None or args.context is not None
        enzyme = normalize_enzyme(args.enzyme)
        context = normalize_context(args.context, enzyme)
        selected_preset = preset_name(enzyme, context)
        if PRESETS[enzyme].get("legacy"):
            print("WARNING: FspEI is retained as a legacy preset for historical/existing datasets.", file=sys.stderr)
        reference_dir = (
            Path(args.reference_dir).expanduser().resolve()
            if args.reference_dir else (Path("reference_reads") / selected_preset).resolve()
        )
        clean_dir = (
            Path(args.clean_dir).expanduser().resolve()
            if args.clean_dir else (Path("clean_reads") / selected_preset).resolve()
        )

        if not reference_dir.is_dir():
            raise FileNotFoundError(f"Reference-read directory not found: {reference_dir}")
        if not clean_dir.is_dir():
            raise FileNotFoundError(f"Clean-read directory not found: {clean_dir}")

        reference_fasta = reference_dir / "reference_reads.fa"
        reference_sites = reference_dir / "sites.tsv.gz"
        if not reference_fasta.is_file():
            raise FileNotFoundError(f"Reference reads not found: {reference_fasta}")
        if not reference_sites.is_file():
            raise FileNotFoundError(f"Reference sites table not found: {reference_sites}")

        reference_cfg = normalize_legacy_config(load_json(reference_dir / "config.json"))
        clean_cfg = normalize_legacy_config(load_json(clean_dir / "config.json"))
        compare_structure_configs(reference_cfg, clean_cfg)

        if explicit_preset:
            if (
                reference_cfg.get("mode") != "preset"
                or reference_cfg.get("enzyme") != enzyme
                or reference_cfg.get("context") != context
            ):
                raise ValueError(
                    f"-e {enzyme} -t {context} does not match the reference-read configuration."
                )
            if (
                clean_cfg.get("mode") != "preset"
                or clean_cfg.get("enzyme") != enzyme
                or clean_cfg.get("context") != context
            ):
                raise ValueError(
                    f"-e {enzyme} -t {context} does not match the clean-read configuration."
                )

        configuration_name = (
            reference_cfg.get("preset_name") or preset_name(
                reference_cfg.get("enzyme") or "MspJI",
                reference_cfg.get("context") or "CpG",
            )
            if reference_cfg.get("mode") == "preset" else "custom"
        )
        output_dir = Path(args.output_dir).expanduser() if args.output_dir else Path("results") / configuration_name
        output_dir.mkdir(parents=True, exist_ok=True)

        for filename in [
            "count_matrix.tsv.gz", "rpm_matrix.tsv.gz", "m_index_matrix.tsv.gz",
            "sample_qc.tsv", "config.json"
        ]:
            path = output_dir / filename
            if path.exists():
                path.unlink()

        samples_output_dir = output_dir / "samples"
        if samples_output_dir.exists():
            shutil.rmtree(samples_output_dir)
        intermediate_dir = output_dir / "intermediate"
        if intermediate_dir.exists():
            shutil.rmtree(intermediate_dir)
        logs_dir = output_dir / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        if args.keep_intermediate:
            intermediate_dir.mkdir(parents=True, exist_ok=True)
        work_dir = output_dir / ".work"
        if work_dir.exists():
            shutil.rmtree(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        index_prefix = ensure_bowtie_index(
            reference_fasta, reference_dir, logs_dir / "bowtie_build.log"
        )

        sample_stats = []
        try:
            for sample_id in sample_ids:
                clean_fastq = clean_dir / f"{sample_id}.fastq.gz"
                if not clean_fastq.is_file():
                    alternative = clean_dir / f"{sample_id}.fastq"
                    if alternative.is_file():
                        clean_fastq = alternative
                    else:
                        raise FileNotFoundError(
                            f"Clean FASTQ not found for sample '{sample_id}': {clean_fastq}"
                        )

                print(f"Mapping {sample_id}", flush=True)
                stats = map_sample(
                    sample_id=sample_id,
                    clean_fastq=clean_fastq,
                    index_prefix=index_prefix,
                    output_dir=output_dir,
                    work_dir=work_dir,
                    threads=args.threads,
                    mismatches=args.mismatches,
                    minimum_count=args.minimum_count,
                    top_fraction=args.top_fraction,
                    depthmax_method=args.depthmax_method,
                    keep_intermediate=args.keep_intermediate,
                )
                sample_stats.append(stats)

            write_final_outputs(
                sample_ids, sample_stats, reference_sites, output_dir, args.minimum_count
            )

            qc_fields = [
                "sample_id", "clean_reads", "total_mapping_fraction",
                "unique_mapping_fraction", "multi_mapping_fraction",
                "detected_sites", "qualified_sites", "depthmax"
            ]
            with (output_dir / "sample_qc.tsv").open("wt", newline="") as out:
                writer = csv.DictWriter(out, fieldnames=qc_fields, delimiter="\t")
                writer.writeheader()
                for stats in sample_stats:
                    writer.writerow({
                        "sample_id": stats["sample_id"],
                        "clean_reads": stats["clean_reads"],
                        "total_mapping_fraction": f"{stats['total_mapping_fraction']:.8f}",
                        "unique_mapping_fraction": f"{stats['unique_mapping_fraction']:.8f}",
                        "multi_mapping_fraction": f"{stats['multi_mapping_fraction']:.8f}",
                        "detected_sites": stats["detected_sites"],
                        "qualified_sites": stats["qualified_sites"],
                        "depthmax": format_metric(stats["depthmax"]),
                    })

            config = {
                "enzyme": reference_cfg.get("enzyme") or "custom",
                "context": reference_cfg.get("context"),
                "preset_name": reference_cfg.get("preset_name") or configuration_name,
                "reference_dir": str(reference_dir),
                "clean_dir": str(clean_dir),
                "samples": str(samples_path),
                "threads": args.threads,
                "mismatches": args.mismatches,
                "minimum_count": args.minimum_count,
                "top_fraction": args.top_fraction,
                "depthmax_method": args.depthmax_method,
                "keep_intermediate": bool(args.keep_intermediate),
                "output_dir": str(output_dir.resolve()),
            }
            with (output_dir / "config.json").open("wt") as out:
                json.dump(config, out, indent=2, ensure_ascii=False)
                out.write("\n")
        finally:
            if work_dir.exists():
                shutil.rmtree(work_dir)

        print(f"Output directory : {output_dir}")
        print(f"Samples processed: {len(sample_ids)}")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
