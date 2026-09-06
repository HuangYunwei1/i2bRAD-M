#!/usr/bin/env python3
"""Shared I/O, validation and enzyme utilities for Serial2BRAD16."""
from __future__ import annotations

import csv
import gzip
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Sequence, TextIO, Tuple


@dataclass(frozen=True)
class FastqRecord:
    header: str
    seq: str
    plus: str
    qual: str

    @property
    def read_id(self) -> str:
        return canonical_read_id(self.header)


def open_text(path: os.PathLike | str, mode: str = "rt") -> TextIO:
    path = str(path)
    if "b" in mode:
        raise ValueError("open_text only accepts text modes")
    if path.lower().endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def canonical_read_id(header: str) -> str:
    value = header.strip().lstrip("@>").split()[0]
    return re.sub(r"/(1|2)$", "", value)


def iter_fastq(path: os.PathLike | str, strict: bool = True) -> Iterator[FastqRecord]:
    with open_text(path, "rt") as handle:
        record_no = 0
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            record_no += 1
            if not (seq and plus and qual):
                raise ValueError(f"FASTQ truncated at record {record_no}: {path}")
            record = FastqRecord(*(x.rstrip("\r\n") for x in (header, seq, plus, qual)))
            if strict:
                if not record.header.startswith("@"):
                    raise ValueError(f"FASTQ header does not start with @ at record {record_no}: {path}")
                if not record.plus.startswith("+"):
                    raise ValueError(f"FASTQ separator does not start with + at record {record_no}: {path}")
                if len(record.seq) != len(record.qual):
                    raise ValueError(f"Sequence/quality length mismatch at record {record_no}: {path}")
            yield record


def write_fastq_record(handle: TextIO, rec: FastqRecord) -> None:
    handle.write(f"{rec.header}\n{rec.seq}\n{rec.plus}\n{rec.qual}\n")


def iter_fasta(path: os.PathLike | str) -> Iterator[Tuple[str, str]]:
    with open_text(path, "rt") as handle:
        header: Optional[str] = None
        parts: List[str] = []
        for line_no, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line, []
            else:
                if header is None:
                    raise ValueError(f"FASTA sequence before header at {path}:{line_no}")
                parts.append(line)
        if header is not None:
            yield header, "".join(parts)


def write_fasta_record(handle: TextIO, header: str, seq: str) -> None:
    handle.write(f">{header.lstrip('@>')}\n{seq}\n")


def strip_suffix(name: str, suffixes: Sequence[str]) -> str:
    low = name.lower()
    for suffix in suffixes:
        if low.endswith(suffix):
            return name[:-len(suffix)]
    return name


def strip_fastq_suffix(name: str) -> str:
    return strip_suffix(name, (".fastq.gz", ".fq.gz", ".fastq", ".fq"))


def strip_fasta_suffix(name: str) -> str:
    return strip_suffix(name, (".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna"))


def safe_filename(value: str) -> str:
    value = value.strip().replace("/", "_").replace("\\", "_")
    value = re.sub(r"[\x00-\x1f<>:\"|?*]+", "_", value)
    return value or "unnamed"


def quality_passes(qual: str, min_q: int, min_percent: float, phred_offset: int) -> bool:
    if not qual:
        return False
    passing = sum((ord(char) - phred_offset) > min_q for char in qual)
    return passing * 100.0 / len(qual) >= min_percent


def default_config_path() -> Path:
    return Path(__file__).resolve().parents[1] / "config" / "enzymes.json"


def load_config(path: os.PathLike | str | None = None) -> dict:
    cfg_path = Path(path) if path else default_config_path()
    with open(cfg_path, "r", encoding="utf-8-sig") as handle:
        config = json.load(handle)
    validate_config(config)
    return config


def validate_config(config: dict) -> None:
    enzymes = config.get("enzymes")
    if not isinstance(enzymes, list) or len(enzymes) != 16:
        raise ValueError("Enzyme configuration must contain exactly 16 enzymes")
    ids = {int(e["id"]) for e in enzymes}
    if ids != set(range(1, 17)):
        raise ValueError("Enzyme IDs must be exactly 1-16")
    for e in enzymes:
        if len(e.get("windows", [])) != 5:
            raise ValueError(f"{e['name']}: exactly five windows are required")
        for pattern in e.get("patterns", []):
            compiled = re.compile(pattern, re.IGNORECASE)
            probe = regex_probe(pattern)
            match = compiled.fullmatch(probe)
            if not match or len(probe) != int(e["tag_length"]):
                raise ValueError(f"{e['name']}: invalid recognition pattern/tag length")
        expected = 5 * int(e["serial_unit_length"]) + 4 * int(config["serial_junction_length"])
        if expected != int(e["expected_assembled_length"]):
            raise ValueError(f"{e['name']}: expected-length derivation is inconsistent")
        margin = 5 * int(config["broad_tolerance_per_tag"])
        if (expected - margin, expected + margin) != (int(e["pear_broad_min"]), int(e["pear_broad_max"])):
            raise ValueError(f"{e['name']}: broad-range derivation is inconsistent")


def regex_probe(pattern: str) -> str:
    value = re.sub(r"\[ACGT\]\{(\d+)\}", lambda m: "A" * int(m.group(1)), pattern)
    value = re.sub(r"\[([ACGT]+)\]", lambda m: m.group(1)[0], value)
    if re.search(r"[\[\]{}()]", value):
        raise ValueError(f"Unsupported pattern syntax: {pattern}")
    return value


def get_enzyme(enzyme_id: int | str, config: dict) -> dict:
    raw = str(enzyme_id).strip()
    if not raw.isdigit() or int(raw) not in range(1, 17):
        raise ValueError("Enzyme must be a numeric preset ID from 1 to 16")
    for enzyme in config["enzymes"]:
        if int(enzyme["id"]) == int(raw):
            return enzyme
    raise ValueError(f"Enzyme ID not found: {raw}")


def resolve_path(value: str, base: os.PathLike | str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = Path(base) / path
    return path.resolve()


def read_manifest(path: os.PathLike | str) -> List[dict]:
    path = Path(path)
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"Manifest is empty: {path}")
    reader = csv.DictReader(lines, delimiter="\t")
    required = ["library_id", "r1", "r2", "t1", "t2", "t3", "t4", "t5"]
    missing = [column for column in required if column not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"Manifest missing columns: {', '.join(missing)}")
    if "enzyme" in (reader.fieldnames or []):
        raise ValueError("Manifest must not contain an enzyme column; select one enzyme globally with -e 1-16")
    rows, seen = [], set()
    for line_no, row in enumerate(reader, 2):
        cleaned = {column: (row.get(column) or "").strip() for column in required}
        library_id = cleaned["library_id"]
        if not library_id:
            raise ValueError(f"Manifest line {line_no}: empty library_id")
        if library_id in seen:
            raise ValueError(f"Manifest line {line_no}: duplicate library_id {library_id}")
        seen.add(library_id)
        if not cleaned["r1"] or not cleaned["r2"]:
            raise ValueError(f"Manifest line {line_no}: r1 and r2 are required")
        for pos in ("t1", "t2", "t3", "t4", "t5"):
            if not cleaned[pos]:
                raise ValueError(f"Manifest line {line_no}: {pos} sample ID is empty")
        cleaned["r1_path"] = resolve_path(cleaned["r1"], path.parent)
        cleaned["r2_path"] = resolve_path(cleaned["r2"], path.parent)
        rows.append(cleaned)
    return rows


def find_pattern_hits(seq: str, patterns: Sequence[re.Pattern]) -> List[Tuple[int, int, str, int]]:
    seq = seq.upper()
    found, seen = [], set()
    for pattern_no, pattern in enumerate(patterns, 1):
        for start in range(len(seq)):
            match = pattern.match(seq, start)
            if not match:
                continue
            key = (match.start(), match.end(), match.group(0).upper())
            if key not in seen:
                seen.add(key)
                found.append((match.start(), match.end(), match.group(0).upper(), pattern_no))
    return sorted(found, key=lambda item: (item[0], item[1], item[3]))


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
FASTA_SUFFIXES = (".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna")


def list_sequence_files(directory: os.PathLike | str, suffixes: Sequence[str]) -> List[Path]:
    """Return matching files below *directory* in deterministic order."""
    root = Path(directory)
    if not root.is_dir():
        raise ValueError(f"Input directory does not exist: {root}")
    suffixes_low = tuple(value.lower() for value in suffixes)
    return sorted(
        (path for path in root.rglob("*") if path.is_file() and path.name.lower().endswith(suffixes_low)),
        key=lambda path: str(path).casefold(),
    )


def infer_fastq_identity(path: os.PathLike | str, stage: str) -> Optional[Tuple[str, int]]:
    """Infer ``(library_prefix, mate)`` from a workflow FASTQ filename.

    ``stage='goodquality'`` accepts GoodQuality files and excludes PoorQuality/filter files.
    ``stage='filtered'`` accepts the paired outputs consumed by PEAR.
    """
    path = Path(path)
    core = strip_fastq_suffix(path.name)
    low = core.lower()
    if stage == "goodquality":
        if "goodquality" not in low or "poorquality" in low or "filter" in low:
            return None
        core = re.sub(r"(?i)(?:\.ATCG)?\.GoodQuality$", "", core)
        core = strip_fastq_suffix(core)
    elif stage == "filtered":
        if "filter" not in low:
            return None
        patterns = (
            r"(?i)(?:\.ATCG)?\.GoodQuality(?:\.(?:fastq|fq))?(?:-filter|\.filter)$",
            r"(?i)(?:\.ATCG)?\.GoodQuality(?:-filter|\.filter)$",
            r"(?i)(?:-filter|\.filter)$",
        )
        original = core
        for pattern in patterns:
            core = re.sub(pattern, "", core)
            if core != original:
                break
        core = strip_fastq_suffix(core)
    else:
        raise ValueError(f"Unsupported pairing stage: {stage}")

    matches = list(re.finditer(r"(?i)(?:^|[._-])(R[12])(?=$|[._-])", core))
    if len(matches) != 1:
        raise ValueError(f"Cannot identify exactly one R1/R2 token in filename: {path.name}")
    match = matches[0]
    mate = int(match.group(1)[1])
    prefix = (core[:match.start()] + core[match.end():]).strip("._-")
    if not prefix:
        raise ValueError(f"Cannot derive a library prefix from filename: {path.name}")
    return safe_filename(prefix), mate


def discover_fastq_pairs(directory: os.PathLike | str, stage: str) -> List[Tuple[str, Path, Path]]:
    """Discover and validate all R1/R2 pairs for a workflow directory."""
    grouped: dict[str, dict[int, Path | str]] = {}
    display_prefix: dict[str, str] = {}
    candidates = list_sequence_files(directory, FASTQ_SUFFIXES)
    for path in candidates:
        identity = infer_fastq_identity(path, stage)
        if identity is None:
            continue
        prefix, mate = identity
        key = prefix.casefold()
        display_prefix.setdefault(key, prefix)
        slot = grouped.setdefault(key, {})
        if mate in slot:
            raise ValueError(
                f"Duplicate R{mate} files for library {display_prefix[key]}: {slot[mate]} ; {path}"
            )
        slot[mate] = path
    if not grouped:
        label = "GoodQuality" if stage == "goodquality" else "filtered paired"
        raise ValueError(f"No {label} FASTQ files were found under: {directory}")
    missing = []
    result = []
    for key in sorted(grouped, key=lambda value: display_prefix[value].casefold()):
        slot = grouped[key]
        if 1 not in slot or 2 not in slot:
            absent = "R1" if 1 not in slot else "R2"
            missing.append(f"{display_prefix[key]} missing {absent}")
            continue
        result.append((display_prefix[key], Path(slot[1]), Path(slot[2])))
    if missing:
        raise ValueError("Unpaired libraries detected: " + "; ".join(missing))
    return result


def list_assembled_fastq(directory: os.PathLike | str) -> List[Path]:
    files = list_sequence_files(directory, FASTQ_SUFFIXES)
    return [path for path in files if strip_fastq_suffix(path.name).lower().endswith(".assembled")]


def list_assembled_fasta(directory: os.PathLike | str) -> List[Path]:
    files = list_sequence_files(directory, FASTA_SUFFIXES)
    return [path for path in files if strip_fasta_suffix(path.name).lower().endswith(".assembled")]


def list_split_fasta(directory: os.PathLike | str) -> List[Path]:
    files = list_sequence_files(directory, FASTA_SUFFIXES)
    return [path for path in files if re.search(r"(?i)_T[1-5]$", strip_fasta_suffix(path.name))]
