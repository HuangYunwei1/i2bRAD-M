#!/usr/bin/env python3
"""Build, update, or validate a MAP2B-compatible species-specific 2b-tag database.

MAP2BDatabaseBuilder v0.3.2

Modes:
  -m 1  build a new database
  -m 2  update an existing database
  -m 3  validate an existing database

The builder supports database construction, update, and validation for 16
Type IIB restriction enzymes and writes databases compatible with the i2bRAD-M
MAP2B workflow.
"""

from __future__ import annotations

import argparse
import collections
import dataclasses
import datetime as dt
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, Iterable, Iterator, List, Mapping, Sequence, Tuple

try:
    import marisa_trie  # type: ignore
except ImportError:
    marisa_trie = None


__version__ = "0.3.2"
DEFAULT_ENZYME_ID = 13
DEFAULT_SHARD_SIZE = 30000
DEFAULT_PROCESSES = 1
ALL_ENZYME_ID = 17
TAXONOMY_FILENAME = "abfh_classify_with_speciename.txt.gz"
ID_MAP_FILENAME = "genome_id_map.tsv"
NEW_ID_MAP_FILENAME = "new_genome_id_map.tsv"


@dataclasses.dataclass(frozen=True)
class EnzymeSpec:
    enzyme_id: int
    name: str
    length: int
    forward_pattern: str
    palindromic: bool = False

    @property
    def regex(self) -> re.Pattern[str]:
        return re.compile(rf"(?=({self.forward_pattern}))")


# The 16 Type IIB patterns follow the MAP2B/2bRAD-M database-building scripts.
# For palindromic patterns we scan one orientation only, preventing one physical
# site from being counted twice when the reverse-complement sequence is scanned.
ENZYMES: Tuple[EnzymeSpec, ...] = (
    EnzymeSpec(1, "CspCI", 33, r"[ACGT]{11}CAA[ACGT]{5}GTGG[ACGT]{10}"),
    EnzymeSpec(2, "AloI", 27, r"[ACGT]{7}GAAC[ACGT]{6}TCC[ACGT]{7}"),
    EnzymeSpec(3, "BsaXI", 27, r"[ACGT]{9}AC[ACGT]{5}CTCC[ACGT]{7}"),
    EnzymeSpec(4, "BaeI", 28, r"[ACGT]{10}AC[ACGT]{4}GTA[CT]C[ACGT]{7}"),
    EnzymeSpec(5, "BcgI", 32, r"[ACGT]{10}CGA[ACGT]{6}TGC[ACGT]{10}"),
    EnzymeSpec(6, "CjeI", 28, r"[ACGT]{8}CCA[ACGT]{6}GT[ACGT]{9}"),
    EnzymeSpec(7, "PpiI", 27, r"[ACGT]{7}GAAC[ACGT]{5}CTC[ACGT]{8}"),
    EnzymeSpec(8, "PsrI", 27, r"[ACGT]{7}GAAC[ACGT]{6}TAC[ACGT]{7}"),
    EnzymeSpec(9, "BplI", 27, r"[ACGT]{8}GAG[ACGT]{5}CTC[ACGT]{8}", True),
    EnzymeSpec(10, "FalI", 27, r"[ACGT]{8}AAG[ACGT]{5}CTT[ACGT]{8}", True),
    EnzymeSpec(11, "Bsp24I", 27, r"[ACGT]{8}GAC[ACGT]{6}TGG[ACGT]{7}"),
    EnzymeSpec(12, "HaeIV", 27, r"[ACGT]{7}GA[CT][ACGT]{5}[AG]TC[ACGT]{9}"),
    EnzymeSpec(13, "CjePI", 27, r"[ACGT]{7}CCA[ACGT]{7}TC[ACGT]{8}"),
    EnzymeSpec(14, "Hin4I", 27, r"[ACGT]{8}GA[CT][ACGT]{5}[ACG]TC[ACGT]{8}"),
    EnzymeSpec(15, "AlfI", 32, r"[ACGT]{10}GCA[ACGT]{6}TGC[ACGT]{10}", True),
    EnzymeSpec(16, "BslFI", 25, r"[ACGT]{6}GGGAC[ACGT]{14}"),
)
ENZYME_BY_ID = {enzyme.enzyme_id: enzyme for enzyme in ENZYMES}
ENZYME_BY_NAME = {enzyme.name.lower(): enzyme for enzyme in ENZYMES}
DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclasses.dataclass
class TaxonomyRecord:
    fields: List[str]
    source_id: str

    @property
    def internal_id(self) -> str:
        return self.fields[0]

    @internal_id.setter
    def internal_id(self, value: str) -> None:
        self.fields[0] = value

    @property
    def species(self) -> str:
        return self.fields[7]

    @property
    def taxon_signature(self) -> str:
        return "\t".join(self.fields[1:8])

    @property
    def genome_path(self) -> str:
        return self.fields[-1]


class BuilderError(RuntimeError):
    pass


def log(message: str) -> None:
    stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{stamp} - MAP2BDatabaseBuilder - {message}", flush=True)


def warn(message: str) -> None:
    stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{stamp} - MAP2BDatabaseBuilder - WARNING - {message}", file=sys.stderr, flush=True)


def require_marisa() -> None:
    if marisa_trie is None:
        raise BuilderError(
            "Python package 'marisa_trie' is unavailable. Activate the MAP2B conda environment first."
        )


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1].upper()


def canonical_tag(tag: str) -> str:
    tag = tag.upper()
    rc = reverse_complement(tag)
    return tag if tag <= rc else rc


def parse_enzyme_selector(value: str) -> int:
    raw = str(value).strip()
    if not raw:
        raise argparse.ArgumentTypeError("Empty enzyme selector")
    lowered = raw.lower()
    if lowered in {"all", "allenzyme", "allenzymes"}:
        return ALL_ENZYME_ID
    if raw.isdigit():
        enzyme_id = int(raw)
        if enzyme_id == ALL_ENZYME_ID or enzyme_id in ENZYME_BY_ID:
            return enzyme_id
    else:
        enzyme = ENZYME_BY_NAME.get(lowered)
        if enzyme is not None:
            return enzyme.enzyme_id
    choices = ", ".join(f"{e.enzyme_id}:{e.name}" for e in ENZYMES)
    raise argparse.ArgumentTypeError(
        f"Unknown enzyme {value!r}. Choose 1-16/name, or 17/AllEnzyme. Supported: {choices}"
    )


def enzyme_specs_from_selector(selector: int) -> List[EnzymeSpec]:
    if selector == ALL_ENZYME_ID:
        return list(ENZYMES)
    return [ENZYME_BY_ID[selector]]


def open_text(path: Path, mode: str = "rt"):
    if "r" in mode and not path.exists():
        raise BuilderError(f"File does not exist: {path}")
    if path.suffix.lower() == ".gz":
        return gzip.open(path, mode, encoding=None if "b" in mode else "utf-8")
    return open(path, mode, encoding=None if "b" in mode else "utf-8")


def read_fasta(path: Path) -> Iterator[str]:
    """Yield uppercase contig sequences from a plain or gzip-compressed FASTA."""
    seen_header = False
    with open_text(path, "rt") as handle:
        seq_parts: List[str] = []
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen_header and seq_parts:
                    yield "".join(seq_parts).upper()
                seen_header = True
                seq_parts = []
            else:
                if not seen_header:
                    raise BuilderError(f"Not a FASTA file (sequence before header): {path}:{line_number}")
                seq_parts.append(line)
        if seen_header and seq_parts:
            yield "".join(seq_parts).upper()
    if not seen_header:
        raise BuilderError(f"No FASTA header found in {path}")


def digest_genome(path: Path, enzyme: EnzymeSpec) -> Dict[str, int]:
    """Return distinct tag -> assembly copy number for one reference assembly."""
    pattern = enzyme.regex
    counts: Dict[str, int] = collections.defaultdict(int)
    for sequence in read_fasta(path):
        for match in pattern.finditer(sequence):
            tag = match.group(1)
            if len(tag) != enzyme.length:
                raise BuilderError(
                    f"Internal enzyme-length mismatch for {enzyme.name}: expected {enzyme.length}, got {len(tag)}"
                )
            normalized = canonical_tag(tag) if enzyme.palindromic else tag.upper()
            counts[normalized] += 1
        if not enzyme.palindromic:
            rc_sequence = reverse_complement(sequence)
            for match in pattern.finditer(rc_sequence):
                tag = match.group(1)
                if len(tag) != enzyme.length:
                    raise BuilderError(
                        f"Internal enzyme-length mismatch for {enzyme.name}: expected {enzyme.length}, got {len(tag)}"
                    )
                counts[tag.upper()] += 1
    return dict(counts)


def read_taxonomy(
    path: Path,
    validate_genomes: bool = True,
) -> Tuple[List[str] | None, List[TaxonomyRecord]]:
    header: List[str] | None = None
    records: List[TaxonomyRecord] = []
    with open_text(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n\r")
            if not line.strip():
                continue
            fields = line.split("\t")
            if line.startswith("#"):
                if header is None:
                    header = fields
                continue
            if len(fields) < 9:
                raise BuilderError(
                    f"Taxonomy line {line_number} has {len(fields)} columns; at least 9 are required: {path}"
                )
            if not fields[0].strip():
                raise BuilderError(f"source_id is empty at taxonomy line {line_number}: {path}")
            if not fields[7].strip() or fields[7] in {"NA", "N/A", "-"}:
                raise BuilderError(f"Species field is missing at taxonomy line {line_number}: {path}")
            if validate_genomes:
                genome = Path(fields[-1]).expanduser()
                if not genome.is_absolute():
                    genome = (path.parent / genome).resolve()
                else:
                    genome = genome.resolve()
                if not genome.exists():
                    raise BuilderError(f"Genome FASTA does not exist at taxonomy line {line_number}: {genome}")
                if not genome.is_file():
                    raise BuilderError(f"Genome path is not a file at taxonomy line {line_number}: {genome}")
                fields[-1] = str(genome)
            records.append(TaxonomyRecord(fields=fields, source_id=fields[0]))
    if not records:
        raise BuilderError(f"No taxonomy records found: {path}")
    validate_source_ids(records, context=str(path))
    validate_species_taxonomy_consistency(records, context=str(path))
    return header, records


def validate_source_ids(records: Sequence[TaxonomyRecord], context: str) -> None:
    source_ids = [record.source_id for record in records]
    duplicates = sorted(source_id for source_id, count in collections.Counter(source_ids).items() if count > 1)
    if duplicates:
        preview = ", ".join(duplicates[:10])
        raise BuilderError(f"Duplicate source_id values in {context}: {preview}")


def validate_species_taxonomy_consistency(records: Sequence[TaxonomyRecord], context: str) -> None:
    species_to_signature: Dict[str, str] = {}
    for record in records:
        previous = species_to_signature.get(record.species)
        if previous is None:
            species_to_signature[record.species] = record.taxon_signature
        elif previous != record.taxon_signature:
            raise BuilderError(
                f"Inconsistent kingdom-to-species taxonomy for species {record.species!r} in {context}"
            )


def harmonize_update_taxonomy(
    old_records: Sequence[TaxonomyRecord],
    new_records: Sequence[TaxonomyRecord],
) -> Tuple[int, int]:
    """Use the existing database lineage for species already present during update."""
    old_taxonomy_by_species = {
        record.species: list(record.fields[1:8])
        for record in old_records
    }
    changed_genomes = 0
    changed_species = set()
    for record in new_records:
        authoritative_taxonomy = old_taxonomy_by_species.get(record.species)
        if authoritative_taxonomy is None or record.fields[1:8] == authoritative_taxonomy:
            continue
        record.fields[1:8] = list(authoritative_taxonomy)
        changed_genomes += 1
        changed_species.add(record.species)

    if changed_genomes:
        species_names = sorted(changed_species)
        preview = ", ".join(species_names[:5])
        if len(species_names) > 5:
            preview += f", ... (+{len(species_names) - 5} more)"
        warn(
            "Harmonized kingdom-to-species taxonomy for "
            f"{changed_genomes} new genome(s) across {len(changed_species)} species "
            f"already present in the source database ({preview}). The existing database "
            "taxonomy was treated as authoritative; new source IDs, strain fields, and "
            "genome paths were preserved."
        )
    return changed_genomes, len(changed_species)


def validate_internal_ids(records: Sequence[TaxonomyRecord]) -> bool:
    ids = [record.internal_id for record in records]
    return len(ids) == len(set(ids)) and all(re.fullmatch(r"\d{8}", value) is not None for value in ids)


def assign_ids(records: Sequence[TaxonomyRecord], start: int) -> List[Tuple[str, str]]:
    if start < 0:
        raise BuilderError("Internal-ID start cannot be negative")
    if start + len(records) > 100_000_000:
        raise BuilderError("The 8-digit MAP2B internal-ID space would be exceeded")
    mapping: List[Tuple[str, str]] = []
    for offset, record in enumerate(records):
        new_id = f"{start + offset:08d}"
        mapping.append((record.source_id, new_id))
        record.internal_id = new_id
    return mapping


def prepare_build_records(
    input_path: Path,
    id_mode: str,
) -> Tuple[List[str] | None, List[TaxonomyRecord], List[Tuple[str, str]]]:
    header, records = read_taxonomy(input_path, validate_genomes=True)
    if id_mode == "preserve":
        if not validate_internal_ids(records):
            raise BuilderError("--id-mode preserve requires unique 8-digit numeric IDs in input column 1")
        mapping = [(record.source_id, record.internal_id) for record in records]
    elif id_mode == "reassign":
        mapping = assign_ids(records, 0)
    else:
        if validate_internal_ids(records):
            mapping = [(record.source_id, record.internal_id) for record in records]
        else:
            mapping = assign_ids(records, 0)
    records.sort(key=lambda record: int(record.internal_id))
    return header, records, mapping


def normalized_output_header(header: List[str] | None) -> List[str] | None:
    if header is None:
        return None
    out = list(header)
    if out:
        out[0] = "#map2b_internal_id"
    return out


def write_taxonomy(path: Path, header: List[str] | None, records: Sequence[TaxonomyRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        out_header = normalized_output_header(header)
        if out_header:
            handle.write("\t".join(out_header) + "\n")
        for record in records:
            handle.write("\t".join(record.fields) + "\n")


def write_id_map(path: Path, mapping: Sequence[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("source_id\tmap2b_internal_id\n")
        for source_id, internal_id in mapping:
            handle.write(f"{source_id}\t{internal_id}\n")


def read_id_map(path: Path) -> List[Tuple[str, str]]:
    if not path.exists():
        raise BuilderError(
            f"Required cumulative source-ID mapping is missing: {path}. "
            "A valid genome_id_map.tsv is required for this database update."
        )
    mapping: List[Tuple[str, str]] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n\r")
            if not line.strip():
                continue
            fields = line.split("\t")
            if line_number == 1 and fields[:2] == ["source_id", "map2b_internal_id"]:
                continue
            if len(fields) < 2:
                raise BuilderError(f"Invalid ID-map line {line_number}: {path}")
            source_id, internal_id = fields[0], fields[1]
            if not source_id or re.fullmatch(r"\d{8}", internal_id) is None:
                raise BuilderError(f"Invalid ID-map entry at line {line_number}: {path}")
            mapping.append((source_id, internal_id))
    if not mapping:
        raise BuilderError(f"No ID mappings found in {path}")
    source_ids = [source_id for source_id, _ in mapping]
    internal_ids = [internal_id for _, internal_id in mapping]
    if len(source_ids) != len(set(source_ids)):
        raise BuilderError(f"Duplicate source_id values in {path}")
    if len(internal_ids) != len(set(internal_ids)):
        raise BuilderError(f"Duplicate MAP2B internal IDs in {path}")
    return mapping


def validate_mapping_against_taxonomy(
    mapping: Sequence[Tuple[str, str]], taxonomy_records: Sequence[TaxonomyRecord], path: Path
) -> None:
    mapped_ids = {internal_id for _, internal_id in mapping}
    taxonomy_ids = {record.internal_id for record in taxonomy_records}
    if mapped_ids != taxonomy_ids:
        missing_map = sorted(taxonomy_ids - mapped_ids)
        extra_map = sorted(mapped_ids - taxonomy_ids)
        detail = []
        if missing_map:
            detail.append(f"taxonomy IDs missing from map: {', '.join(missing_map[:5])}")
        if extra_map:
            detail.append(f"map IDs absent from taxonomy: {', '.join(extra_map[:5])}")
        raise BuilderError(f"{path} is not cumulative/consistent with taxonomy ({'; '.join(detail)})")


def load_or_bootstrap_update_mapping(
    source_db: Path, taxonomy_records: Sequence[TaxonomyRecord]
) -> List[Tuple[str, str]]:
    """Load an ID map or initialize one from an existing MAP2B database."""
    mapping_path = source_db / ID_MAP_FILENAME
    if mapping_path.exists():
        mapping = read_id_map(mapping_path)
        validate_mapping_against_taxonomy(mapping, taxonomy_records, mapping_path)
        return mapping

    if not validate_internal_ids(taxonomy_records):
        raise BuilderError(
            f"{ID_MAP_FILENAME} is missing and the existing taxonomy does not contain "
            "unique 8-digit MAP2B internal IDs; the database cannot be updated"
        )

    mapping = [
        (f"legacy_map2b_{record.internal_id}", record.internal_id)
        for record in taxonomy_records
    ]
    validate_mapping_against_taxonomy(mapping, taxonomy_records, mapping_path)
    warn(
        f"{ID_MAP_FILENAME} is not present; an ID map will be initialized from the existing MAP2B taxonomy."
    )
    return mapping


def shard_boundary_for_id(internal_id: int, shard_size: int) -> int:
    return ((internal_id // shard_size) + 1) * shard_size


def list_master_shards(database_dir: Path, enzyme: EnzymeSpec) -> List[Path]:
    prefix = f"{enzyme.name}.species.marisa"
    shards: List[Tuple[int, Path]] = []
    for path in database_dir.glob(prefix + "*"):
        suffix = path.name[len(prefix) :]
        if suffix.isdigit():
            shards.append((int(suffix), path))
    if not shards:
        raise BuilderError(f"No master shards matching {prefix}* were found in {database_dir}")
    return [path for _, path in sorted(shards)]


def master_shard_boundaries(database_dir: Path, enzyme: EnzymeSpec) -> List[int]:
    return [int(path.name.split(".marisa")[-1]) for path in list_master_shards(database_dir, enzyme)]


def detect_database_enzymes(database_dir: Path) -> List[EnzymeSpec]:
    detected: List[EnzymeSpec] = []
    for enzyme in ENZYMES:
        prefix = f"{enzyme.name}.species.marisa"
        if any(path.name[len(prefix) :].isdigit() for path in database_dir.glob(prefix + "*")):
            detected.append(enzyme)
    return detected


def infer_shard_size(database_dir: Path, enzyme: EnzymeSpec) -> int:
    manifest_path = database_dir / f"{enzyme.name}.species.database.json"
    if manifest_path.exists():
        try:
            with open(manifest_path, "r", encoding="utf-8") as handle:
                manifest = json.load(handle)
            value = int(manifest["shard_size"])
        except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
            raise BuilderError(f"Cannot read shard_size from manifest: {manifest_path}") from exc
        if value <= 0:
            raise BuilderError(f"Invalid shard_size in manifest: {manifest_path}")
        return value
    boundaries = master_shard_boundaries(database_dir, enzyme)
    inferred = 0
    for boundary in boundaries:
        inferred = math.gcd(inferred, boundary)
    if inferred <= 0:
        raise BuilderError(
            f"Cannot infer shard size for legacy database {database_dir}; provide -s explicitly"
        )
    warn(
        f"Manifest missing for {enzyme.name}; inferred shard size {inferred} from master-shard filenames"
    )
    return inferred


def resolve_update_shard_size(
    database_dir: Path,
    enzymes: Sequence[EnzymeSpec],
    requested_size: int | None,
) -> int:
    inferred_by_enzyme = {enzyme.name: infer_shard_size(database_dir, enzyme) for enzyme in enzymes}
    inferred_values = set(inferred_by_enzyme.values())
    if len(inferred_values) != 1:
        details = ", ".join(f"{name}={size}" for name, size in sorted(inferred_by_enzyme.items()))
        raise BuilderError(
            f"Detected enzyme databases do not share one shard size ({details}); safe shared-taxonomy update is refused"
        )
    existing_size = next(iter(inferred_values))
    if requested_size is not None and requested_size != existing_size:
        raise BuilderError(
            f"Requested -s {requested_size} does not match existing database shard size {existing_size}"
        )
    return existing_size


def save_main_shard(path: Path, pairs: Iterable[Tuple[str, bytes]], force: bool = False) -> None:
    require_marisa()
    if path.exists() and not force:
        raise BuilderError(f"Refusing to overwrite existing master shard: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_name(path.name + ".tmp")
    if temp_path.exists():
        temp_path.unlink()
    trie = marisa_trie.BytesTrie(pairs)
    trie.save(str(temp_path))
    os.replace(temp_path, path)


def build_main_shards(
    records: Sequence[TaxonomyRecord],
    enzyme: EnzymeSpec,
    output_dir: Path,
    shard_size: int,
    force: bool = False,
) -> List[Path]:
    prefix = output_dir / f"{enzyme.name}.species"
    grouped: Dict[int, List[TaxonomyRecord]] = collections.defaultdict(list)
    for record in records:
        grouped[shard_boundary_for_id(int(record.internal_id), shard_size)].append(record)

    output_files: List[Path] = []
    for boundary in sorted(grouped):
        pairs: List[Tuple[str, bytes]] = []
        log(f"Building master shard {enzyme.name}.species.marisa{boundary}")
        for record in sorted(grouped[boundary], key=lambda item: int(item.internal_id)):
            tag_counts = digest_genome(Path(record.genome_path), enzyme)
            log(f"  digested {record.internal_id} ({record.species}): {len(tag_counts):,} distinct tags")
            for tag, copy_number in tag_counts.items():
                if copy_number <= 0:
                    raise BuilderError(f"Internal non-positive copy number for genome {record.internal_id}: {copy_number}")
                key = f"{record.internal_id}{copy_number:04d}"
                pairs.append((key, tag.encode("ascii")))
        output = Path(f"{prefix}.marisa{boundary}")
        save_main_shard(output, pairs, force=force)
        output_files.append(output)
        log(f"Saved {output} with {len(pairs):,} genome-tag records")
    return output_files


def tuple_chars_to_text(value: Sequence[bytes]) -> str:
    return b"".join(value).decode("ascii")


def detect_memory_limit_bytes() -> int | None:
    """Return the smallest visible physical/cgroup memory limit, when detectable."""
    candidates: List[int] = []
    try:
        page_size = int(os.sysconf("SC_PAGE_SIZE"))
        page_count = int(os.sysconf("SC_PHYS_PAGES"))
        if page_size > 0 and page_count > 0:
            candidates.append(page_size * page_count)
    except (AttributeError, OSError, TypeError, ValueError):
        pass

    for limit_path in (
        Path("/sys/fs/cgroup/memory.max"),
        Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"),
    ):
        try:
            value = limit_path.read_text(encoding="ascii").strip()
            if value != "max":
                limit = int(value)
                # Some cgroup v1 installations use a near-2^63 value for unlimited.
                if 0 < limit < (1 << 60):
                    candidates.append(limit)
        except (OSError, UnicodeError, ValueError):
            pass
    return min(candidates) if candidates else None


def automatic_sort_memory() -> str:
    """Choose a conservative GNU-sort buffer without exposing another CLI option."""
    memory_limit = detect_memory_limit_bytes()
    if memory_limit is None:
        return "1G"
    mib = 1024 * 1024
    # Use at most 5% of visible memory, bounded to 64 MiB..4 GiB. This is the
    # total sort buffer, not a per-process allocation.
    selected = max(64 * mib, min(memory_limit // 20, 4 * 1024 * mib))
    return f"{selected // mib}M"


def unique_temp_parent(database_dir: Path) -> Path:
    """Use an explicitly configured TMPDIR, otherwise keep temp files by the DB."""
    configured = os.environ.get("TMPDIR")
    if configured:
        candidate = Path(configured).expanduser()
        if candidate.is_dir() and os.access(candidate, os.W_OK):
            return candidate.resolve()
        warn(f"TMPDIR is not a writable directory; using the database directory instead: {candidate}")
    return database_dir


def find_gnu_sort() -> str:
    """Return GNU sort executable or raise before creating large temporary files."""
    executable = shutil.which("sort")
    if executable is None:
        raise BuilderError(
            "GNU sort was not found on PATH; it is required for low-memory unique-database rebuilding"
        )
    try:
        result = subprocess.run(
            [executable, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise BuilderError(f"Cannot execute sort: {executable}") from exc
    if result.returncode != 0 or b"GNU" not in result.stdout:
        raise BuilderError(
            "The sort command on PATH is not GNU coreutils sort; "
            "GNU sort is required for bounded-memory rebuilding"
        )
    return executable


def build_global_unique_database(
    database_dir: Path,
    taxonomy_records: Sequence[TaxonomyRecord],
    enzyme: EnzymeSpec,
    output_stem: str | None = None,
    processes: int = DEFAULT_PROCESSES,
) -> Tuple[Path, Path, int, int]:
    """Build the global unique database with bounded memory and external sorting.

    Returns (marisa_path, stat_path, unique_genome_tag_records, distinct_unique_tags).
    """
    require_marisa()
    id_to_taxon = {record.internal_id: record.taxon_signature for record in taxonomy_records}
    id_unique_counts: Dict[str, int] = collections.defaultdict(int)
    distinct_unique_tags = 0
    unique_pair_count = 0

    stem = output_stem or f"{enzyme.name}.species.uniq"
    marisa_path = database_dir / f"{stem}.marisa"
    stat_path = database_dir / f"{stem}.stat.xls"
    temp_marisa = marisa_path.with_name(marisa_path.name + ".tmp")
    if temp_marisa.exists():
        temp_marisa.unlink()

    sort_executable = find_gnu_sort()
    sort_memory = automatic_sort_memory()
    temp_parent = unique_temp_parent(database_dir)
    log(
        f"Rebuilding {enzyme.name} global unique database with GNU sort: "
        f"processes={processes}, sort_buffer={sort_memory}, temp_parent={temp_parent}"
    )

    try:
        with tempfile.TemporaryDirectory(
            prefix=f".{enzyme.name}.uniq.", dir=str(temp_parent)
        ) as temp_name:
            temp_dir = Path(temp_name)
            raw_path = temp_dir / "tag_occurrences.tsv"
            sorted_path = temp_dir / "tag_occurrences.sorted.tsv"
            master_record_count = 0

            with open(raw_path, "wb", buffering=1024 * 1024) as output_handle:
                for shard in list_master_shards(database_dir, enzyme):
                    log(f"Reading master shard {shard.name}")
                    try:
                        trie = marisa_trie.BytesTrie().mmap(str(shard))
                    except Exception as exc:
                        raise BuilderError(
                            f"Cannot read {shard} as a {enzyme.name} master shard. "
                            "The shard may have an incompatible value format."
                        ) from exc
                    # iteritems() streams records. items()/keys() would build a
                    # potentially huge Python list and defeat the memory reduction.
                    for sid, tag in trie.iteritems():
                        if len(sid) < 12 or not sid.isdigit():
                            raise BuilderError(f"Invalid master-db key in {shard}: {sid!r}")
                        genome_id = sid[:8]
                        if genome_id not in id_to_taxon:
                            raise BuilderError(
                                f"Master-db genome ID {genome_id} is absent from taxonomy: {shard}"
                            )
                        if len(tag) != enzyme.length or tag.strip(b"ACGT"):
                            raise BuilderError(
                                f"Invalid tag recovered from {shard}: {tag!r}; "
                                f"expected {enzyme.length}-nt A/C/G/T"
                            )
                        sid_bytes = sid.encode("ascii")
                        output_handle.write(tag + b"\t" + sid_bytes + b"\n")
                        master_record_count += 1
                        if master_record_count % 10_000_000 == 0:
                            log(f"  streamed {master_record_count:,} genome-tag records")

            raw_size = raw_path.stat().st_size
            log(
                f"Wrote {master_record_count:,} genome-tag records "
                f"({raw_size / (1024 ** 3):.2f} GiB); sorting by tag"
            )
            environment = dict(os.environ)
            environment["LC_ALL"] = "C"
            sort_command = [
                sort_executable,
                "--field-separator=\t",
                "--key=1,1",
                f"--buffer-size={sort_memory}",
                f"--parallel={processes}",
                f"--temporary-directory={temp_dir}",
                f"--output={sorted_path}",
                str(raw_path),
            ]
            try:
                subprocess.run(sort_command, check=True, env=environment)
            except (OSError, subprocess.CalledProcessError) as exc:
                return_code = getattr(exc, "returncode", "unknown")
                raise BuilderError(
                    f"GNU sort failed (exit code {return_code}). "
                    "Check temporary-disk space and scheduler limits."
                ) from exc
            raw_path.unlink()
            log("Sorting finished; scanning equal-tag groups and building the unique MARISA database")

            def unique_pair_stream() -> Iterator[Tuple[str, bytes]]:
                nonlocal distinct_unique_tags, unique_pair_count
                current_tag: bytes | None = None
                current_taxon: str | None = None
                multiple_taxa = False
                single_copy_ids: List[bytes] = []

                def flush_group() -> Iterator[Tuple[str, bytes]]:
                    nonlocal distinct_unique_tags, unique_pair_count
                    if current_tag is None or multiple_taxa or not single_copy_ids:
                        return
                    distinct_unique_tags += 1
                    tag_text = current_tag.decode("ascii")
                    for genome_id_bytes in single_copy_ids:
                        genome_id = genome_id_bytes.decode("ascii")
                        id_unique_counts[genome_id] += 1
                        unique_pair_count += 1
                        yield tag_text, genome_id_bytes

                with open(sorted_path, "rb", buffering=1024 * 1024) as input_handle:
                    for line_number, line in enumerate(input_handle, 1):
                        try:
                            tag, sid = line.rstrip(b"\n").split(b"\t", 1)
                        except ValueError as exc:
                            raise BuilderError(
                                f"Malformed temporary sort record at line {line_number}"
                            ) from exc
                        if tag != current_tag:
                            if current_tag is not None:
                                yield from flush_group()
                            current_tag = tag
                            current_taxon = None
                            multiple_taxa = False
                            single_copy_ids = []

                        genome_id_bytes = sid[:8]
                        genome_id = genome_id_bytes.decode("ascii")
                        taxon = id_to_taxon[genome_id]
                        if current_taxon is None:
                            current_taxon = taxon
                        elif taxon != current_taxon:
                            multiple_taxa = True
                        if sid[8:] == b"0001":
                            single_copy_ids.append(genome_id_bytes)
                        if line_number % 10_000_000 == 0:
                            log(f"  scanned {line_number:,} sorted genome-tag records")

                    if current_tag is not None:
                        yield from flush_group()

            unique_trie = marisa_trie.BytesTrie(unique_pair_stream())
            unique_trie.save(str(temp_marisa))
    except Exception:
        temp_marisa.unlink(missing_ok=True)
        raise

    os.replace(temp_marisa, marisa_path)

    taxon_to_counts: Dict[str, List[int]] = collections.defaultdict(list)
    for genome_id, count in id_unique_counts.items():
        taxon_to_counts[id_to_taxon[genome_id]].append(count)
    taxon_average = {
        taxon: round(sum(values) / len(values), 4) for taxon, values in taxon_to_counts.items()
    }

    temp_stat = stat_path.with_name(stat_path.name + ".tmp")
    with open(temp_stat, "w", encoding="utf-8") as handle:
        for genome_id in sorted(id_unique_counts, key=int):
            taxon = id_to_taxon[genome_id]
            # 10 columns: ID + 7 taxonomy fields + assembly count + species average.
            handle.write(
                f"{genome_id}\t{taxon}\t{id_unique_counts[genome_id]}\t{taxon_average[taxon]}\n"
            )
    os.replace(temp_stat, stat_path)

    log(
        f"Saved global unique database with {unique_pair_count:,} genome-tag records "
        f"and {distinct_unique_tags:,} distinct unique tags"
    )
    return marisa_path, stat_path, unique_pair_count, distinct_unique_tags

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_manifest(
    database_dir: Path,
    enzyme: EnzymeSpec,
    shard_size: int,
    taxonomy_path: Path,
    mode: str,
    genome_count: int,
) -> None:
    manifest = {
        "builder": "MAP2BDatabaseBuilder.py",
        "builder_version": __version__,
        "created_at": dt.datetime.now().isoformat(timespec="seconds"),
        "mode": mode,
        "enzyme_id": enzyme.enzyme_id,
        "enzyme": enzyme.name,
        "tag_length": enzyme.length,
        "shard_size": shard_size,
        "genome_count": genome_count,
        "taxonomy_file": taxonomy_path.name,
        "taxonomy_sha256": sha256_file(taxonomy_path),
        "stat_columns": 10,
        "master_key_format": "8-digit internal_id + copy_number zero-padded to at least 4 digits",
        "source_mapping_file": ID_MAP_FILENAME,
    }
    manifest_path = database_dir / f"{enzyme.name}.species.database.json"
    temp_path = manifest_path.with_name(manifest_path.name + ".tmp")
    with open(temp_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    os.replace(temp_path, manifest_path)


def clean_generated_database_files(output_dir: Path) -> None:
    for enzyme in ENZYMES:
        for path in output_dir.glob(f"{enzyme.name}.species.marisa*"):
            if path.is_file():
                path.unlink()
        for suffix in (".species.uniq.marisa", ".species.uniq.stat.xls", ".species.database.json"):
            path = output_dir / f"{enzyme.name}{suffix}"
            if path.exists() and path.is_file():
                path.unlink()
    for name in (TAXONOMY_FILENAME, ID_MAP_FILENAME, NEW_ID_MAP_FILENAME):
        path = output_dir / name
        if path.exists() and path.is_file():
            path.unlink()


def resolve_build_enzymes(selector: int | None) -> List[EnzymeSpec]:
    return enzyme_specs_from_selector(DEFAULT_ENZYME_ID if selector is None else selector)


def resolve_existing_enzymes(
    database_dir: Path,
    selector: int | None,
    mode_name: str,
) -> List[EnzymeSpec]:
    detected = detect_database_enzymes(database_dir)
    if not detected:
        raise BuilderError(f"No recognized MAP2B master enzyme database was detected in {database_dir}")
    detected_ids = {enzyme.enzyme_id for enzyme in detected}
    if selector is None:
        return detected
    requested = enzyme_specs_from_selector(selector)
    requested_ids = {enzyme.enzyme_id for enzyme in requested}
    if not requested_ids.issubset(detected_ids):
        missing = [ENZYME_BY_ID[eid].name for eid in sorted(requested_ids - detected_ids)]
        raise BuilderError(f"Requested enzyme database(s) not found in {database_dir}: {', '.join(missing)}")
    if mode_name == "update" and requested_ids != detected_ids:
        detected_names = ", ".join(enzyme.name for enzyme in detected)
        requested_names = ", ".join(enzyme.name for enzyme in requested)
        raise BuilderError(
            "Partial update of a multi-enzyme database is refused because all enzyme databases share one taxonomy. "
            f"Detected [{detected_names}], requested [{requested_names}]. Omit -e to update all detected enzymes."
        )
    return requested


def command_build(args: argparse.Namespace) -> None:
    require_marisa()
    input_path = Path(args.input).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    shard_size = DEFAULT_SHARD_SIZE if args.shard_size is None else args.shard_size
    enzymes = resolve_build_enzymes(args.enzyme)

    if output_dir.exists() and any(output_dir.iterdir()):
        if not args.force:
            raise BuilderError(
                f"Output directory is not empty: {output_dir}; use a new directory or --force for a deliberate rebuild"
            )
        log(f"--force enabled: removing existing MAP2B database artifacts from {output_dir}")
        clean_generated_database_files(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    header, records, mapping = prepare_build_records(input_path, args.id_mode)
    taxonomy_path = output_dir / TAXONOMY_FILENAME
    write_taxonomy(taxonomy_path, header, records)
    write_id_map(output_dir / ID_MAP_FILENAME, mapping)
    new_map = output_dir / NEW_ID_MAP_FILENAME
    if new_map.exists():
        new_map.unlink()

    for enzyme in enzymes:
        build_main_shards(records, enzyme, output_dir, shard_size, force=args.force)
        build_global_unique_database(output_dir, records, enzyme, processes=args.processes)
        write_manifest(output_dir, enzyme, shard_size, taxonomy_path, "build", len(records))

    log(
        f"Database construction finished: {output_dir} "
        f"({len(records):,} genomes; enzymes: {', '.join(e.name for e in enzymes)}; shard_size={shard_size})"
    )


def copy_database_to_workdir(source: Path, destination: Path) -> Tuple[Path, Path]:
    """Copy source to a temporary sibling; return (workdir, final_destination)."""
    if destination.exists():
        raise BuilderError(f"Update output directory already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    workdir = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.building-", dir=str(destination.parent))
    )
    # mkdtemp creates a directory; copytree with dirs_exist_ok copies into it.
    shutil.copytree(source, workdir, dirs_exist_ok=True)
    return workdir, destination


def finalize_workdir(workdir: Path, destination: Path) -> None:
    os.replace(workdir, destination)


def command_update(args: argparse.Namespace) -> None:
    require_marisa()
    source_db = Path(args.database).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    if not source_db.is_dir():
        raise BuilderError(f"Existing database directory does not exist: {source_db}")

    enzymes = resolve_existing_enzymes(source_db, args.enzyme, mode_name="update")
    shard_size = resolve_update_shard_size(source_db, enzymes, args.shard_size)

    old_tax_path = source_db / TAXONOMY_FILENAME
    old_header, old_records = read_taxonomy(old_tax_path, validate_genomes=False)
    if not validate_internal_ids(old_records):
        raise BuilderError("Existing MAP2B taxonomy must contain unique 8-digit numeric IDs")

    old_mapping = load_or_bootstrap_update_mapping(source_db, old_records)
    old_source_ids = {source_id for source_id, _ in old_mapping}
    existing_internal_ids = {record.internal_id for record in old_records}

    input_path = Path(args.input).expanduser().resolve()
    new_header, new_records = read_taxonomy(input_path, validate_genomes=True)
    new_source_ids = {record.source_id for record in new_records}
    overlap = sorted(old_source_ids & new_source_ids)
    if overlap:
        raise BuilderError(
            "New input contains source_id values already present in the cumulative database mapping: "
            + ", ".join(overlap[:10])
        )

    # The source database defines taxonomy for species it already contains. This avoids
    # treating nomenclature differences (for example, GTDB versus a custom annotation)
    # as distinct species-level taxa during global unique-tag reconstruction.
    harmonize_update_taxonomy(old_records, new_records)

    # All shared-taxonomy enzyme databases must start new genomes at the same unused shard interval.
    maxima = {enzyme.name: max(master_shard_boundaries(source_db, enzyme)) for enzyme in enzymes}
    if len(set(maxima.values())) != 1:
        details = ", ".join(f"{name}={boundary}" for name, boundary in sorted(maxima.items()))
        raise BuilderError(
            f"Detected enzyme databases have different maximum shard boundaries ({details}); safe update is refused"
        )
    start_id = next(iter(maxima.values()))
    new_mapping = assign_ids(new_records, start_id)
    if any(record.internal_id in existing_internal_ids for record in new_records):
        raise BuilderError("Assigned update IDs collide with the existing database")
    new_records.sort(key=lambda record: int(record.internal_id))
    combined_records = sorted(old_records + new_records, key=lambda record: int(record.internal_id))
    validate_species_taxonomy_consistency(combined_records, context="combined old+new taxonomy")
    cumulative_mapping = list(old_mapping) + list(new_mapping)

    in_place = bool(args.in_place)
    if in_place:
        if output_dir != source_db:
            raise BuilderError("--in-place requires -o/--output-dir to be the same directory as -d/--database")
        workdir = source_db
        final_destination = source_db
        warn("--in-place modifies the existing database directly; interruption can leave newly added master shards behind")
    else:
        workdir, final_destination = copy_database_to_workdir(source_db, output_dir)

    try:
        # Build new immutable master shards first.
        for enzyme in enzymes:
            build_main_shards(new_records, enzyme, workdir, shard_size, force=False)

        # Recompute global specificity against the complete old+new master set.
        rebuilt_outputs: List[Tuple[EnzymeSpec, Path, Path]] = []
        for enzyme in enzymes:
            temp_stem = f"{enzyme.name}.species.uniq.rebuild"
            temp_marisa, temp_stat, _, _ = build_global_unique_database(
                workdir, combined_records, enzyme, output_stem=temp_stem, processes=args.processes
            )
            rebuilt_outputs.append((enzyme, temp_marisa, temp_stat))

        # Replace active unique databases only after every enzyme rebuilt successfully.
        for enzyme, temp_marisa, temp_stat in rebuilt_outputs:
            os.replace(temp_marisa, workdir / f"{enzyme.name}.species.uniq.marisa")
            os.replace(temp_stat, workdir / f"{enzyme.name}.species.uniq.stat.xls")

        taxonomy_path = workdir / TAXONOMY_FILENAME
        write_taxonomy(taxonomy_path, old_header or new_header, combined_records)
        write_id_map(workdir / NEW_ID_MAP_FILENAME, new_mapping)
        write_id_map(workdir / ID_MAP_FILENAME, cumulative_mapping)
        for enzyme in enzymes:
            write_manifest(workdir, enzyme, shard_size, taxonomy_path, "update", len(combined_records))

        if not in_place:
            finalize_workdir(workdir, final_destination)
            workdir = final_destination

    except Exception:
        if not in_place and workdir.exists():
            shutil.rmtree(workdir, ignore_errors=True)
        raise

    log(
        f"Database update finished: {workdir} "
        f"(+{len(new_records):,} genomes; total {len(combined_records):,}; "
        f"enzymes: {', '.join(e.name for e in enzymes)}; shard_size={shard_size})"
    )


def count_and_validate_master_records(
    database_dir: Path,
    taxonomy_ids: set[str],
    enzyme: EnzymeSpec,
    shard_size: int,
) -> Tuple[int, int]:
    require_marisa()
    fmt = f"{enzyme.length}c"
    shard_count = 0
    record_count = 0
    for shard in list_master_shards(database_dir, enzyme):
        boundary = int(shard.name.split(".marisa")[-1])
        if boundary <= 0 or boundary % shard_size != 0:
            raise BuilderError(
                f"Master shard boundary {boundary} is incompatible with shard_size {shard_size}: {shard}"
            )
        try:
            trie = marisa_trie.RecordTrie(fmt).mmap(str(shard))
        except Exception as exc:
            raise BuilderError(f"Cannot open master shard with RecordTrie('{fmt}'): {shard}") from exc
        unique_keys = set(trie.keys())
        shard_count += 1
        for sid in unique_keys:
            if re.fullmatch(r"\d{12,}", sid) is None:
                raise BuilderError(f"Invalid master-db key in {shard}: {sid!r}")
            genome_id = sid[:8]
            copy_number = int(sid[8:])
            if genome_id not in taxonomy_ids:
                raise BuilderError(f"Master-db genome ID {genome_id} is absent from taxonomy: {shard}")
            if copy_number < 1:
                raise BuilderError(f"Invalid copy number in master-db key {sid!r}: {shard}")
            if shard_boundary_for_id(int(genome_id), shard_size) != boundary:
                raise BuilderError(f"Master-db key {sid!r} is stored in the wrong shard {shard.name}")
            values = trie[sid]
            record_count += len(values)
            for value in values:
                tag = tuple_chars_to_text(value)
                if len(tag) != enzyme.length or re.fullmatch(r"[ACGT]+", tag) is None:
                    raise BuilderError(f"Invalid {enzyme.name} tag in {shard}: {tag!r}")
    return shard_count, record_count


def validate_unique_database(
    database_dir: Path,
    taxonomy_by_id: Mapping[str, TaxonomyRecord],
    enzyme: EnzymeSpec,
) -> Tuple[int, int, Dict[str, int]]:
    require_marisa()
    uniq_path = database_dir / f"{enzyme.name}.species.uniq.marisa"
    if not uniq_path.exists():
        raise BuilderError(f"Global unique database is missing: {uniq_path}")
    try:
        uniq = marisa_trie.BytesTrie().mmap(str(uniq_path))
    except Exception as exc:
        raise BuilderError(f"Cannot open global unique database: {uniq_path}") from exc
    distinct_tags = set(uniq.keys())
    relation_count = 0
    per_genome: Dict[str, int] = collections.defaultdict(int)
    for tag in distinct_tags:
        if len(tag) != enzyme.length or re.fullmatch(r"[ACGT]+", tag) is None:
            raise BuilderError(f"Invalid tag in global unique database: {tag!r}")
        ids: List[str] = []
        for genome_id_bytes in uniq[tag]:
            try:
                genome_id = genome_id_bytes.decode("ascii")
            except AttributeError:
                genome_id = bytes(genome_id_bytes).decode("ascii")
            if re.fullmatch(r"\d{8}", genome_id) is None:
                raise BuilderError(f"Invalid genome ID in global unique database: {genome_id!r}")
            if genome_id not in taxonomy_by_id:
                raise BuilderError(f"Unique-db genome ID {genome_id} is absent from taxonomy")
            ids.append(genome_id)
            relation_count += 1
            per_genome[genome_id] += 1
        taxa = {taxonomy_by_id[genome_id].taxon_signature for genome_id in ids}
        if len(taxa) > 1:
            raise BuilderError(f"Tag {tag!r} in unique database maps to more than one species taxonomy")
    return relation_count, len(distinct_tags), dict(per_genome)


def validate_stat_file(
    database_dir: Path,
    taxonomy_by_id: Mapping[str, TaxonomyRecord],
    enzyme: EnzymeSpec,
    unique_relation_count: int,
    unique_per_genome: Mapping[str, int],
) -> int:
    stat_path = database_dir / f"{enzyme.name}.species.uniq.stat.xls"
    if not stat_path.exists():
        raise BuilderError(f"Statistics file is missing: {stat_path}")
    seen_ids: set[str] = set()
    parsed_counts: Dict[str, int] = {}
    parsed_averages: Dict[str, float] = {}
    taxon_to_counts: Dict[str, List[int]] = collections.defaultdict(list)

    with open(stat_path, "r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n\r")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) not in (10, 11):
                raise BuilderError(
                    f"Invalid statistics line {line_number}: expected 10 or 11 columns, got {len(fields)}"
                )
            genome_id = fields[0]
            if genome_id in seen_ids:
                raise BuilderError(f"Duplicate genome ID in statistics file at line {line_number}: {genome_id}")
            if genome_id not in taxonomy_by_id:
                raise BuilderError(f"Statistics genome ID is absent from taxonomy at line {line_number}: {genome_id}")
            expected_taxon_fields = taxonomy_by_id[genome_id].fields[1:8]
            if fields[1:8] != expected_taxon_fields:
                raise BuilderError(f"Statistics taxonomy mismatch for genome ID {genome_id} at line {line_number}")
            try:
                count = int(fields[8])
                average = float(fields[9])
                legacy_extra = float(fields[10]) if len(fields) == 11 else None
            except ValueError as exc:
                raise BuilderError(f"Non-numeric statistics value at line {line_number}: {stat_path}") from exc
            if count < 0 or average < 0 or (legacy_extra is not None and legacy_extra < 0):
                raise BuilderError(f"Negative statistics value at line {line_number}: {stat_path}")
            seen_ids.add(genome_id)
            parsed_counts[genome_id] = count
            parsed_averages[genome_id] = average
            taxon_to_counts[taxonomy_by_id[genome_id].taxon_signature].append(count)

    if sum(parsed_counts.values()) != unique_relation_count:
        raise BuilderError(
            f"statistics/unique-database count mismatch for {enzyme.name}: "
            f"stat column 9 sum={sum(parsed_counts.values())}, unique genome-tag records={unique_relation_count}"
        )
    if parsed_counts != dict(unique_per_genome):
        raise BuilderError(f"Statistics per-genome counts do not match the unique database for {enzyme.name}")

    expected_taxon_average = {
        taxon: round(sum(values) / len(values), 4) for taxon, values in taxon_to_counts.items()
    }
    for genome_id, observed_average in parsed_averages.items():
        taxon = taxonomy_by_id[genome_id].taxon_signature
        if abs(observed_average - expected_taxon_average[taxon]) > 1e-4:
            raise BuilderError(
                f"Statistics species-average mismatch for genome ID {genome_id}: "
                f"observed {observed_average}, expected {expected_taxon_average[taxon]}"
            )
    return len(seen_ids)


def validate_manifest(
    database_dir: Path,
    taxonomy_path: Path,
    taxonomy_count: int,
    enzyme: EnzymeSpec,
    shard_size: int,
) -> None:
    manifest_path = database_dir / f"{enzyme.name}.species.database.json"
    if not manifest_path.exists():
        warn(f"Manifest not found for {enzyme.name}; structural validation continues: {manifest_path}")
        return
    try:
        with open(manifest_path, "r", encoding="utf-8") as handle:
            manifest = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise BuilderError(f"Cannot read manifest: {manifest_path}") from exc
    checks = {
        "enzyme_id": enzyme.enzyme_id,
        "enzyme": enzyme.name,
        "tag_length": enzyme.length,
        "shard_size": shard_size,
        "genome_count": taxonomy_count,
        "taxonomy_file": taxonomy_path.name,
    }
    for key, expected in checks.items():
        if key in manifest and manifest[key] != expected:
            raise BuilderError(
                f"Manifest mismatch for {enzyme.name}: {key}={manifest[key]!r}, expected {expected!r}"
            )
    if "taxonomy_sha256" in manifest:
        observed = sha256_file(taxonomy_path)
        if manifest["taxonomy_sha256"] != observed:
            raise BuilderError(
                f"Manifest taxonomy SHA-256 mismatch for {enzyme.name}: "
                f"{manifest['taxonomy_sha256']} != {observed}"
            )


def validate_one_enzyme(database_dir: Path, records: Sequence[TaxonomyRecord], enzyme: EnzymeSpec) -> None:
    taxonomy_by_id = {record.internal_id: record for record in records}
    shard_size = infer_shard_size(database_dir, enzyme)
    shard_count, master_records = count_and_validate_master_records(
        database_dir, set(taxonomy_by_id), enzyme, shard_size
    )
    unique_records, distinct_unique_tags, per_genome = validate_unique_database(
        database_dir, taxonomy_by_id, enzyme
    )
    stat_rows = validate_stat_file(
        database_dir, taxonomy_by_id, enzyme, unique_records, per_genome
    )
    validate_manifest(
        database_dir,
        database_dir / TAXONOMY_FILENAME,
        len(records),
        enzyme,
        shard_size,
    )
    log(
        f"Validation passed for {enzyme.name}: {len(records):,} genomes, {shard_count} master shards, "
        f"{master_records:,} master records, {unique_records:,} unique genome-tag records, "
        f"{distinct_unique_tags:,} distinct unique tags, {stat_rows:,} stat rows"
    )


def command_validate(args: argparse.Namespace) -> None:
    require_marisa()
    database_dir = Path(args.database).expanduser().resolve()
    if not database_dir.is_dir():
        raise BuilderError(f"Database directory does not exist: {database_dir}")
    taxonomy_path = database_dir / TAXONOMY_FILENAME
    _, records = read_taxonomy(taxonomy_path, validate_genomes=False)
    if not validate_internal_ids(records):
        raise BuilderError("Taxonomy does not contain unique 8-digit MAP2B internal IDs")
    enzymes = resolve_existing_enzymes(database_dir, args.enzyme, mode_name="validate")

    mapping_path = database_dir / ID_MAP_FILENAME
    if mapping_path.exists():
        mapping = read_id_map(mapping_path)
        validate_mapping_against_taxonomy(mapping, records, mapping_path)
    else:
        warn(f"{ID_MAP_FILENAME} is not present; validating the MAP2B database from its taxonomy and tag files")

    for enzyme in enzymes:
        validate_one_enzyme(database_dir, records, enzyme)


def print_enzyme_table() -> None:
    print("ID\tEnzyme\tTag_length")
    for enzyme in ENZYMES:
        print(f"{enzyme.enzyme_id}\t{enzyme.name}\t{enzyme.length}")
    print(f"{ALL_ENZYME_ID}\tAllEnzyme\t-")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build, update, or validate a MAP2B-compatible species-specific 2b-tag database.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--list-enzymes", action="store_true", help="List enzyme IDs/names/tag lengths and exit.")
    parser.add_argument(
        "-m",
        "--mode",
        type=int,
        choices=[1, 2, 3],
        help="Run mode: 1=build, 2=update, 3=validate.",
    )
    parser.add_argument("-i", "--input", help="Taxonomy/genome input table (plain text or .gz).")
    parser.add_argument("-d", "--database", help="Existing database directory for mode 2/3.")
    parser.add_argument("-o", "--output-dir", help="Output database directory for mode 1/2.")
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=DEFAULT_PROCESSES,
        help="CPU threads used by global unique-database sorting (default: 1).",
    )
    parser.add_argument(
        "-e",
        "--enzyme",
        type=parse_enzyme_selector,
        help=(
            "Enzyme ID/name. 17/AllEnzyme means all 16 separately. "
            "Mode 1 default: 13(CjePI). Mode 2/3 default: auto-detect existing database enzyme(s)."
        ),
    )
    parser.add_argument(
        "-s",
        "--shard-size",
        type=int,
        default=None,
        help=(
            f"Master-shard internal-ID interval. Mode 1 default: {DEFAULT_SHARD_SIZE}. "
            "Mode 2 default: read manifest/infer existing value."
        ),
    )
    parser.add_argument(
        "--id-mode",
        choices=["auto", "preserve", "reassign"],
        default="auto",
        help="Mode 1 internal-ID handling (default: auto).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Mode 1 only: remove existing recognized MAP2B database artifacts before deliberate rebuild.",
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Mode 2 only: update the existing database directory directly; not recommended.",
    )
    return parser


def require_arg(parser: argparse.ArgumentParser, args: argparse.Namespace, attr: str, flag: str) -> None:
    if not getattr(args, attr):
        parser.error(f"{flag} is required for mode {args.mode}")


def validate_mode_arguments(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if args.mode is None:
        parser.error("-m/--mode is required: 1=build, 2=update, 3=validate")
    if args.shard_size is not None and args.shard_size <= 0:
        parser.error("-s/--shard-size must be positive")
    if args.processes <= 0:
        parser.error("-p/--processes must be positive")

    if args.mode == 1:
        require_arg(parser, args, "input", "-i/--input")
        require_arg(parser, args, "output_dir", "-o/--output-dir")
        if args.database:
            parser.error("-d/--database is not used in mode 1")
        if args.in_place:
            parser.error("--in-place is only valid in mode 2")
    elif args.mode == 2:
        require_arg(parser, args, "database", "-d/--database")
        require_arg(parser, args, "input", "-i/--input")
        require_arg(parser, args, "output_dir", "-o/--output-dir")
        if args.force:
            parser.error("--force is only valid in mode 1")
        if args.id_mode != "auto":
            parser.error("--id-mode is only used in mode 1")
    elif args.mode == 3:
        require_arg(parser, args, "database", "-d/--database")
        if args.input:
            parser.error("-i/--input is not used in mode 3")
        if args.output_dir:
            parser.error("-o/--output-dir is not used in mode 3")
        if args.force or args.in_place:
            parser.error("--force/--in-place are not valid in mode 3")
        if args.id_mode != "auto":
            parser.error("--id-mode is only used in mode 1")


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.list_enzymes:
        print_enzyme_table()
        return 0
    validate_mode_arguments(parser, args)

    try:
        if args.mode == 1:
            command_build(args)
        elif args.mode == 2:
            command_update(args)
        else:
            command_validate(args)
    except BuilderError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
