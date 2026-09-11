#!/usr/bin/env python3
"""Non-invasive integrity and interface checks for the i2bRAD-M repository."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import py_compile
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_READER_SHA256 = "e6e035e06767425765eaeee29d243b1de097a43415eff87fcd571342f003ae06"
EXPECTED_DATABASE_BUILDER_SHA256 = "0c2393561da099d8e7abd8d0816d270e903d6b5ef22c28001d683c8ebd3b1c87"
PYTHON_PROGRAMS = [
    ROOT / "software/map2b/bin/MAP2B.py",
    ROOT / "software/map2b/bin/MAP2B-Cross-domain.py",
    ROOT / "software/map2b/scripts/split_abundance.py",
    ROOT / "software/database_builder/MAP2BDatabaseBuilder.py",
    ROOT / "software/absolute_quantification/MAP2BAbsoluteQuantifier.py",
    ROOT / "software/methylation/build_reference_reads.py",
    ROOT / "software/methylation/extract_fixed_reads.py",
    ROOT / "software/methylation/map_and_quantify_reads.py",
    ROOT / "software/serial_tag/serial2brad.py",
    ROOT / "software/serial_tag/scripts/extractN.py",
    ROOT / "software/serial_tag/scripts/goodquality.py",
    ROOT / "software/serial_tag/scripts/common_find.py",
    ROOT / "software/serial_tag/scripts/pear_merge.py",
    ROOT / "software/serial_tag/scripts/fq2fasta.py",
    ROOT / "software/serial_tag/scripts/split5.py",
    ROOT / "software/serial_tag/scripts/recover_sample_tags.py",
    ROOT / "software/metatranscriptome/MTKrakenProfiler.py",
    ROOT / "software/metatranscriptome/RNADNARatio.py",
]
PERL_PROGRAMS = sorted((ROOT / "software/host_snp").glob("*.pl"))
failures: list[str] = []
passes = 0
warnings = 0


def report(ok: bool, label: str, detail: str = "") -> None:
    global passes
    suffix = f": {detail}" if detail else ""
    print(f"{'PASS' if ok else 'FAIL'}  {label}{suffix}")
    if ok:
        passes += 1
    else:
        failures.append(label + suffix)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


for relative in [
    "README.md",
    "LICENSE",
    "environment/install.sh",
    "environment/env/i2bRAD-M-linux-64.validated.lock.txt",
    "software/map2b/bin/MAP2B.py",
    "software/database_builder/MAP2BDatabaseBuilder.py",
    "software/absolute_quantification/MAP2BAbsoluteQuantifier.py",
    "software/host_snp/ENZYME_TABLE.tsv",
    "software/methylation/PRESET_TABLE.tsv",
    "docs/multi_omics/METHYLATION.md",
    "docs/edna/SERIAL_TAG.md",
    "docs/multi_omics/METATRANSCRIPTOME.md",
    "software/serial_tag/config/enzymes.json",
    "examples/standard_map2b/map2b_samples.tsv",
    "examples/edna/serial_tag_samples.tsv",
    "examples/multi_omics/metatranscriptome_samples.tsv",
]:
    report((ROOT / relative).is_file(), f"required file {relative}")

with tempfile.TemporaryDirectory(prefix="i2bradm-check-") as tmp:
    for index, program in enumerate(PYTHON_PROGRAMS):
        try:
            py_compile.compile(str(program), cfile=str(Path(tmp) / f"{index}.pyc"), doraise=True)
            report(True, f"Python syntax {program.relative_to(ROOT)}")
        except Exception as exc:
            report(False, f"Python syntax {program.relative_to(ROOT)}", str(exc))

env = dict(os.environ)
env["PYTHONDONTWRITEBYTECODE"] = "1"
for program in PYTHON_PROGRAMS:
    completed = subprocess.run(
        [sys.executable, str(program), "--help"],
        cwd=ROOT,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    report(
        completed.returncode == 0 and "usage:" in completed.stdout.lower(),
        f"help smoke {program.relative_to(ROOT)}",
        f"exit={completed.returncode}" if completed.returncode else "",
    )

perl = shutil.which("perl")
if perl:
    for program in PERL_PROGRAMS:
        completed = subprocess.run(
            [perl, "-c", str(program)],
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        last = completed.stdout.strip().splitlines()[-1] if completed.stdout.strip() else ""
        report(completed.returncode == 0, f"Perl syntax {program.relative_to(ROOT)}", last)
else:
    warnings += 1
    print("WARN  Perl unavailable; host-SNP syntax checks skipped")

reader = ROOT / "software/map2b/scripts/CalculateRelativeAbundance_Single2bEnzyme.py"
actual = sha256(reader)
report(actual == EXPECTED_READER_SHA256, "MAP2B abundance-reader integrity", actual)

builder = ROOT / "software/database_builder/MAP2BDatabaseBuilder.py"
builder_actual = sha256(builder)
report(
    builder_actual == EXPECTED_DATABASE_BUILDER_SHA256,
    "MAP2B database-builder integrity",
    builder_actual,
)

serial_config = ROOT / "software/serial_tag/config/enzymes.json"
try:
    import json
    enzyme_config = json.loads(serial_config.read_text(encoding="utf-8-sig"))
    enzyme_ids = {int(item["id"]) for item in enzyme_config.get("enzymes", [])}
    report(
        len(enzyme_config.get("enzymes", [])) == 16 and enzyme_ids == set(range(1, 17)),
        "serial-tag 16-enzyme configuration",
    )
except Exception as exc:
    report(False, "serial-tag 16-enzyme configuration", str(exc))

for manifest, expected_header in [
    (ROOT / "examples/edna/serial_tag_samples.tsv", "library_id\tr1\tr2\tt1\tt2\tt3\tt4\tt5"),
    (ROOT / "examples/multi_omics/metatranscriptome_samples.tsv", "sample_id\tread1\tread2"),
]:
    first_line = manifest.read_text(encoding="utf-8-sig").splitlines()[0] if manifest.is_file() else ""
    report(first_line == expected_header, f"manifest header {manifest.relative_to(ROOT)}")

map2b = ROOT / "software/map2b"
report(
    all((map2b / name).exists() for name in ("bin", "scripts", "config", "database")),
    "cross-domain relative MAP2B layout",
)

cache_items = [p for p in ROOT.rglob("*") if p.name == "__pycache__" or p.suffix == ".pyc"]
report(
    not cache_items,
    "no Python cache artifacts",
    ", ".join(str(p.relative_to(ROOT)) for p in cache_items[:5]),
)
large = [
    p
    for p in ROOT.rglob("*")
    if p.is_file() and ".git" not in p.parts and p.stat().st_size > 100 * 1024 * 1024
]
report(
    not large,
    "no files over GitHub 100 MiB limit",
    ", ".join(str(p.relative_to(ROOT)) for p in large),
)

print(f"\nSummary: {passes} passed, {len(failures)} failed, {warnings} warning(s).")
if failures:
    raise SystemExit(1)