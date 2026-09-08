# i2bRAD-M

Computational companion repository for the Nature Protocols article:

*An Integrated 2bRAD-M (i2bRAD-M) Approach for High-Resolution, Quantitative and Cross-Domain Microbial Profiling in Challenging Environmental and Biomedical Samples*

## Overview

This repository contains the software, environment specification, input templates, and usage documentation for the computational components of the i2bRAD-M protocol. The available analysis branches are:

1. species-level microbial profiling with MAP2B;
2. construction, extension, and validation of MAP2B-compatible databases;
3. host-referenced absolute abundance estimation;
4. host SNP and host methylation analysis;
5. serial-tag preprocessing for five-tag libraries;
6. metatranscriptomic species profiling and RNA/DNA integration; and
7. cross-domain profiling of microbial and macro-organism eDNA.

Use the accompanying article for experimental design, library preparation, and the recommended analysis sequence. Use this repository for the executable commands and software interfaces.

## Repository layout

```text
i2bRAD-M/
├── README.md
├── environment/                  reproducible Linux Conda environment
├── software/
│   ├── map2b/                    MAP2B and cross-domain analysis
│   ├── database_builder/         database construction and validation
│   ├── absolute_quantification/  host-referenced abundance conversion
│   ├── host_snp/                 Type IIB host SNP workflow
│   ├── methylation/              host methylation workflow
│   ├── serial_tag/               five-tag serial-library preprocessing
│   └── metatranscriptome/        RNA profiling and RNA/DNA integration
├── docs/                         workflow-specific documentation
├── examples/manifests/           tab-delimited input templates
├── data/                         data-availability guidance
├── tests/                        repository and interface checks
└── SHA256SUMS.txt                file-integrity checksums
```

The root README is the main entry point. Detailed instructions are organized by workflow under `docs/` rather than repeated in every software directory.

## Installation

The supplied environment targets **Linux x86_64 with glibc >= 2.17**. From the repository root:

```bash
cd environment
bash install.sh
conda activate i2bRAD-M
./check_environment.sh
cd ..
```

`environment/install.sh` uses the exact Linux package lock in `environment/env/` and installs the included SOAP2 executables after checksum verification. The readable `environment/environment.yml` is provided for dependency inspection; the exact lock is the reproducible installation route. See [`environment/README.md`](environment/README.md).

## Workflow guides

| Analysis | Main program(s) | Documentation |
|---|---|---|
| Microbial relative abundance | `software/map2b/bin/MAP2B.py` | [`docs/INPUT_OUTPUT.md`](docs/INPUT_OUTPUT.md) |
| Database construction and validation | `software/database_builder/MAP2BDatabaseBuilder.py` | [`docs/DATABASE_BUILDER.md`](docs/DATABASE_BUILDER.md) |
| Absolute abundance | `software/absolute_quantification/MAP2BAbsoluteQuantifier.py` | [`docs/ABSOLUTE_QUANTIFICATION.md`](docs/ABSOLUTE_QUANTIFICATION.md) |
| Host SNPs | Perl programs in `software/host_snp/` | [`docs/HOST_SNP.md`](docs/HOST_SNP.md) |
| Host methylation | programs in `software/methylation/` | [`docs/METHYLATION.md`](docs/METHYLATION.md) |
| Serial-tag preprocessing | programs in `software/serial_tag/` | [`docs/SERIAL_TAG.md`](docs/SERIAL_TAG.md) |
| Metatranscriptome profiling | `MTKrakenProfiler.py` and `RNADNARatio.py` | [`docs/METATRANSCRIPTOME.md`](docs/METATRANSCRIPTOME.md) |
| Cross-domain eDNA | `software/map2b/bin/MAP2B-Cross-domain.py` | [`docs/CROSS_DOMAIN.md`](docs/CROSS_DOMAIN.md) |

## Core analysis sequence

Run commands from the repository root after activating `i2bRAD-M`. Replace example paths and enzyme identifiers with values appropriate for the library and reference database.

### 1. Build and validate a database

```bash
python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 1 -i /path/genome_manifest.tsv.gz -e <enzyme_id> -o /path/database

python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 3 -d /path/database

# Optional: update an existing database; -p accelerates external sorting
python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 2 -d /path/base_database -i /path/additional_genomes.tsv.gz \
  -o /path/updated_database -p 8
```

### 2. Run MAP2B

```bash
python3 software/map2b/bin/MAP2B.py \
  -i examples/manifests/map2b_samples.tsv \
  -e <enzyme_id> -s /absolute/path/database \
  -o results/map2b -p 4
```

The enzyme identifier must match both the library preparation and database construction. The i2bRAD-M MAP2B interface supports enzyme IDs `1`-`16`, and all 16 enzymes can be analyzed using G-score filtering with `-g` and a matching database. For BcgI (`5`) and CjePI (`13`), i2bRAD-M additionally provides a Random Forest-based filtering workflow. To facilitate standard analyses, pre-built database download lists are provided for BsaXI (`3`), BcgI (`5`), and CjePI (`13`); use `-s` to provide a matching custom database for the other enzymes. Follow the article when selecting the G-score threshold.

### 3. Optional absolute quantification

```bash
python3 software/absolute_quantification/MAP2BAbsoluteQuantifier.py \
  -i results/map2b/Abundance.xls \
  -m examples/manifests/absolute_quantification_metadata.tsv \
  -s <host_species_exactly_as_in_Species_column> \
  -o results/absolute_quantification
```

### 4. Methylation analysis

The protocol route uses the default methylation preset, so the default enzyme option is omitted:

```bash
python3 software/methylation/build_reference_reads.py \
  -r /path/reference.fa.gz

python3 software/methylation/extract_fixed_reads.py \
  -i examples/manifests/methylation_samples.tsv

python3 software/methylation/map_and_quantify_reads.py \
  -i examples/manifests/methylation_samples.tsv
```

Additional supported options are documented in [`docs/METHYLATION.md`](docs/METHYLATION.md).

### 5. Optional serial-tag preprocessing

```bash
python3 software/serial_tag/serial2brad.py \
  -i examples/manifests/serial_tag_samples.tsv \
  -e <enzyme_id> -o results/serial_tag -p 4
```

The serial-tag workflow supports 16 Type IIB restriction enzymes. See [`docs/SERIAL_TAG.md`](docs/SERIAL_TAG.md) for the stepwise protocol route and output assignment.

### 6. Optional metatranscriptome analysis

```bash
python3 software/metatranscriptome/MTKrakenProfiler.py \
  -i examples/manifests/metatranscriptome_samples.tsv \
  -d /path/kraken2_database -o results/RNA_profile -p 8

python3 software/metatranscriptome/RNADNARatio.py \
  -d results/map2b/Abundance.xls \
  -r results/RNA_profile/RNA_species_abundance.tsv \
  -o results/RNA_DNA_ratio
```

Kraken2 is required for RNA taxonomic profiling. Quality control and rRNA removal are upstream preprocessing steps and may be performed with the tools selected for the study.

### 7. Cross-domain analysis

```bash
python3 software/map2b/bin/MAP2B-Cross-domain.py \
  -i examples/manifests/map2b_samples.tsv \
  -e <enzyme_id> -s /absolute/path/cross_domain_database \
  -o results/cross_domain -p 4
```

## Input templates and data

Tab-delimited input templates are provided in `examples/manifests/`. Replace the example paths and identifiers with values for the current analysis. Large sequencing datasets and reference databases are maintained outside this source-code repository; consult the article's data-availability statement and [`data/README.md`](data/README.md).

## Repository checks

```bash
python3 tests/run_checks.py
```

The checks verify Python syntax and command-line interfaces, check Perl syntax when Perl is available, confirm key file integrity, and detect common repository artifacts. Use the workflow guides for analysis of study-specific input data.

## Integrity verification

From the repository root on Linux:

```bash
sha256sum -c SHA256SUMS.txt
```

## Citation and use

Please cite the accompanying Nature Protocols article when using this workflow. The repository contents are provided under the terms in [`LICENSE`](LICENSE); bundled third-party components remain subject to their respective terms as described in [`environment/THIRD_PARTY.md`](environment/THIRD_PARTY.md).
