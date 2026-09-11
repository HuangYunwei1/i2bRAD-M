# i2bRAD-M

This repository provides the computational resources accompanying the i2bRAD-M protocol:

*Integrated 2bRAD-M Approach for Comprehensive and Cost-efficient Metagenomic Profiling of Challenging Environmental and Biomedical Specimens*

## Overview

This repository contains the software, environment specification, input templates, and usage documentation for the computational components of the i2bRAD-M protocol. Following the organization of the accompanying article, the workflow comprises:

1. **Reference database construction**
   - construction or download of a microbial MAP2B reference database;
   - extension of a microbial database with host or other macro-organism genomes for holo-genome, absolute-quantification, or eDNA analysis.
2. **Standard MAP2B analysis**
   - species-level microbial profiling;
   - optional host-referenced absolute quantification.
3. **Extended module A: Multi-omics analysis**
   - metatranscriptomic species profiling and RNA/DNA integration;
   - epigenome analysis (MethylRAD-based).
4. **Extended module B: Holo-genome analysis**
   - host-SNP genotype analysis;
   - microbiome identification with MAP2B.
5. **Extended module C: Serial sequencing analysis for eDNA profiling**
   - cross-domain profiling of microbial and macro-organism eDNA;
   - direct analysis of single-tag libraries and preprocessing of five-tag serial libraries.

Use the accompanying article for experimental design, library preparation, and the recommended analysis sequence. Use this repository for the executable commands and software interfaces.

## Repository layout

```text
i2bRAD-M/
├── README.md
├── environment/                  reproducible Linux Conda environment
├── software/
│   ├── database_builder/         database construction and validation
│   ├── map2b/                    MAP2B and cross-domain analysis
│   ├── absolute_quantification/  host-referenced abundance conversion
│   ├── metatranscriptome/        RNA profiling and RNA/DNA integration
│   ├── methylation/              MethylRAD-based epigenome workflow
│   ├── host_snp/                 Type IIB host SNP workflow
│   └── serial_tag/               five-tag serial-library preprocessing
├── docs/                         module-based workflow documentation
│   ├── 01_Reference_database_construction/
│   ├── 02_Standard_MAP2B/
│   │   ├── MAP2B/
│   │   └── Host-referenced_absolute_quantification/
│   ├── 03_Multi-omics/
│   │   ├── Metatranscriptome/
│   │   └── Epigenome/
│   ├── 04_Holo-genome/
│   │   └── Host-SNP_genotype_analysis/
│   └── 05_eDNA/
│       ├── Cross-domain_profiling/
│       └── Serial_sequencing/
├── examples/                     module-based templates and example tables
│   ├── 01_Reference_database_construction/
│   ├── 02_Standard_MAP2B/
│   │   ├── MAP2B/
│   │   └── Host-referenced_absolute_quantification/
│   ├── 03_Multi-omics/
│   │   ├── Metatranscriptome/
│   │   └── Epigenome/
│   ├── 04_Holo-genome/
│   │   └── Host-SNP_genotype_analysis/
│   └── 05_eDNA/
│       ├── Cross-domain_profiling/
│       └── Serial_sequencing/
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

`environment/install.sh` creates the unified environment from the supplied Linux package lock and installs the included SOAP2 executables. The readable `environment/environment.yml` lists the direct dependencies. See [`environment/README.md`](environment/README.md).

## Workflow guides

| Article section or module | Analysis | Main program(s) | Documentation |
|---|---|---|---|
| Reference database construction | Build, update, and validate MAP2B databases | `software/database_builder/MAP2BDatabaseBuilder.py` | [`docs/01_Reference_database_construction/DATABASE_BUILDER.md`](docs/01_Reference_database_construction/DATABASE_BUILDER.md) |
| Standard MAP2B analysis | Microbial relative abundance | `software/map2b/bin/MAP2B.py` | [`docs/02_Standard_MAP2B/MAP2B/MAP2B.md`](docs/02_Standard_MAP2B/MAP2B/MAP2B.md) |
| (Optional) Host-referenced absolute quantification | Host-referenced absolute quantification | `software/absolute_quantification/MAP2BAbsoluteQuantifier.py` | [`docs/02_Standard_MAP2B/Host-referenced_absolute_quantification/ABSOLUTE_QUANTIFICATION.md`](docs/02_Standard_MAP2B/Host-referenced_absolute_quantification/ABSOLUTE_QUANTIFICATION.md) |
| Extended module A: Multi-omics | Metatranscriptome analysis | programs in `software/metatranscriptome/` | [`docs/03_Multi-omics/Metatranscriptome/METATRANSCRIPTOME.md`](docs/03_Multi-omics/Metatranscriptome/METATRANSCRIPTOME.md) |
| Extended module A: Multi-omics | Epigenome analysis (MethylRAD-based) | programs in `software/methylation/` | [`docs/03_Multi-omics/Epigenome/EPIGENOME.md`](docs/03_Multi-omics/Epigenome/EPIGENOME.md) |
| Extended module B: Holo-genome | Host-SNP genotype analysis and microbiome identification | programs in `software/host_snp/` and `software/map2b/bin/MAP2B.py` | [`docs/04_Holo-genome/Host-SNP_genotype_analysis/HOST_SNP.md`](docs/04_Holo-genome/Host-SNP_genotype_analysis/HOST_SNP.md), [`docs/02_Standard_MAP2B/MAP2B/MAP2B.md`](docs/02_Standard_MAP2B/MAP2B/MAP2B.md) |
| Extended module C: eDNA | Cross-domain profiling of single-tag libraries | `software/map2b/bin/MAP2B-Cross-domain.py` | [`docs/05_eDNA/Cross-domain_profiling/CROSS_DOMAIN.md`](docs/05_eDNA/Cross-domain_profiling/CROSS_DOMAIN.md) |
| Extended module C: eDNA | Serial sequencing analysis for eDNA profiling | programs in `software/serial_tag/` | [`docs/05_eDNA/Serial_sequencing/SERIAL_TAG.md`](docs/05_eDNA/Serial_sequencing/SERIAL_TAG.md) |

## Computational workflow following the protocol

Run commands from the repository root after activating `i2bRAD-M`. Replace example paths and enzyme identifiers with values appropriate for the library and reference database.

### 1. Reference database construction

#### 1.1 Construct and validate a microbial reference database

```bash
python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 1 -i /path/genome_manifest.tsv.gz -e <enzyme_id> -o /path/database

python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 3 -d /path/database
```

#### 1.2 Extend a database for holo-genome, absolute-quantification, or eDNA analysis

```bash
python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 2 -d /path/base_database -i /path/additional_genomes.tsv.gz \
  -o /path/updated_database -p 8
```

### 2. Standard MAP2B analysis

```bash
python3 software/map2b/bin/MAP2B.py \
  -i examples/02_Standard_MAP2B/MAP2B/map2b_samples.tsv \
  -e <enzyme_id> -s /absolute/path/database \
  -o results/map2b -p 4
```

i2bRAD-M supports 16 Type IIB restriction enzymes through G-score-based profiling with matching databases. A Random Forest-based workflow is additionally available for BcgI (`5`) and CjePI (`13`). The enzyme identifier must match the library preparation and database; see [`docs/02_Standard_MAP2B/MAP2B/MAP2B.md`](docs/02_Standard_MAP2B/MAP2B/MAP2B.md) for database selection and filtering options.

### 3. Optional host-referenced absolute quantification

```bash
python3 software/absolute_quantification/MAP2BAbsoluteQuantifier.py \
  -i results/map2b/Abundance.xls \
  -m examples/02_Standard_MAP2B/Host-referenced_absolute_quantification/absolute_quantification_metadata.tsv \
  -s <host_species_exactly_as_in_Species_column> \
  -o results/absolute_quantification
```

### 4. Extended module A: Multi-omics analysis

#### 4.1 Metatranscriptome analysis

```bash
python3 software/metatranscriptome/MTKrakenProfiler.py \
  -i examples/03_Multi-omics/Metatranscriptome/metatranscriptome_samples.tsv \
  -d /path/kraken2_database -o results/RNA_profile -p 8

python3 software/metatranscriptome/RNADNARatio.py \
  -d results/map2b/Abundance.xls \
  -r results/RNA_profile/RNA_species_abundance.tsv \
  -o results/RNA_DNA_ratio
```

Provide a compatible Kraken2 database with `-d`. Quality control and rRNA removal are upstream preprocessing steps and may be performed with the tools selected for the study.

#### 4.2 Epigenome analysis (MethylRAD-based)

The protocol route uses the default methylation preset, so the default enzyme option is omitted:

CpG and CHG methylation analyses are performed separately using the corresponding `-t` option; see [`docs/03_Multi-omics/Epigenome/EPIGENOME.md`](docs/03_Multi-omics/Epigenome/EPIGENOME.md).

```bash
python3 software/methylation/build_reference_reads.py \
  -r /path/reference.fa.gz

python3 software/methylation/extract_fixed_reads.py \
  -i examples/03_Multi-omics/Epigenome/epigenome_samples.tsv

python3 software/methylation/map_and_quantify_reads.py \
  -i examples/03_Multi-omics/Epigenome/epigenome_samples.tsv
```

Additional supported options are documented in [`docs/03_Multi-omics/Epigenome/EPIGENOME.md`](docs/03_Multi-omics/Epigenome/EPIGENOME.md).

### 5. Extended module B: Holo-genome analysis

#### 5.1 Host-SNP genotype analysis

Host SNP genotyping uses the stepwise Perl workflow provided in `software/host_snp/`. See [`docs/04_Holo-genome/Host-SNP_genotype_analysis/HOST_SNP.md`](docs/04_Holo-genome/Host-SNP_genotype_analysis/HOST_SNP.md) for workspace preparation, host-reference tag construction, mapping, and genotype calling.

#### 5.2 Microbiome identification

Use the integrated host-microbial reference database and follow the standard MAP2B workflow in Section 2.

### 6. Extended module C: Serial sequencing analysis for eDNA profiling

#### 6.1 Single-tag cross-domain analysis

Single-tag libraries can be analyzed directly:

```bash
python3 software/map2b/bin/MAP2B-Cross-domain.py \
  -i examples/05_eDNA/Cross-domain_profiling/cross_domain_samples.tsv \
  -e <enzyme_id> -s /absolute/path/cross_domain_database \
  -o results/cross_domain -p 4
```

#### 6.2 Five-tag serial-library preprocessing

For five-tag serial libraries, first demultiplex and recover the corresponding single-tag datasets:

```bash
python3 software/serial_tag/serial2brad.py \
  -i examples/05_eDNA/Serial_sequencing/serial_tag_samples.tsv \
  -e <enzyme_id> -o results/serial_tag -p 4
```

Prepare the MAP2B sample manifest from the recovered single-tag files and run the same cross-domain command. The serial-tag workflow supports 16 Type IIB restriction enzymes; see [`docs/05_eDNA/Serial_sequencing/SERIAL_TAG.md`](docs/05_eDNA/Serial_sequencing/SERIAL_TAG.md) and [`docs/05_eDNA/Cross-domain_profiling/CROSS_DOMAIN.md`](docs/05_eDNA/Cross-domain_profiling/CROSS_DOMAIN.md).

## Input templates and data

Input manifests and small example files are provided under `examples/`. Replace the example paths and identifiers with values for the current analysis. Large sequencing datasets and reference databases are maintained outside this source-code repository; consult the article's data-availability statement and [`data/README.md`](data/README.md).

## Verification

From the repository root on Linux:

```bash
python3 tests/run_checks.py
sha256sum -c SHA256SUMS.txt
```

These commands check the repository interfaces and file integrity.

## Citation and use

Please cite the accompanying Nature Protocols article when using this workflow. The repository contents are provided under the terms in [`LICENSE`](LICENSE); bundled third-party components remain subject to their respective terms as described in [`environment/THIRD_PARTY.md`](environment/THIRD_PARTY.md).
