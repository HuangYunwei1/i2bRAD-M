# Metatranscriptome analysis

This module profiles species-level signals from paired metatranscriptomic reads with Kraken2 and optionally integrates the RNA relative-abundance table with the DNA-based MAP2B `Abundance.xls` table.

## Requirements and input stage

`MTKrakenProfiler.py` requires Python 3.9 or later, a Kraken2 executable available in `PATH`, and a compatible Kraken2 database. The Kraken2 database is maintained outside this repository and is supplied with `-d`.

The script expects paired reads that have already undergone the preprocessing appropriate for the study, such as quality control and rRNA removal. Tools such as fastp and SortMeRNA may be used for these upstream steps, but they are not called by the programs in this module.

## RNA sample manifest

Use a tab-separated file with the required header:

```text
sample_id  read1  read2
```

Paths may be absolute or relative to the manifest. Sample identifiers must contain only letters, numbers, `.`, `_`, or `-`. See `examples/manifests/metatranscriptome_samples.tsv`.

## Species-level RNA profiling

From the repository root:

```bash
python3 software/metatranscriptome/MTKrakenProfiler.py \
  -i examples/manifests/metatranscriptome_samples.tsv \
  -d /path/to/kraken2_database \
  -o results/RNA_profile -p 8
```

Principal outputs include:

- `RNA_species_counts.tsv`: species-level fragment-count matrix;
- `RNA_species_abundance.tsv`: species-level relative-abundance matrix;
- `sample_qc.tsv`: per-sample classification and quality-control summary;
- `kraken_reports/` and `kraken_output/`: Kraken2 report and classification files;
- `logs/<sample>.log`: per-sample Kraken2 run log;
- `config.json`: analysis parameters and software information.

## RNA and DNA integration

Use the MAP2B species abundance table and the RNA abundance table as input files:

```bash
python3 software/metatranscriptome/RNADNARatio.py \
  -d results/map2b/Abundance.xls \
  -r results/RNA_profile/RNA_species_abundance.tsv \
  -o results/RNA_DNA_ratio \
  -a 0.001 -b 0.001
```

`-a` and `-b` specify the minimum DNA and RNA relative abundances. Setting either threshold to `0` disables that threshold.

The principal output is `RNA_DNA_ratio.tsv`, with one row per matched sample-species combination and the following fields:

```text
Sample  Species  DNA_relative_abundance  RNA_relative_abundance  RNA_DNA_ratio  log2_ratio  Status
```

Supporting outputs include `sample_summary.tsv`, `unmatched_samples.tsv`, `unmatched_species.tsv`, `run.log`, and `config.json`, which record sample matching, species matching, filtering summaries, and the exact configuration.

## Example tables

Small example DNA and RNA abundance tables are provided under `examples/metatranscriptome/` for checking table structure and the integration command.
