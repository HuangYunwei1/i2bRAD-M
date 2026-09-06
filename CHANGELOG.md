# Changelog

## Unreleased — 2026-09-05

- Added a serial-tag preprocessing module that supports 16 Type IIB restriction enzymes and produces five position-resolved fixed-length tags for downstream profiling.
- Added metatranscriptome species profiling with Kraken2 and optional RNA/DNA relative-abundance integration.
- Added reviewer-facing workflow documentation and input examples for the serial-tag and metatranscriptome modules.
- Extended the MAP2B and MAP2B-Cross-domain single-enzyme workflows to accept all 16 supported Type IIB enzymes.
- Added the complete 16-enzyme digestion patterns and quantitative-database tag lengths; CspCI uses the 33-bp definition, and BplI/FalI/AlfI are scanned in one direction to avoid duplicate extraction.
- Retained random-forest use for BsaXI, BcgI and CjePI; the other enzymes require a matching custom database and G-score filtering.
- Updated MAP2BDatabaseBuilder to v0.3.2.
- Reworked global unique-database reconstruction to use streaming and GNU external sort, substantially reducing peak memory use on large databases; added optional `-p/--processes` sort parallelism.
- Added safe update support for official or legacy MAP2B databases without `genome_id_map.tsv`, while preserving existing internal IDs and master shards.
- Harmonized added genomes to the base-database lineage when an existing species name has a different taxonomy annotation.
- Extended validation to accept official 11-column statistics files as well as builder-generated 10-column files.
- Recorded the successful GTDB-plus-five-genomes update benchmark (30 min 11 s; approximately 15.4 GiB peak memory).

## 0.1.0 — 2026-08-20

- Updated MAP2BDatabaseBuilder to v0.3.1, allowing tag copy numbers above 9999 while preserving the fixed 8-digit genome-ID prefix and downstream database formats.
- Added a unified Linux Conda environment named `i2bRAD-M` with an exact package lock and post-installation validation.
- Added MAP2B microbial profiling and cross-domain analysis workflows.
- Added database construction, absolute quantification, host-SNP, and methylation modules.
- Added workflow documentation, tab-delimited input templates, integrity checksums, and repository checks.
- Standardized repository paths and documentation for use with the accompanying Nature Protocols article.
