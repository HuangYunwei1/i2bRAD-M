# Changelog

## Unreleased — 2026-09-05

- Added a serial-tag preprocessing module that supports 16 Type IIB restriction enzymes and produces five position-resolved fixed-length tags for downstream profiling.
- Added metatranscriptome species profiling with Kraken2 and optional RNA/DNA relative-abundance integration.
- Added workflow documentation and input examples for the serial-tag and metatranscriptome modules.
- Extended the MAP2B and MAP2B-Cross-domain workflows to support 16 Type IIB restriction enzymes with matching databases and G-score filtering.
- Added complete 16-enzyme digestion definitions and quantitative-database tag lengths.
- Added a Random Forest-based filtering workflow for BcgI and CjePI.
- Updated MAP2BDatabaseBuilder to v0.3.2 with memory-efficient database updates, optional parallel sorting, compatibility with existing MAP2B databases, and expanded validation.

## 0.1.0 — 2026-08-20

- Updated MAP2BDatabaseBuilder to v0.3.1, allowing tag copy numbers above 9999 while preserving the fixed 8-digit genome-ID prefix and downstream database formats.
- Added a unified Linux Conda environment named `i2bRAD-M` with an exact package lock and post-installation validation.
- Added MAP2B microbial profiling and cross-domain analysis workflows.
- Added database construction, absolute quantification, host-SNP, and methylation modules.
- Added workflow documentation, tab-delimited input templates, integrity checksums, and repository checks.
- Standardized repository paths and documentation for use with the accompanying Nature Protocols article.
