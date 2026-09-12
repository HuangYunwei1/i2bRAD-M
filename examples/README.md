# Input templates and example tables

This directory follows the analysis modules in the accompanying protocol:

| Module | Files |
|---|---|
| `01_Reference_database_construction/` | Genome manifest for MAP2B database construction |
| `02_Standard_MAP2B/MAP2B/` | Standard MAP2B sample manifest |
| `02_Standard_MAP2B/Host-referenced_absolute_quantification/` | Host-concentration metadata templates, with optional library volumes |
| `03_Multi-omics/Metatranscriptome/` | Metatranscriptome manifest and compact RNA/DNA integration tables |
| `03_Multi-omics/Epigenome/` | Epigenome analysis manifest for MethylRAD data |
| `04_Holo-genome/Host-SNP_genotype_analysis/` | Host-SNP working-directory example |
| `05_eDNA/Cross-domain_profiling/` | Cross-domain sample manifest |
| `05_eDNA/Serial_sequencing/` | Five-tag serial-library manifest |

Copy the relevant template to a working directory, replace the example values and paths, and preserve tab delimiters. Store large input and output files outside the Git repository.

For a compact executable workflow demonstration with selected expected outputs, download the [i2bRAD-M Demo from Zenodo](https://doi.org/10.5281/zenodo.22719104).
