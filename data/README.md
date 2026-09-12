# Data and database availability

The datasets supporting the protocol are available through the following NCBI BioProjects:

- [`PRJNA689204`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA689204): published 2bRAD-M datasets from the ATCC MSA 1002 mock community and FFPE samples;
- [`PRJNA1033794`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1033794) and [`PRJNA1284320`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1284320): whole-metagenome sequencing data from the Gut Puzzle cohort;
- [`PRJNA1022766`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1022766): published larval *Patinopecten yessoensis* 2bRAD-M data; and
- [`PRJNA1521654`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1521654): newly generated milk 2bRAD-M, scallop MethylRAD, and eDNA 2bRAD-M datasets.

A compact workflow demonstration package with runnable inputs and selected expected outputs is available from Zenodo: [`10.5281/zenodo.22719104`](https://doi.org/10.5281/zenodo.22719104).

Large reference genomes and MAP2B databases are maintained outside this source-code repository. Pre-built MAP2B databases can be downloaded as described in the root [`README.md`](../README.md); custom databases can be prepared with the database builder. The MAP2B `-s` option should point to the database directory used for the analysis.
