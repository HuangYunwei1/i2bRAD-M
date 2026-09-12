# i2bRAD-M Demo

The compact i2bRAD-M Demo is archived on Zenodo:

- DOI: [`10.5281/zenodo.22719104`](https://doi.org/10.5281/zenodo.22719104)
- File: `i2bRAD-M_Demo_v1.0.tar.gz`
- SHA256: `f7975ecfb766d55a8d58c13f3d2582392cc57999388f43206c3ffa258bb3b6fc`

## Download and preparation

Activate the `i2bRAD-M` environment and run:

```bash
wget "https://zenodo.org/records/22719104/files/i2bRAD-M_Demo_v1.0.tar.gz?download=1" \
  -O i2bRAD-M_Demo_v1.0.tar.gz
tar -xzf i2bRAD-M_Demo_v1.0.tar.gz
export I2BRADM_ROOT=/path/to/i2bRAD-M
cd i2bRAD-M_Demo_v1.0
```

Each module contains compact input files, a `run_demo.sh` entry point, and selected results under `expected_output/`.

## Demo modules

Run the examples independently in the protocol order:

| Module | Run directory | Main expected output |
|---|---|---|
| Reference database construction | `01_Reference_database_construction/` | MAP2B database metadata and taxonomy tables |
| Standard MAP2B | `02_Standard_MAP2B/MAP2B/` | `Abundance.xls`, `Coverage.xls` |
| Host-referenced absolute quantification | `02_Standard_MAP2B/Host-referenced_absolute_quantification/` | `AbsoluteConcentration.xls` |
| Metatranscriptome | `03_Multi-omics/Metatranscriptome/` | RNA species profiles and RNA/DNA ratios |
| Epigenome | `03_Multi-omics/Epigenome/` | CpG and CHG matrices and sample QC |
| Holo-genome Host-SNP | `04_Holo-genome/Host-SNP/` | `filter_of_all_codom` |
| Holo-genome MAP2B | `04_Holo-genome/MAP2B/` | `Abundance.xls`, `Coverage.xls` |
| eDNA cross-domain profiling | `05_eDNA/Cross-domain_profiling/` | combined, microbial, and macro-organism abundance tables |
| eDNA serial sequencing | `05_eDNA/Serial_sequencing/` | recovered tag-count summary |

For example:

```bash
cd 02_Standard_MAP2B/MAP2B
bash run_demo.sh
```

The scripts create a `work/` directory and do not overwrite a previous run. After all modules have been run, return to the Demo root and validate the key results:

```bash
python3 validate_outputs.py
```

Module-specific input formats and full workflow instructions remain documented under this repository's `docs/` directory.
