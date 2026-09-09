# Environment validation summary

## Clean installation

The exact lock file `env/i2bRAD-M-linux-64.validated.lock.txt` was used to create a new Conda environment named `i2bRAD-M`. SOAP2 2.19 was then installed into the same environment, and `check_environment.sh` completed successfully.

The check confirmed that `python`, `perl`, `soap`, `2bwt-builder`, `bowtie`, `bowtie-build`, `pear`, `kraken2`, `kraken2-build`, and `kraken2-inspect` resolved from the active environment; required Python imports and Perl `Parallel::ForkManager` loaded successfully; expected tool versions were reported; and required shared libraries were available.

## Functional smoke tests

1. **SOAP2 2.19:** `2bwt-builder` constructed a synthetic reference index and SOAP2 aligned an exact synthetic read.
2. **Bowtie 1.3.0:** `bowtie-build` constructed a synthetic reference index and Bowtie produced the expected alignment.
3. **PEAR 0.9.6:** a synthetic paired-end test merged 20 of 20 read pairs.
4. **Kraken2 2.1.2:** `MTKrakenProfiler.py` successfully classified a 100,000-fragment paired-end test dataset with an external Kraken2 database and generated the species-count, relative-abundance, sample-QC, report, classification, log, and configuration outputs.
5. **Cross-module compatibility:** all host-SNP Perl programs passed syntax checks, and the MAP2B and database-builder command-line interfaces remained operational in the combined environment.
6. **Shared libraries:** no missing libraries were reported for the dynamically linked executables checked during validation.

`check_environment.sh` is the routine post-installation check; the functional smoke tests above document the validated package set.

## Package record

The locked Kraken2 artifact is `kraken2-2.1.2-pl5262h7d875b9_0`, selected to retain the validated Perl 5.26.2 runtime used by the host-SNP workflow. The locked pandas artifact is `pandas-2.3.1-py39h1b6b32d_0`, and `pandas.__version__` reports `2.3.1` in the clean validated environment.
