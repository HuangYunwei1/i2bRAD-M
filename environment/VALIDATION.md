# Environment validation summary

## Clean installation

The exact lock file `env/i2bRAD-M-linux-64.validated.lock.txt` was used to create a new Conda environment named `i2bRAD-M`. SOAP2 2.19 was then installed into the same environment, and `check_environment.sh` completed successfully.

The check confirmed that `python`, `perl`, `soap`, `2bwt-builder`, `bowtie`, `bowtie-build`, and `pear` resolved from the active environment; required Python imports and Perl `Parallel::ForkManager` loaded successfully; expected tool versions were reported; and required shared libraries were available.

## Functional smoke tests

1. **SOAP2 2.19:** `2bwt-builder` constructed a synthetic reference index and SOAP2 aligned an exact synthetic read.
2. **Bowtie 1.3.0:** `bowtie-build` constructed a synthetic reference index and Bowtie produced the expected alignment.
3. **PEAR 0.9.6:** a synthetic paired-end test merged 20 of 20 read pairs.
4. **Shared libraries:** no missing libraries were reported for the dynamically linked executables checked during validation.

`check_environment.sh` is the routine post-installation check; the functional smoke tests above document the validated package set.

## Package notes

The locked Conda artifact for pandas is `pandas-2.3.1-py39h1b6b32d_0`, while `pandas.__version__` reports `2.3.2` in the validated environment. The artifact identity in the exact lock is used for reproduction.

scikit-learn 0.24.1 may emit a `pkg_resources` deprecation warning with the validated Setuptools version. This warning does not indicate an installation failure.
