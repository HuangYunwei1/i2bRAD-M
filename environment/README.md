# i2bRAD-M computational environment

This directory defines the single Conda environment used by all computational modules in this repository. The environment name is `i2bRAD-M`.

## Installation

From the repository root:

```bash
cd environment
bash install.sh
conda activate i2bRAD-M
./check_environment.sh
cd ..
```

To replace an existing environment with the same name:

```bash
cd environment
bash install.sh --force
```

## Reproducibility files

- `env/i2bRAD-M-linux-64.validated.lock.txt` is the authoritative Linux x86_64 package set used by `install.sh`.
- `environment.yml` is the readable declaration of direct dependencies.
- `scripts/install_soap2.sh` installs the included SOAP2 2.19 executables after SHA256 verification.
- `check_environment.sh` checks executable resolution, package imports, versions, and shared-library availability.

The exact lock avoids dependency re-solving and is the recommended installation route.

## Platform requirements

- Linux x86_64;
- glibc 2.17 or later;
- Conda;
- `bash`, `gzip`, `wget`, `sha256sum`, and `ldd` available on the host.

Native Windows and macOS installations have not been validated.

## Core validated tools

- Python 3.9.23
- Perl 5.26.2
- Parallel::ForkManager 2.02
- SOAP2 2.19 (`soap`, `2bwt-builder`)
- Bowtie 1.3.0
- PEAR 0.9.6
- NumPy 1.26.4
- scikit-learn 0.24.1
- joblib 1.5.1
- marisa-trie 0.7.7
- pysam 0.23.3

See [`VALIDATION.md`](VALIDATION.md) for the validation summary and [`THIRD_PARTY.md`](THIRD_PARTY.md) for third-party software notices.
