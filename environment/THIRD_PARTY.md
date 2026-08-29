# Third-party software notices

Third-party programs and packages used by i2bRAD-M remain subject to the terms of their respective authors and distributors. Inclusion in this repository or installation through the environment package does not relicense those components.

## SOAPaligner / SOAP2 2.19

The host-SNP workflow uses SOAP2 2.19. The Linux x86_64 executables `soap` and `2bwt-builder` are stored in `third_party/soap2/`, verified against `third_party/soap2/SHA256SUMS.txt`, and installed into the active `i2bRAD-M` environment by `scripts/install_soap2.sh`. SOAP2 source history is maintained by the BGI SOAP2 project.

## PEAR 0.9.6

PEAR is installed from the exact Conda package URL recorded in the environment lock. No PEAR executable is stored directly in this repository.

## Bowtie 1.3.0 and supporting libraries

Bowtie and its runtime dependencies are installed from the package URLs recorded in the exact Conda lock. Their package files are not stored directly in this repository.

## Conda packages

All packages retrieved by Conda retain their upstream licenses and notices. Package names, versions, builds, channels, and exact artifact URLs are recorded in `env/i2bRAD-M-linux-64.validated.lock.txt`.
