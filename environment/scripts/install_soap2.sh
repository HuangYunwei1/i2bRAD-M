#!/usr/bin/env bash
set -euo pipefail

EXPECTED_ENV="i2bRAD-M"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SOAP_DIR="$ROOT_DIR/third_party/soap2"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "ERROR: no Conda environment is active." >&2
  exit 1
fi

if [[ "$(basename "$CONDA_PREFIX")" != "$EXPECTED_ENV" ]]; then
  echo "ERROR: expected active environment '$EXPECTED_ENV'; got '$CONDA_PREFIX'." >&2
  exit 1
fi

if [[ "$(uname -s)" != "Linux" || "$(uname -m)" != "x86_64" ]]; then
  echo "ERROR: validated SOAP2 2.19 binaries target Linux x86_64." >&2
  exit 1
fi

for f in soap 2bwt-builder; do
  [[ -x "$SOAP_DIR/$f" ]] || {
    echo "ERROR: missing executable: $SOAP_DIR/$f" >&2
    exit 1
  }
done

printf '%s  %s\n' \
  '32ea145ca2ee9efd93312b1de89578a3a7fd524355b372c018c91ec0e9b67194' "$SOAP_DIR/soap" \
  '4a75e5c9e0149b89102199dfbeac24d3465f3a9ac9b7119eac5b7c56acc41982' "$SOAP_DIR/2bwt-builder" \
  | sha256sum -c -

install -m 0755 "$SOAP_DIR/soap" "$CONDA_PREFIX/bin/soap"
install -m 0755 "$SOAP_DIR/2bwt-builder" "$CONDA_PREFIX/bin/2bwt-builder"

echo "SOAP2 2.19 installed into: $CONDA_PREFIX/bin"
