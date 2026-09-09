#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="i2bRAD-M"
LOCK_FILE="$ROOT_DIR/env/i2bRAD-M-linux-64.validated.lock.txt"
FORCE=0

usage() {
  cat <<USAGE
Usage: bash install.sh [--force]

Creates the Linux x86_64 Conda environment named i2bRAD-M from the explicit
lock file, installs SOAP2 2.19, and runs the environment checker.

  --force   remove an existing i2bRAD-M environment before reinstalling
USAGE
}

case "${1:-}" in
  "") ;;
  --force) FORCE=1 ;;
  -h|--help) usage; exit 0 ;;
  *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
esac

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is not available in PATH." >&2
  exit 1
fi

if [[ "$(uname -s)" != "Linux" || "$(uname -m)" != "x86_64" ]]; then
  echo "ERROR: this package lock targets Linux x86_64 only." >&2
  exit 1
fi

if [[ ! -r "$LOCK_FILE" ]]; then
  echo "ERROR: missing lock file: $LOCK_FILE" >&2
  exit 1
fi

# Host utilities required by the i2bRAD-M workflows.
for x in gzip wget sha256sum ldd; do
  if ! command -v "$x" >/dev/null 2>&1; then
    echo "ERROR: required host utility not found: $x" >&2
    exit 1
  fi
done

# Check the minimum supported glibc version without requiring Python.
if command -v getconf >/dev/null 2>&1; then
  libc_line="$(getconf GNU_LIBC_VERSION 2>/dev/null || true)"
  libc_ver="${libc_line#glibc }"
  if [[ "$libc_line" == glibc\ * ]]; then
    IFS=. read -r libc_major libc_minor _ <<<"$libc_ver"
    if (( libc_major < 2 || (libc_major == 2 && libc_minor < 17) )); then
      echo "ERROR: glibc >=2.17 is required; detected $libc_ver" >&2
      exit 1
    fi
  fi
fi

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

exists=0
while IFS= read -r prefix; do
  [[ "$(basename "$prefix")" == "$ENV_NAME" ]] && exists=1
 done < <(conda env list | awk 'NF && $1 !~ /^#/ {print $NF}')

if (( exists )); then
  if (( FORCE )); then
    echo "Removing existing environment: $ENV_NAME"
    conda env remove -n "$ENV_NAME" -y
  else
    echo "ERROR: Conda environment '$ENV_NAME' already exists." >&2
    echo "Run 'conda activate $ENV_NAME && ./check_environment.sh' to check it," >&2
    echo "or reinstall explicitly with: bash install.sh --force" >&2
    exit 2
  fi
fi

echo "===== CREATE $ENV_NAME FROM PACKAGE LOCK ====="
conda create -n "$ENV_NAME" --file "$LOCK_FILE" -y

# Keep the environment isolated from packages installed under ~/.local.
echo "===== ISOLATE PYTHON PACKAGES ====="
conda env config vars set -n "$ENV_NAME" PYTHONNOUSERSITE=1

echo "===== ACTIVATE ====="
conda activate "$ENV_NAME"
echo "CONDA_PREFIX=$CONDA_PREFIX"

echo "===== INSTALL SOAP2 2.19 ====="
bash "$ROOT_DIR/scripts/install_soap2.sh"
hash -r

echo "===== CHECK ENVIRONMENT ====="
bash "$ROOT_DIR/check_environment.sh"

echo
echo "============================================"
echo " i2bRAD-M installation PASSED"
echo " Activate later with: conda activate i2bRAD-M"
echo "============================================"
