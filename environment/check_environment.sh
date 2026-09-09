#!/usr/bin/env bash
set -uo pipefail

EXPECTED_ENV="i2bRAD-M"
fail=0
ok()  { printf '[OK]   %s\n' "$*"; }
bad() { printf '[FAIL] %s\n' "$*" >&2; fail=1; }

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "ERROR: activate i2bRAD-M first: conda activate i2bRAD-M" >&2
  exit 1
fi

if [[ "$(basename "$CONDA_PREFIX")" == "$EXPECTED_ENV" ]]; then
  ok "active environment: $CONDA_PREFIX"
else
  bad "expected environment '$EXPECTED_ENV'; active prefix is $CONDA_PREFIX"
fi

if [[ "${PYTHONNOUSERSITE:-}" == "1" ]]; then
  ok "Python user-site packages are isolated"
else
  bad "PYTHONNOUSERSITE=1 is not active; reactivate the environment or reinstall it with install.sh"
fi

# Programs that must resolve from the active Conda environment.
for x in python perl soap 2bwt-builder bowtie bowtie-build pear kraken2 kraken2-build kraken2-inspect; do
  p="$(command -v "$x" 2>/dev/null || true)"
  if [[ -z "$p" ]]; then
    bad "$x not found"
  elif [[ "$p" == "$CONDA_PREFIX"/* ]]; then
    ok "$x -> $p"
  else
    bad "$x is outside the active environment: $p"
  fi
done

# Project scripts also use these standard Linux utilities, so availability
# (rather than prefix ownership) is the relevant release check.
for x in gzip wget; do
  p="$(command -v "$x" 2>/dev/null || true)"
  if [[ -n "$p" ]]; then
    ok "host utility $x -> $p"
  else
    bad "required host utility not found: $x"
  fi
done

# Check Python imports and runtime versions. Package identities are checked
# from conda-meta below; marisa-trie does not expose a Python __version__.
if python - <<'EOF_PY'
import sys, warnings
warnings.filterwarnings(
    "ignore",
    message=r"pkg_resources is deprecated as an API.*",
    category=UserWarning,
)
import joblib, marisa_trie, numpy, pandas, pysam, sklearn

expected = {
    "python": "3.9.23",
    "numpy": "1.26.4",
    "sklearn": "0.24.1",
    "joblib": "1.5.1",
    "pysam": "0.23.3",
    "pandas": "2.3.1",
}
actual = {
    "python": ".".join(map(str, sys.version_info[:3])),
    "numpy": numpy.__version__,
    "sklearn": sklearn.__version__,
    "joblib": joblib.__version__,
    "pysam": pysam.__version__,
    "pandas": pandas.__version__,
}
for key, want in expected.items():
    got = actual[key]
    if got != want:
        raise SystemExit(f"{key}: expected {want}, got {got}")

print("[OK]   Python imports and runtime versions")
for key in ("python", "numpy", "sklearn", "joblib", "pysam", "pandas"):
    print(f"       {key}: {actual[key]}")
print("       marisa_trie: import OK")
EOF_PY
then
  :
else
  bad "Python import/version check"
fi

# Conda artifacts recorded in the package lock.
expected_meta=(
  'python-3.9.23-hc30ae73_0_cpython.json'
  'perl-5.26.2-h36c2ea0_1008.json'
  'perl-parallel-forkmanager-2.02-pl526_0.json'
  'joblib-1.5.1-pyhd8ed1ab_0.json'
  'marisa-trie-0.7.7-py39h5a03fae_2.json'
  'numpy-1.26.4-py39h474f0d3_0.json'
  'pandas-2.3.1-py39h1b6b32d_0.json'
  'pysam-0.23.3-py39hbe927e6_3.json'
  'scikit-learn-0.24.1-py39h4dfa638_0.json'
  'scipy-1.13.1-py39haf93ffa_0.json'
  'setuptools-80.9.0-pyhff2d567_0.json'
  'bowtie-1.3.0-py39h176da8b_2.json'
  'tbb-2020.2-h4bd325d_4.json'
  'pear-0.9.6-h9d449c0_10.json'
  'libzlib-1.3.2-h25fd6f3_3.json'
  'kraken2-2.1.2-pl5262h7d875b9_0.json'
)
meta_missing=0
for f in "${expected_meta[@]}"; do
  if [[ ! -f "$CONDA_PREFIX/conda-meta/$f" ]]; then
    bad "locked Conda artifact metadata missing: $f"
    meta_missing=1
  fi
done
if (( ! meta_missing )); then
  ok "locked Conda artifact set"
  ok "pandas Conda artifact: 2.3.1 py39h1b6b32d_0"
fi

perl_ver="$(perl -e 'printf "%vd", $^V' 2>/dev/null || true)"
if [[ "$perl_ver" == "5.26.2" ]]; then
  ok "Perl version 5.26.2"
else
  bad "Perl version check (got: ${perl_ver:-unknown})"
fi

pfm_ver="$(perl -MParallel::ForkManager -e 'print $Parallel::ForkManager::VERSION' 2>/dev/null || true)"
if [[ "$pfm_ver" == "2.02" ]]; then
  ok "Parallel::ForkManager 2.02"
else
  bad "Parallel::ForkManager check (got: ${pfm_ver:-unknown})"
fi

soap_out="$(soap 2>&1 || true)"
if grep -Fq 'Version: 2.19' <<<"$soap_out"; then
  ok "SOAP2 version 2.19"
else
  bad "SOAP2 version check"
fi

bowtie_out="$(bowtie --version 2>&1 || true)"
if grep -Fq 'version 1.3.0' <<<"$bowtie_out"; then
  ok "Bowtie version 1.3.0"
else
  bad "Bowtie version check"
fi

pear_out="$(pear -h 2>&1 || true)"
if grep -Fq 'PEAR v0.9.6' <<<"$pear_out"; then
  ok "PEAR version 0.9.6"
else
  bad "PEAR version check"
fi

kraken_out="$(kraken2 --version 2>&1 || true)"
if grep -Fq 'Kraken version 2.1.2' <<<"$kraken_out"; then
  ok "Kraken2 version 2.1.2"
else
  bad "Kraken2 version check"
fi

# Shared-library sanity checks. SOAP2 2.19 may legitimately be non-dynamic.
for x in soap 2bwt-builder bowtie-align-s bowtie-build-s pear; do
  p="$CONDA_PREFIX/bin/$x"
  [[ -x "$p" ]] || { bad "$x executable missing at $p"; continue; }
  ldd_out="$(ldd "$p" 2>&1 || true)"
  if grep -Fq 'not found' <<<"$ldd_out"; then
    bad "$x has missing shared libraries"
  elif grep -Eq 'not a dynamic executable|statically linked' <<<"$ldd_out"; then
    ok "$x shared-library check (static/non-dynamic)"
  else
    ok "$x shared-library check"
  fi
done

if (( fail )); then
  echo "Environment check FAILED." >&2
  exit 1
fi

echo "Environment check PASSED."
