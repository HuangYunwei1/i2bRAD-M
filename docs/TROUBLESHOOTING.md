# Troubleshooting

## Commands or dependencies are unavailable

Activate the unified environment and run its checker:

```bash
cd environment
conda activate i2bRAD-M
bash check_environment.sh
cd ..
```

If the environment has not been installed, run `bash environment/install.sh` from the repository root. The exact lock used by the installer is the supported reproduction route.

## MAP2B cannot find scripts, configuration, or database files

Keep `MAP2B.py` and `MAP2B-Cross-domain.py` in `software/map2b/bin/`, because their companion scripts and configuration are resolved from the MAP2B directory layout. Supply an absolute path to the prepared database with `-s`.

## Database validation fails

Check the following:

- the input file is tab-delimited;
- genome paths exist and FASTA headers are valid;
- source identifiers are unique;
- taxonomy fields are consistent;
- species is in the required column;
- the genome path is in the final column;
- the enzyme and tag length match the intended library; and
- update identifiers do not collide with existing records.

Run the database-builder validation mode before profiling.

## Absolute quantification cannot match samples or host rows

Match metadata `sample_id` values exactly to the abundance-table sample columns. The value supplied with `-s` must exactly match a value in the `Species` column. Use a five-column metadata table only with the corresponding concentration/volume option described in [`ABSOLUTE_QUANTIFICATION.md`](02_Standard_MAP2B/Host-referenced_absolute_quantification/ABSOLUTE_QUANTIFICATION.md).

## Epigenome analysis stages cannot find earlier outputs

Run all stages from one working directory or supply matching reference and cleaned-read directories explicitly. Use the same methylation preset throughout a run. For the standard article workflow, the default preset is used and no explicit `-e` option is required.

## Host-SNP mapping fails

From the workflow working directory, verify that the required `proc_data/`, `ref/HQ_ref_codom`, and `soap/ref` paths exist. Confirm that `soap` and `2bwt-builder` resolve from the active `i2bRAD-M` environment and that sample names are consistent between input files and intermediate outputs.

## Path-related problems

The validated target is Linux x86_64. Use a Linux filesystem and, where possible, simple paths without spaces or shell-special characters. Use absolute paths for external databases and study data.

## Verify the repository checkout

```bash
python3 tests/run_checks.py
sha256sum -c SHA256SUMS.txt
```

If a checksum fails, restore the affected file from the tagged repository release before continuing.
