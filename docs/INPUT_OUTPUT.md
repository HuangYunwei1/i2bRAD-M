# Standard MAP2B input and output

Exact optional arguments:

```bash
python3 software/map2b/bin/MAP2B.py --help
```

## Sample manifest

MAP2B accepts a tab-separated file with no required header. Lines beginning with `#` are ignored.

```text
# single-end
sample01<TAB>/absolute/path/sample01.fastq.gz
# paired-end
sample02<TAB>/absolute/path/sample02_R1.fastq.gz<TAB>/absolute/path/sample02_R2.fastq.gz
```

Use unique sample IDs and preferably absolute paths. A template is at `examples/manifests/map2b_samples.tsv`.

## Required consistency

The same Type IIB enzyme must be used for wet-lab library construction, database construction, and MAP2B `-e`. The i2bRAD-M MAP2B interface supports enzyme IDs `1`-`16` (default `13`, CjePI), and all 16 enzymes can be analyzed using G-score filtering with `-g` and a matching database. For BcgI (`5`) and CjePI (`13`), i2bRAD-M additionally provides a Random Forest-based filtering workflow. To facilitate standard analyses, pre-built database download lists are provided for BsaXI (`3`), BcgI (`5`), and CjePI (`13`); for the other enzymes, provide a matching custom database directory with `-s`.

## Standard command

```bash
python3 software/map2b/bin/MAP2B.py \
  -i /path/sample.list -e <enzyme_id> \
  -s /absolute/path/database -o results/map2b -p 4
```

For low-biomass samples, the manuscript recommends evaluating `-g 5`. Record the selected threshold with the analysis.

## Principal output

`Abundance.xls` is the species-level relative-abundance table. Despite its extension, it is tab-separated text. Taxonomy columns precede sample columns; the `Species` column is used by the absolute quantifier. Retain the full MAP2B result directory for reproducibility and downstream analysis.

## Database location

A full biological database is not bundled. Keep generated/downloaded databases outside Git (recommended) or locally under `software/map2b/database/`, which is ignored except for its README.

## Run record

Record repository commit/version, environment version, database version/checksum, exact command, enzyme, sample-manifest checksum, and non-default thresholds.
