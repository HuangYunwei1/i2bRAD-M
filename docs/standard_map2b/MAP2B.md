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

Use unique sample IDs and preferably absolute paths. A template is at `examples/standard_map2b/map2b_samples.tsv`.

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

| Output file | Format | Description |
|---|---|---|
| `Abundance.xls` | Tab-separated text | Species-level relative-abundance matrix |

The table begins with `#Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, and `Species`, followed by one abundance column per sample. The `Species` column is used by downstream modules such as absolute quantification. Retain the full MAP2B result directory for reproducibility and downstream analysis.

## Database location

A full biological database is not bundled. Store prepared MAP2B databases in a suitable local directory and provide the database path with `-s`.

## Run record

Record the software and database versions, enzyme, exact command, and non-default thresholds used for the analysis.
