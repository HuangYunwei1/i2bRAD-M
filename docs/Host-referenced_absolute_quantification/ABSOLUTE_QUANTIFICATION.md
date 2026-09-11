# Host-referenced absolute quantification

Program: `software/absolute_quantification/MAP2BAbsoluteQuantifier.py` (v0.2.0).

## Prerequisites

1. Use a database containing both microbes and the selected host.
2. Run MAP2B against that database.
3. Measure host-cell concentration for every sample.
4. Confirm the exact host string in the `Species` column of `Abundance.xls`.

## Metadata

Concentration only:

```text
sample_id	ch
sample01	2000000
```

`ch` is host cells/mL. Each sample ID must exactly match a sample column in `Abundance.xls`.

With library cell equivalents (`-c`):

```text
sample_id	ch	ve	vd	vl
sample01	2000000	0.5	100	10
```

`ve` is original sample volume used for extraction (mL), `vd` is final DNA elution volume, and `vl` is DNA volume used for library construction. Keep `vd`/`vl` units consistent.

## Commands

```bash
python3 software/absolute_quantification/MAP2BAbsoluteQuantifier.py \
  -i results/map2b/Abundance.xls \
  -m examples/Host-referenced_absolute_quantification/absolute_quantification_metadata.tsv \
  -s <host_name_exactly_as_in_Species> \
  -o results/absolute_quantification
```

Add `-c` and use the five-column template to calculate library equivalents.

## Principal outputs

| Output file | Format | Description |
|---|---|---|
| `AbsoluteConcentration.xls` | Tab-separated text | Host-referenced concentration estimates in cells/mL |
| `LibraryCellEquivalent.xls` | Tab-separated text | Library cell-equivalent estimates; generated only with `-c` |

Both tables retain the MAP2B taxonomy columns followed by sample columns.

The result is host-referenced, not an independent direct count. Report host measurement, database composition, naming, and assumptions with the output.
