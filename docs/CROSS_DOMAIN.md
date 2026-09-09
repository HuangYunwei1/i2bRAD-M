# Cross-domain eDNA profiling

Runner: `software/map2b/bin/MAP2B-Cross-domain.py`.

Run the command from the repository root; the runner uses the accompanying scripts, configuration, and database resources under `software/map2b/`.

## Database

Start from a microbial database and use database-builder mode 2 to add biologically relevant macro-organism genomes. See [`DATABASE_BUILDER.md`](DATABASE_BUILDER.md). Record database composition/genome versions and validate the expanded database.

## Run

```bash
python3 software/map2b/bin/MAP2B-Cross-domain.py \
  -i examples/manifests/map2b_samples.tsv \
  -e <enzyme_id> -s /absolute/path/cross_domain_database \
  -o results/cross_domain -p 4
```

The cross-domain workflow follows the same enzyme, database-selection, and filtering rules as standard MAP2B. The enzyme must match the library and database; see [`INPUT_OUTPUT.md`](INPUT_OUTPUT.md).

## Outputs

By default:

- `Abundance.micro.xls`: Bacteria, Archaea, Fungi
- `Abundance.macro.xls`: remaining represented kingdoms

Add `-b` to retain the total `Abundance.xls`. An existing total table can also be split separately:

```bash
python3 software/map2b/scripts/split_abundance.py \
  -i results/map2b/Abundance.xls -o results/cross_domain_split --keep-total
```

Macro-organism eDNA values depend on database scope, genome representation, DNA shedding, extraction, and library efficiency. Interpret them as comparative detection signals rather than direct biomass unless separately calibrated.
