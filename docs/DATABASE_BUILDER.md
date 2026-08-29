# MAP2B database builder

Program: `software/database_builder/MAP2BDatabaseBuilder.py` (v0.3.1). It builds, updates, or validates species-level MAP2B-compatible tag databases.

## Input manifest

Use a tab-separated plain-text or gzip-compressed table with at least nine columns. Column 8 is `species`; the final column is the genome FASTA path.

```text
#source_id	kingdom	phylum	class	order	family	genus	species	genome_path
GCF_000005845.2	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Escherichia	Escherichia_coli	/path/GCF_000005845.2.fna.gz
```

`source_id` must be unique; each FASTA must exist and contain a header. Taxonomy must be consistent across assemblies assigned to the same species. See `examples/manifests/database_genomes.tsv`.

## Build (`-m 1`)

```bash
python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 1 -i /path/genome_manifest.tsv.gz \
  -e <enzyme_id> -o /path/new_database
```

Mode 1 defaults to CjePI (`-e 13`) and a shard size of 30000, but the enzyme should be specified explicitly for reproducible analysis. Use a new output directory. `--force` deliberately removes recognized database artifacts from the selected output directory and should be used only for an intentional rebuild.

The master-database key contains an 8-digit internal genome ID followed by the tag copy number, zero-padded to at least four digits. Version 0.3.1 supports copy numbers above 9999 without truncation; for example, copy number 13085 is stored as `0000000013085`.

## Validate (`-m 3`)

```bash
python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 3 -d /path/database
```

With `-e` omitted, validation detects all represented enzymes. It checks taxonomy IDs, shard placement, tag lengths and composition, copy-number keys, unique-tag consistency, statistics, and the database manifest.

## Update (`-m 2`)

```bash
python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 2 -d /path/base_database \
  -i /path/additional_genomes.tsv.gz \
  -o /path/updated_database
```

Mode 2 normally detects the enzyme and shard size from the base database. New source IDs must not collide with existing records. Prefer a new output directory over `--in-place`.

## Cross-domain database

Start from the microbial database and use mode 2 to add biologically relevant host, animal, plant, dietary, or other genomes. Validate the completed database before analysis. The enzyme must remain identical across the library, base database, added genomes, and downstream profiling.

## Files to retain with an analysis

Retain the complete database directory, builder version and checksum, input manifest, validation output, enzyme setting, and database checksums. The checksum of the repository copy of the builder is recorded in `software/database_builder/SHA256SUMS.txt`.
