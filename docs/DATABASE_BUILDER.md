# MAP2B database builder

Program: `software/database_builder/MAP2BDatabaseBuilder.py` (v0.3.2). It builds, updates, or validates species-level MAP2B-compatible tag databases.

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

The master-database key contains an 8-digit internal genome ID followed by the tag copy number, zero-padded to at least four digits. Version 0.3.2 supports copy numbers above 9999 without truncation; for example, copy number 13085 is stored as `0000000013085`.

## Validate (`-m 3`)

```bash
python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 3 -d /path/database
```

With `-e` omitted, validation detects all represented enzymes. It checks taxonomy IDs, shard placement, tag lengths and composition, copy-number keys, unique-tag consistency, statistics, and the database manifest. Both the builder's 10-column statistics format and the official MAP2B 11-column format are accepted.

## Update (`-m 2`)

```bash
python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 2 -d /path/base_database \
  -i /path/additional_genomes.tsv.gz \
  -o /path/updated_database \
  -p 8
```

Mode 2 normally detects the enzyme and shard size from the base database. New source IDs must not collide with existing records. `-p/--processes` controls GNU sort parallelism during global unique-database reconstruction and defaults to 1. Prefer a new output directory over `--in-place`.

Version 0.3.2 reconstructs the global unique database as a streamed external-sort workflow rather than retaining the complete tag map in a Python dictionary. Temporary files are removed automatically after normal completion or an ordinary error. The sort buffer is selected automatically (approximately 5% of visible memory, capped at 4 GiB); `TMPDIR` may be set to a large, fast scratch location when available.

Official or legacy MAP2B databases without `genome_id_map.tsv` can be updated when their taxonomy contains unique 8-digit internal IDs. The existing IDs and master shards remain unchanged, and the updated output receives a cumulative mapping file with synthetic identifiers for legacy records. If an added genome belongs to a species already present but uses a different kingdom-to-species lineage, the base-database lineage is treated as authoritative; the added source ID, strain information, and genome path are retained. Existing mapping files, when present, remain subject to strict consistency checks.

## Tested update performance

A BsaXI update of the 259,389-genome GTDB database with five additional microbial genomes was completed on September 2, 2026 using `-p 20`. The update processed 316,888,595 genome-tag records and completed in 30 min 11 s. Peak resident memory was 16,171,960 kB (approximately 15.4 GiB), no swap was used, and the completed 259,394-genome database occupied approximately 2.5 GB. Runtime depends on database size, storage performance, available memory, and server load.

## Cross-domain database

Start from the microbial database and use mode 2 to add biologically relevant host, animal, plant, dietary, or other genomes. Validate the completed database before analysis. The enzyme must remain identical across the library, base database, added genomes, and downstream profiling.

## Files to retain with an analysis

Retain the complete database directory, builder version and checksum, input manifest, validation output, enzyme setting, and database checksums. The checksum of the repository copy of the builder is recorded in `software/database_builder/SHA256SUMS.txt`.
