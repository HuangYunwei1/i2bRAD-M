# MAP2B database builder

Program: `software/database_builder/MAP2BDatabaseBuilder.py` (v0.3.2). It builds, updates, or validates species-level MAP2B-compatible tag databases for 16 Type IIB restriction enzymes.

## Input manifest

Use a tab-separated plain-text or gzip-compressed table with at least nine columns. Column 8 is `species`; the final column is the genome FASTA path.

```text
#source_id	kingdom	phylum	class	order	family	genus	species	genome_path
GCF_000005845.2	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Escherichia	Escherichia_coli	/path/GCF_000005845.2.fna.gz
```

`source_id` must be unique; each FASTA must exist and contain a header. Taxonomy must be consistent across assemblies assigned to the same species. See `examples/manifests/database_genomes.tsv`.

## Enzyme selection

Specify one enzyme identifier from `1` to `16` with `-e`, or use `-e 17` (`AllEnzyme`) to process all 16 enzymes separately. The enzyme must match the wet-lab library and downstream MAP2B analysis.

## Build (`-m 1`)

```bash
python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 1 -i /path/genome_manifest.tsv.gz \
  -e <enzyme_id> -o /path/new_database
```

Mode 1 defaults to CjePI (`-e 13`) and a shard size of 30000. Specify the enzyme explicitly for reproducible analysis and use a new output directory. The optional `--force` setting is intended for a deliberate rebuild of recognized database files.

## Validate (`-m 3`)

```bash
python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 3 -d /path/database
```

With `-e` omitted, validation detects all represented enzymes. It checks the taxonomy, tag database, statistics, and database metadata. Both Builder-generated and official MAP2B database layouts are supported.

## Update (`-m 2`)

```bash
python3 software/database_builder/MAP2BDatabaseBuilder.py \
  -m 2 -d /path/base_database \
  -i /path/additional_genomes.tsv.gz \
  -o /path/updated_database \
  -p 8
```

Mode 2 normally detects the enzyme and shard size from the base database. Existing MAP2B databases can serve as the base database when their taxonomy table contains valid unique identifiers. Use unique source IDs for the added genomes. `-p/--processes` controls parallel sorting during global unique-database reconstruction and defaults to 1. A new output directory is recommended, and temporary working files are managed automatically.

## Cross-domain database

Start from a microbial database and use mode 2 to add biologically relevant host, animal, plant, dietary, or other genomes. Validate the completed database before analysis. The enzyme must remain identical across the library, base database, added genomes, and downstream profiling.

## Files to retain

Retain the complete database directory, Builder version, input manifest, enzyme setting, validation record, and database checksums. The checksum of the repository copy of the Builder is recorded in `software/database_builder/SHA256SUMS.txt`.
