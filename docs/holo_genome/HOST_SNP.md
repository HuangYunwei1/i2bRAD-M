# Host SNP workflow

Programs in `software/host_snp/` implement tag extraction, reference preparation, SOAP2 mapping, and genotype calling. The manuscript example uses BsaXI (`-e 3`). Use the same enzyme ID for sample and reference.

| Program | Purpose |
|---|---|
| `sample_extraction.pl` | extract sample Type IIB tags |
| `Extract_cut_site.pl` | extract reference tags and padded SOAP2 reference |
| `mk_HQ_ref_codom.pl` | convert `ref_tag` to `HQ_ref_codom` |
| `reads_map.pl` | map tags with SOAP2 |
| `codom_calling.pl` | call/filter codominant genotypes |

The enzyme table is `software/host_snp/ENZYME_TABLE.tsv`.

## Environment

```bash
conda activate i2bRAD-M
```

If needed, run `environment/check_environment.sh` from the repository to verify the installed tools.

## Suggested layout

```text
host_snp_project/
├── rawdata/
├── genome/
├── sample_tags/
├── ref_build/
└── work/
    ├── proc_data/
    ├── ref/
    └── soap/
```

```bash
SCRIPTS=/absolute/path/to/i2bRAD-M/software/host_snp
```

## 1. Extract sample tags

```bash
perl "$SCRIPTS/sample_extraction.pl" \
  -i rawdata/sample01.fastq.gz -e 3 \
  -od sample_tags -op sample01 -qc yes -fm fa
```

## 2. Prepare reference

```bash
perl "$SCRIPTS/Extract_cut_site.pl" \
  -r genome/reference.fa.gz -e 3 -o ref_build
perl "$SCRIPTS/mk_HQ_ref_codom.pl" \
  -i ref_build/ref_tag -o ref_build/HQ_ref_codom
```

## 3. Arrange and map

```bash
mkdir -p work/proc_data work/ref work/soap
cp sample_tags/sample01.BsaXI.fa work/proc_data/
cp ref_build/ref work/soap/
cp ref_build/HQ_ref_codom work/ref/
cd work
perl "$SCRIPTS/reads_map.pl" -e 3
find reads_mapping -maxdepth 1 -type f -printf '%f\n' | LC_ALL=C sort > sample_order.txt
```

Confirm the exact sample filename before copying. `reads_map.pl` expects this relative working layout.

## 4. Call genotypes

```bash
perl "$SCRIPTS/codom_calling.pl" -c 0.3 -m 0.05
```

`-c` is call-rate (0–1); `-m` is minor-allele frequency (0–0.5). Preserve `sample_order.txt`; it defines sample-column order.

## Principal outputs

| Output | Description |
|---|---|
| `sample_order.txt` | Sample-column order used in genotype tables |
| `reads_mapping/<sample>.Bs` | Per-sample mapping result |
| `genotype/all_codom` | All called codominant loci |
| `genotype/filter_of_all_codom` | Loci retained after call-rate and minor-allele-frequency filtering |

The two genotype tables contain tag ID, reference-tag sequence, position within the tag, reference base, and sample genotypes in the order recorded by `sample_order.txt`. Missing calls are represented by `--`.
