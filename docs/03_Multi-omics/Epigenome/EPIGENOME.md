# Epigenome analysis (MethylRAD-based)

The methylation module at `software/methylation/` supports reference-read construction, fixed-length read extraction, and site-level methylation quantification for CpG and CHG contexts.

1. `build_reference_reads.py` builds methylation-associated reference reads.
2. `extract_fixed_reads.py` extracts matching sample reads.
3. `map_and_quantify_reads.py` maps and generates raw count, RPM, and M-index results.

## Core manuscript route

The default is MspJI (`-e 1`) with CpG (`-t 1`). Because `-e 1` is already default, the basic commands do not repeat it.

```bash
python3 software/methylation/build_reference_reads.py -r /path/reference.fa.gz
python3 software/methylation/extract_fixed_reads.py \
  -i examples/03_Multi-omics/Epigenome/epigenome_samples.tsv
python3 software/methylation/map_and_quantify_reads.py \
  -i examples/03_Multi-omics/Epigenome/epigenome_samples.tsv
```

Run all three from the same project working directory so default paths align:

```text
reference_reads/MspJI_CpG/
clean_reads/MspJI_CpG/
results/MspJI_CpG/
```

## Sample manifest

```text
sample_id	fastq
sample01	/absolute/path/sample01_R1.fastq.gz
```

Use the same two-column manifest for extraction and mapping, with one single-end or R1 read file per sample.

## Main outputs

| Output | Format | Description |
|---|---|---|
| `samples/<sample_id>.sites.tsv.gz` | Gzip-compressed TSV | Per-site coordinates, sequence, raw count, RPM, and M-index |
| `count_matrix.tsv.gz` | Gzip-compressed TSV | Raw-count matrix across samples |
| `rpm_matrix.tsv.gz` | Gzip-compressed TSV | RPM matrix across samples |
| `m_index_matrix.tsv.gz` | Gzip-compressed TSV | M-index matrix across samples |
| `sample_qc.tsv` | TSV | Per-sample processing summary |
| `config.json` | JSON | Run configuration |

The per-sample site table contains `site_id`, `seqname`, `start_1based`, `sequence`, `raw_count`, `rpm`, and `m_index`.

The site-level file is the principal output for downstream analysis.

## Advanced options

Other presets/custom motifs remain available. See `software/methylation/PRESET_TABLE.tsv` and each program's `--help`. For a non-default context/enzyme, use the same `-e` and `-t` across all three stages and do not mix output directories.

```bash
# Advanced MspJI CHG route
python3 software/methylation/build_reference_reads.py -r reference.fa.gz -t 2
python3 software/methylation/extract_fixed_reads.py -i samples.tsv -t 2
python3 software/methylation/map_and_quantify_reads.py -i samples.tsv -t 2
```

The mapping stage requires Bowtie from the supplied environment.
