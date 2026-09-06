# Serial-tag preprocessing

The serial-tag workflow reconstructs five position-resolved Type IIB restriction tags from paired-end serial-tag libraries. It supports 16 Type IIB restriction enzymes and produces fixed-length FASTA files suitable for MAP2B or MAP2B-Cross-domain analysis.

## Requirements

- Python 3.9 or later;
- PEAR 0.9.6 or a compatible `pear` executable;
- paired-end FASTQ files, optionally gzip-compressed.

PEAR is included in the supplied i2bRAD-M environment. Use the same enzyme identifier throughout merging, candidate-window generation, tag recovery, database construction, and downstream profiling.

## Enzyme identifiers

```bash
python3 software/serial_tag/serial2brad.py --list-enzymes
```

The `-e` option accepts one identifier from `1` to `16` per run. The complete preset definitions are stored in `software/serial_tag/config/enzymes.json` and summarized in `software/serial_tag/config/enzyme_presets.tsv`.

## Article-following directory workflow

Place the paired FASTQ files in `rawdata/`. File names must contain an identifiable R1/R2 marker and a shared library prefix. From the repository root:

```bash
cd software/serial_tag/scripts

python3 extractN.py -d rawdata
python3 goodquality.py -d extractNdata
python3 common_find.py -d goodqualitydata
python3 pear_merge.py -d goodquality-filter -e <enzyme_id>
python3 fq2fasta.py -d 5tagData
python3 split5.py -d 5tagfasta -e <enzyme_id>
python3 recover_sample_tags.py -d split5Data -e <enzyme_id>
```

The default directory chain is:

```text
rawdata
-> extractNdata
-> goodqualitydata
-> goodquality-filter
-> 5tagData
-> 5tagfasta
-> split5Data
-> tagData
```

`split5.py` generates five enzyme-specific overlapping candidate windows. `recover_sample_tags.py` then recovers the final enzyme-defined fixed-length tags. For BsaXI libraries, the recovered tags are 27 bp.

## Manifest controller

For multiple libraries, the complete workflow can also be run with the controller:

```bash
cd software/serial_tag
python3 serial2brad.py \
  -i ../../examples/manifests/serial_tag_samples.tsv \
  -e <enzyme_id> -o ../../results/serial_tag -p 4
```

The manifest must contain the following tab-separated header:

```text
library_id  r1  r2  t1  t2  t3  t4  t5
```

`t1` to `t5` assign the five recovered positions to their downstream sample identifiers.

## Main outputs

- `tagData/`: enzyme-defined fixed-length FASTA files and tag-recovery QC tables from the directory workflow;
- `final_tags/`: sample-level FASTA files from the manifest controller;
- `map2b_manifests/`: two-column sample manifests for downstream MAP2B analysis;
- step-specific and batch-level QC tables generated alongside the intermediate files.

Use unique library and sample identifiers and keep the complete QC output with the analysis record.
