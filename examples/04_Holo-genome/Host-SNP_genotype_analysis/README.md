# Host-SNP genotype analysis example layout

The Host-SNP workflow accepts sample and reference files directly rather than through a tab-separated manifest. Arrange the input and generated files using the working-directory structure below:

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

See [`HOST_SNP.md`](../../../docs/04_Holo-genome/Host-SNP_genotype_analysis/HOST_SNP.md) for the commands and file-placement steps.
