# Input and example files

This directory provides tab-delimited templates and compact example tables for the computational workflows:

- `manifests/map2b_samples.tsv`: MAP2B sample list;
- `manifests/methylation_samples.tsv`: methylation sample list;
- `manifests/database_genomes.tsv`: database-builder genome manifest;
- `manifests/absolute_quantification_metadata.tsv`: host-concentration metadata;
- `manifests/absolute_quantification_metadata_with_volumes.tsv`: optional library-equivalent metadata;
- `manifests/serial_tag_samples.tsv`: serial-tag library and T1-T5 sample assignments;
- `manifests/metatranscriptome_samples.tsv`: paired metatranscriptome read manifest; and
- `metatranscriptome/`: compact DNA and RNA abundance tables for the integration step.

Copy the relevant template to a working directory, replace the example values and paths, and preserve tab delimiters. Store large input and output files outside the Git repository.
