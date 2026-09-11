# Data and database availability

Large sequencing datasets, reference genomes, and MAP2B databases are not stored in this source-code repository. Obtain the study datasets and associated identifiers from the data-availability statement of the accompanying article and obtain third-party reference sequences from their original repositories.

For local analyses:

- store large data outside the Git checkout;
- use absolute paths in the module-specific tab-delimited templates under `examples/`;
- use the same Type IIB enzyme setting for library preparation, database construction, and profiling;
- record reference-database versions, source accessions, retrieval dates, and file checksums; and
- do not commit controlled-access or personally identifiable data.

The MAP2B `-s` option must point to the prepared database directory used for the analysis.
