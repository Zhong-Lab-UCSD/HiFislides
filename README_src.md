# Workflow

Scripts can be found at [./bin/src](./bin/src).

- `hifi_barcode_processing.sh`: deduplicate the barcode reads and generate bwa index of their sequence.

```
-b : Directory of the scripts.
-F : Flowcell ID.
-f : Flowcell type (NextSeq, MiniSeq).
-d : Directory of the barcode fastq files.
-N : Suffix of the fastq file (example R1_001.fastq.gz).
-t : Max CPU threads for parallelized processing, at least 4 (default 8).
-o : Parent output directory.
```

- `hifi_wrapper.sh`: main processing pipeline.
- `hifi_extract_roi.sh`: subset final output data to extract only tiles under ROI.


