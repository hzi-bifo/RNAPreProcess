# RNA Seq Pre-Processing Pipeline

Preforms quality metrics and alignment of raw rna seq fastq files.

To Run:

RNAPreProcess is currently configured for use on BIFO Cluster, however paths can be updates via the config file. Prior to running the pipeline the configuration file will need to be updated.

geneLength	A file of geneLengths for FPKM normalization, formatted as output of gtftools with columns names: gene, mean, median, longest_isoform, merged

```bash
bash RNAprocess.sh -c <path to config file>
```

For additional arguments and use options use: 

```bash
bash RNAprprocess.sh -h 
```

Arguments:



Software Requirements:

STAR - 2.78a
SRA Tool Kit - 2.10.9
edgeR

# Methods

