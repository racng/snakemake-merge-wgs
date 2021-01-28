# snakemake-merge-wgs
Snakemake pipeline for filtering, lifting over to hg38, and merging WGS gvcf files.

## Requirements 
Python dependencies are managed by snakemake vvia conda environment, but missing linux libraries may cause error during installing of CrossMap (ex. libcurl4-openssl-dev)
Install them with `sudo apt-get install` command before continuing running snakemake pipeline.

## 