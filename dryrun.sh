#!/bin/bash
CONDA_ENV=conda
snakemake --cores 8 \
	--use-conda --conda-prefix=$CONDA_ENV \
	--configfile=config.yaml -n -r




