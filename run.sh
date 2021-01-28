#!/bin/bash
CONDA_ENV=conda
nohup snakemake --cores 2 \
	--use-conda --conda-prefix=$CONDA_ENV \
	--configfile=config.yaml > log.txt & echo $! > pid.txt



