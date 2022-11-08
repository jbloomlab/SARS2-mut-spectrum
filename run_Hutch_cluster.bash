#!/bin/bash
#
#SBATCH -c 32

snakemake -j 4 --use-conda --rerun-incomplete
