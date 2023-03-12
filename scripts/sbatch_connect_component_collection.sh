#!/bin/bash

#SBATCH --job-name=AF90AvsA
#SBATCH --qos=1day
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G            # --mem=1024G will send to bigmem /// 300G for Uniref50
#SBATCH --output=slurm_output/UR50_comp_output%A.out
#SBATCH --error=slurm_output/UR50_comp_error%A.err

DATABASES="databases"

python3 get_connected_components.py ../data_generated/AFDBv4_90.fasta ../data_generated/AFDB90v4_all-gainst-all.m8 ../data_generated