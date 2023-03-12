#!/bin/bash

#SBATCH --job-name=AF90cmgr
#SBATCH --qos=1week
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G            # --mem=1024G will send to bigmem /// 300G for Uniref50
#SBATCH --output=slurm_output/UR50_comp_output%A.out
#SBATCH --error=slurm_output/UR50_comp_error%A.err

python3 make_communities_map.py ../data_generated