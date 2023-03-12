#!/bin/bash

#SBATCH --job-name=AF90coms
#SBATCH --qos=1day
#SBATCH --cpus-per-task=128
#SBATCH --mem=50G            # --mem=1024G will send to bigmem 
#SBATCH --output=slurm_output/AF90_communities_output%A.out
#SBATCH --error=slurm_output/AF90_communities_error%A.err

python3 get_communities_summary.py data_generated/AFDB90v4_cc_data.csv data_generated/uniprot_community_taxonomy_map.csv 128