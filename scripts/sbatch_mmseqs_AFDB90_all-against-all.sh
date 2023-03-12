#!/bin/bash

#SBATCH --job-name=UR50AvsA
#SBATCH --qos=1day
#SBATCH --cpus-per-task=60
#SBATCH --mem=100G            # --mem=1024G will send to bigmem
#SBATCH --output=slurm_output/UR50_mmseqs_output%A.out
#SBATCH --error=slurm_output/UR50_mmseqs_error%A.err

DATABASES="../databases"
DBFOLDER="${DATABASES}/mmseqsDBs"

# load MMseqs
ml MMseqs2

# # # create database for mmseqs2
# # mkdir $DBFOLDER
mmseqs createdb AFDBv4_90.fasta ${DBFOLDER}/AFDB90v4

# run mmseqs
mmseqs easy-search AFDBv4_90.fasta ${DBFOLDER}/AFDB90v4 ../data_generated/AFDB90v4_all-gainst-all.m8 tmp -e 1e-4 --threads 60
