#!/bin/bash
##SBATCH --account={{resources.account | default("")}}
#SBATCH --time={{resources.walltime | default("07:00:00")}}
#SBATCH --mem=64G
#SBATCH -n 1
#SBATCH --cpus-per-task=8 # Crucial: Map Snakemake rule's 'threads' resource to cpus-per-task
