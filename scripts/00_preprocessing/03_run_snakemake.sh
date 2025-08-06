#!/bin/bash
#SBATCH --job-name=isoseq                                                              
#SBATCH --time=10:00:00                               
#SBATCH --mem=2G
#SBATCH -n 1 # threaded 
#SBATCH --cpus-per-task=1
#SBATCH -o slurm.isoseq.job.%j.out
#SBATCH -e slurm.isoseq.job.%j.err

# activate conda environment
source $HOME/.bash_profile
module load python3
conda activate isoseq

# change directory to where Snakefile is located
CWD="/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/00_preprocessing"
cd $CWD
snakemake --nolock -s Snakefile --jobs 6 --executor slurm --profile slurm_profile --rerun-incomplete --default-resources mem_mb=120000 ntasks=1 threads=8 runtime=550 cpus_per_task=8