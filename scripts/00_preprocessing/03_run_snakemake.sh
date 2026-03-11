#!/bin/bash
#SBATCH --job-name=flair                                                             
#SBATCH --time=4:00:00                               
#SBATCH --mem=1G
#SBATCH -n 1 # threaded 
#SBATCH -o slurm.flair.job.%j.out
#SBATCH -e slurm.flair.job.%j.err

# activate conda environment
source $HOME/.bash_profile
#module load python3
conda activate flair

# change directory to where Snakefile is located
CWD="/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/00_preprocessing"
cd $CWD

snakemake --nolock -s Snakefile.Flair --jobs 6 --executor slurm --profile slurm_profile --rerun-incomplete --default-resources mem_mb=64000 ntasks=1 threads=8 runtime=120 cpus_per_task=8
# When submitting Snakemake:
