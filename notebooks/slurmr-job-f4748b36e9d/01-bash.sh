#!/bin/sh
#SBATCH --job-name=slurmr-job-f4748b36e9d
#SBATCH --output=/mnt/scratcha/fmlab/baker02/tidyscreen/notebooks/slurmr-job-f4748b36e9d/02-output-%A-%a.out
#SBATCH --array=1-20
#SBATCH --time=1-00:00:00
#SBATCH --mem=6G
#SBATCH --job-name=slurmr-job-f4748b36e9d
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
export OMP_NUM_THREADS=1
/home/baker02/miniconda3/envs/tidyscreen/lib/R/bin/Rscript  /mnt/scratcha/fmlab/baker02/tidyscreen/notebooks/slurmr-job-f4748b36e9d/00-rscript.r
