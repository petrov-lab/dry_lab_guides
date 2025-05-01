#! /bin/bash
#
#SBATCH --job-name=snakemake
#SBATCH --time=7-00:00
#SBATCH --partition=dpetrov,hns
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

. /home/users/jahemker/.bashrc
conda activate snakemake

snakemake --profile profile/
