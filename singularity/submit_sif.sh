#!/bin/bash
#
#SBATCH --job-name=sif
#SBATCH --partition=dpetrov,hns,normal
#SBATCH --mail-type=BEGIN,END,FAIL
#
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

outdir="$SCRATCH/singularity_images"
mkdir -p $outdir
sif="$outdir/bcftools.sif"
dockerhub_link="docker://staphb/bcftools:1.22"

singularity build $sif $dockerhub_link
