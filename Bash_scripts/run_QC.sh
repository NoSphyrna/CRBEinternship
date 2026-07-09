#!/bin/bash
#SBATCH -J run_basecalling
#SBATCH -o /home/%u/work/job_logs/dorado/output_%j.out
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH -c 16

#Load modules
module purge

module load bioinfo/FastQC/0.12.1
module load bioinfo/pycoQC/2.5.2

#default
basecalled="$HOME/work/Nanopore/run1/basecalled_sup/"
stats="$HOME/work/Nanopore/run1/stats/"

#Charge config file (a liitle trick to make sure it's form the same directory as the script)
source "$SLURM_SUBMIT_DIR/config_nanopore.cfg"

# Merge all fastq files in one fastq
find "$basecalled" -name "*.fastq" -not -name "merged.fastq" -exec cat {} + >"$basecalled/merged.fastq"

# PycoQC
pycoQC -f "$basecalled/sequencing_summary.txt" -o "$stats/pycoQC.html"

# FastQC
fastqc -o "$stats" -memory 64G -t 16 "$basecalled/merged.fastq"
