#!/bin/bash
#SBATCH -J run_verify
#SBATCH -o $HOME/work/output_verify.out
#SBATCH --mem=8G
#SBATCH -c 4

module purge

module load compilers/gcc/15.1.0
module load statistics/R/4.5.0

#default
input_clean="$HOME/save/cleaned_OTU_tables"
outut_pq_verify="$HOME/work/data_taxinfo/verified_pq_default"

source ./config.cfg

Rscript ../R_scripts/verify_gna_from_clean_data.R "$input_clean" "$outut_pq_verify"
