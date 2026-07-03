#!/bin/bash
#SBATCH -J run_basecalling
#SBATCH -o /home/%u/work/job_logs/dorado/output_%j.out
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:nvidia_a100:1
#SBATCH -t 24:00:00
#SBATCH --mem=16G
#SBATCH -c 8

#Load modules
module purge

module load bioinfo/Dorado/2.0.1
module load nvidia/CudaToolkit/12.4

#default
input="/home/bperez/save/MITI/Nanopore_run1"
output="$HOME/work/output_dorado_sup/"
model="$HOME/work/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0/"
kit_name="SQK-NBD114-24"

dorado basecaller --recursive --device cuda:all --emit-fastq --kit-name "$kit_name" --emit-summary -o "$output" "$model" "$input"
