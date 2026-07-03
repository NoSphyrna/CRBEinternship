#!/bin/bash
#SBATCH -J run_basecalling
#SBATCH -o /home/%u/work/job_logs/dorado/output_%j.out
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:nvidia_a100:1
#SBATCH -t 00:30:00
#SBATCH --mem=128G
#SBATCH -c 32

#Load modules
module purge

module load bioinfo/Dorado/2.0.1
module load nvidia/CudaToolkit/12.4

#default
input="$HOME/work/test_pod5/"
output="$HOME/work/test_output/"
model="$HOME/work/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0/"
kit_name="SQK-NBD114-24"

echo "************ Check GPU ************"
nvidia-smi

echo "******** Basecalling test *********"
time dorado basecaller --recursive --device cuda:all --emit-fastq --kit-name "$kit_name" --emit-summary -o "$output" "$model" "$input"
