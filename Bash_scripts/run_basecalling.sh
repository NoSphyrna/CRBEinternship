#!/bin/bash
#SBATCH -J run_basecalling
#SBATCH -o /home/%u/work/job_logs/dorado/output_%j.out
#SBATCH -e /home/%u/work/job_logs/dorado/error_%j.out
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:nvidia_a100:1
#SBATCH -t 24:00:00
#SBATCH --mem=16G
#SBATCH -c 8

#Load modules
module purge

module load bioinfo/Dorado/2.0.1
module load nvidia/CudaToolkit/12.4 # To use gpus

#default
pod5="/home/bperez/save/MITI/Nanopore_run1/"
basecalled="$HOME/work/Nanopore/run1/basecalled_sup/"
model_dir="$HOME/work/Nanopore/dorado_models/"
model_name="dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
model="$model_dir/$model_name"
kit_name="SQK-NBD114-24"

#Charge config file (a liitle trick to make sure it's form the same directory as the script)
source "$SLURM_SUBMIT_DIR/config_nanopore.cfg"

# Checks for directories and files

if [ ! -d "$pod5" ]; then
	echo "Input pod5 folder not found : $pod5 doesn't exist"
	exit 1
fi

if [ ! -d "$basecalled" ]; then
	mkdir -p "$basecalled"
fi

if [ ! -d "$model_dir" ]; then
	mkdir -p "$model_dir"
fi

# If model isn't downloaded already, then download it
if [ ! -f "$model" ]; then
	dorado download --model "$model_name" --directory "$model_dir"
fi

# Basecalling :
# --recursive allows to treat all pod5 files given in the input folder even when it's in different folders
# --device cuda:all allows to make sure we use all available gpus when runnong the basecall
# --emit-fastq chage the output from bam files to fastq files
# --emit-moves for the polishing with dorado
# --kit-name the kit used to sequence the samples
# --emit-summary to have a summary for pycoQC
# if you want to treat separatly pass and fail add a --min-qscore for the treshold between pass and fail
# -o output folder where the arborescence of fastq files whlie be placed
# model : the model used for the basecalling (by default you can choose fast, hac or sup and it will download the corresponding one)
# but here we downloaded first the model with dorado download then used the downloaded model.
# input : the input folder containg all pod5 files we want to treat (as we use recursive the pod5 can be in lower folder
# e.g. input/flow1/pod5/file.pod5 input/flow2/pod5/file2/pod5)
dorado basecaller --recursive --device cuda:all --emit-fastq --emit-moves --kit-name "$kit_name" --emit-summary -o "$basecalled" "$model" "$pod5"
