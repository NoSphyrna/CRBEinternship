#!/bin/bash
#SBATCH -J train_hitac
#SBATCH -o /home/%u/work/job_logs/hitac/output_%j.out
#SBATCH -e /home/%u/work/job_logs/hitac/error_%j.out
#SBATCH -t 48:00:00
#SBATCH --mem=256G
#SBATCH -c 128

#Load modules
module purge

module load devel/Miniconda/Miniconda3
# module load nvidia/CudaToolkit/12.4 # To use gpus

#default
databases="$HOME/work/database"
original="$databases/original"
trimmed="$databases/trimmed"
unite_utax="$trimmed/utax_reference_dataset_all_19.02.2025_full_ITS_ITS5_ITS4_trimmed.fasta"
euk_utax="$trimmed/SINTAX_EUK_ITS_v2.0_full_ITS_ITS5_ITS4_trimmed.fasta"

hitac_models="$HOME/work/HiTaC_models/"
hitac="$HOME/work/conda/envs/hitac"
#Charge config file (a liitle trick to make sure it's form the same directory as the script)

source "$SLURM_SUBMIT_DIR/config_databases.cfg"
source "$SLURM_SUBMIT_DIR/config_nanopore.cfg"

# we activate the conda environment
# TODO :Check if conda env is here and create it if not
source activate "$hitac"
# Checks for directories and files

if [ ! -d "$original" ]; then
	echo "Invalid argument: $original doesn't exist"
	exit 1
fi

if [ ! -d "$trimmed" ]; then
	echo "Invalid argument: $trimmed doesn't exist"
	exit 1
fi

# If model isn't downloaded already, then download it
if [ ! -d "$hitac_models" ]; then
	mkdir -p "$hitac_models"
fi

# We fit the model and register it in the hitac_models folder
# The threads argument is important because auto
# evaluation my think there ar more threads available
hitac-fit \
	--reference "$unite_utax" \
	--threads 128 \
	--classifier "$hitac_models/unite_utax_euk_trimmed_2025_6mers_hitac.pkl"
