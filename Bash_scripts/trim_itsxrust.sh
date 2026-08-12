#!/bin/bash
#SBATCH -J trim_its
#SBATCH -o /home/%u/work/job_logs/databases/output_%j.out
#SBATCH -e /home/%u/work/job_logs/databases/error_%j.out
#SBATCH -t 48:00:00
#SBATCH --mem=64G
#SBATCH -c 32

#Load modules
module purge

module load devel/Miniconda/Miniconda3

source activate ~/work/conda/envs/itsxrust
#default
databases="$HOME/work/database"
original="$databases/original"
itsx="$databases/itsxrust/"

#Charge config file (a liitle trick to make sure it's form the same directory as the script)
source "$SLURM_SUBMIT_DIR/config_databases.cfg"

# Checkings (particularly import to do this because Cutadapt doesn't handle well missing directories)

if [ ! -d "$original" ]; then
	echo "Invalid argument: $original doesn't exist"
	exit 1
fi

if [ ! -d "$itsx" ]; then
	mkdir -p "$itsx"
fi

#### ITS1
file="$original/sh_general_release_dynamic_s_all_19.02.2025_dev.fasta"
echo "File is ok: $file"
NAME_DB="$(basename "$file" ".fasta")_ITS1.fasta"
echo "Output name : $NAME_DB"

# Trimm with itsxrust
itsxrust extract \
	--input "$file" \
	--hmm ~/work/ITSx_HMM/F.hmm \
	--output "$itsx/$NAME_DB" \
	--region its1 \
	--hmmer-cpu 32 \
	--explain 10 \
	--write-skipped "$itsx/tmp_skipped.fasta"

# then add the skipped reads to the final file to keep all reads
cat "$itsx/tmp_skipped.fasta" >>"$itsx/$NAME_DB"
