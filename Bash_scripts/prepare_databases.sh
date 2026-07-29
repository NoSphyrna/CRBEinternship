#!/bin/bash
#SBATCH -J prepare_databases
#SBATCH -o /home/%u/work/job_logs/databases/output_%j.out
#SBATCH -e /home/%u/work/job_logs/databases/error_%j.out
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH -c 16

# This script is adapted from the step00_database_ITS.sh
# For Cutadapt limiting factor is CPUs not ram

#Load modules
module purge

module load bioinfo/Cutadapt/5.0
module load bioinfo/VSEARCH/2.29.3

#default
databases="$HOME/work/database"
original="$databases/original"
trimmed="$databases/trimmed"
kraken="$databases/kraken"
dnabarcoder="$databases/dnabarcoder"

data_script="$HOME/work/Metabarcoding/Illumina/00_prepare_database.sh"
#Charge config file (a liitle trick to make sure it's form the same directory as the script)
source "$SLURM_SUBMIT_DIR/config_databases.cfg"

# Checkings (particularly import to do this because Cutadapt doesn't handle well missing directories)

TRIM=false
KRAKEN=false
DNABAR=false
usage() {
	echo "Usage: $0 [-h][-t][-k][-d]"
	echo "	-t trimm databases placed in the original folder"
	echo "	-k adapt databases placed in the original folder to kraken2"
	echo "	-d adapt databases placed in the original folder to dnabarcoding"
	exit 1
}

while getopts "tkdh" opt; do
	case $opt in
	t)
		TRIM=true
		;;
	k)
		KRAKEN=true
		;;
	d)
		DNABAR=true
		;;
	h)
		usage
		;;
	\?)
		echo "Invalid option" >&2
		usage
		;;
	esac
done

if [ "$TRIM" = false ] && [ "$KRAKEN" = false ] && [ "$DNABAR" = false ]; then
	TRIM=true
	KRAKEN=true
	DNABAR=true
fi

if [ ! -d "$original" ]; then
	echo "Invalid argument: $original doesn't exist"
	exit 1
fi

if [ ! -d "$trimmed" ]; then
	mkdir -p "$trimmed"
fi

if [ ! -d "$kraken" ]; then
	mkdir -p "$kraken"
fi

if [ ! -d "$dnabarcoder" ]; then
	mkdir -p "$dnabarcoder"
fi

if [ "$TRIM" = true ]; then
	#### ITS1
	FILES=("$original"/*)
	echo "Files were ok: ${FILES[@]}"
	for file in "${FILES[@]}"; do
		NAME_DB="$(basename "$file" ".fasta")_ITS1_fung02"
		echo "First name : $NAME_DB"

		bash "$data_script" \
			-i "$file" \
			-f GGAAGTAAAAGTCGTAACAAGG \
			-r CAAGAGATCCGTTGYTGAAAGTK \
			-o "$trimmed" \
			-n "$NAME_DB" \
			-x false \
			-y false \
			-l 80

		#### full ITS

		NAME_DB="$(basename "$file" ".fasta")_full_ITS_ITS5_ITS4"

		echo "Second name : $NAME_DB"
		# forward: ITS5 (Fung02_F)
		# reverse: ITS4ngsUni

		bash "$data_script" \
			-i "$file" \
			-f GGAAGTAAAAGTCGTAACAAGG \
			-r CGCCTSCSCTTANTDATATGC \
			-o "$trimmed" \
			-n "$NAME_DB" \
			-x false \
			-y false \
			-l 300
	done
fi

# if [ "$KRAKEN" = true ]; then
#
# 	module load bioinfo/Kraken2/2.17.1
#
# 	FILES=("$trimmed"*)
# 	for file in "${FILES[@]}"; do
#
# 	done
# fi
