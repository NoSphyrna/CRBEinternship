#!/bin/bash
#SBATCH -J run_basecalling
#SBATCH -o /home/%u/work/job_logs/dorado/output_%j.out
#SBATCH -e /home/%u/work/job_logs/dorado/error_%j.out
#SBATCH -t 24:00:00
#SBATCH --mem=16G
#SBATCH -c 8

# Get the config file as input
usage() {
	echo "Usage: $0 [-h] <config_file>"
	exit 1
}

while getopts "h" opt; do
	case $opt in
	h)
		usage
		;;
	\?)
		echo "Invalid option" >&2
		usage
		;;
	esac
done

if [ $# != 1 ]; then
	echo "Invalid number of arguments : $#"
	usage
fi

if [ ! -f "$1" ] || ! CONFIG="$1"; then
	echo "Invalid config file : $1"
	usage
fi

#Load modules
module purge

module load bioinfo/FastQC/0.12.1
module load bioinfo/pycoQC/2.5.2

#default
basecalled="$HOME/work/Nanopore/run1/basecalled_sup/"
stats="$HOME/work/Nanopore/run1/stats/"

merge_fastq="$basecalled/merged.fastq"
#Charge config file (a liitle trick to make sure it's form the same directory as the script)
source "$CONFIG"

if [ ! -d "$basecalled" ]; then
	echo "Couldn't find the basecall folder : $basecalled doesn't exist"
	exit 1
fi
if [ ! -d "$stats" ]; then
	mkdir -p "$stats"
fi

name=$(basename "$merge_fastq")
# Merge all fastq files in one fastq
find "$basecalled" -name "*.fastq" -not -name "$name" -exec cat {} + >"$merge_fastq"
# Count reads
echo -e "Basecalling\t$(($(wc -l <$merge_fastq) / 4))" >"$stats/nb_reads.tsv"

# PycoQC
pycoQC -f "$basecalled/sequencing_summary.txt" -o "$stats/pycoQC.html"

# FastQC
fastqc -o "$stats" -memory 16G -t 8 "$basecalled/merged.fastq"
