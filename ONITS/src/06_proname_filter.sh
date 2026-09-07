#!/bin/bash
#SBATCH -J run_proname_filter
#SBATCH -o /home/%u/work/job_logs/proname/output_%j.out
#SBATCH -e /home/%u/work/job_logs/proname/error_%j.out
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH -c 32

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

module load containers/Apptainer/1.4.1

#default
working_dir="$HOME/work/Nanopore/"
run="$working_dir/run1/"
demux="$run/demux/"
stats="$run/stats/"
proname_dir="$working_dir/proname/run1/"
proname="$proname_dir/../proname_v2.3.0-amd64.sif"
#Charge config file (a litle trick to make sure it's form the same directory as the script)
source "$CONFIG"

# Checkings (particularly import to do this because Cutadapt doesn't handle well missing directories)

if [ ! -d "$demux" ]; then
	echo "Input files not found : $demux doesn't exist" >&2
	exit 1
fi

if [ ! -d "$stats" ]; then
	mkdir -p "$stats"
fi

if [ ! -d "$proname_dir" ]; then
	mkdir -p "$proname_dir"
fi

if [ ! -f "$proname" ]; then
	cd "$proname_dir" || exit 1
	apptainer pull docker://benn888/proname:v2.3.0-amd64
fi

apptainer exec --bind "$proname_dir":/data --pwd /data "$proname" proname_filter \
	--datatype simplex \
	--filtminlen 100 \
	--filtmaxlen 2000 \
	--filtminqual 15 \
	--threads 32 \
	--verbose \
	--inputpath "$demux"

# Count remaing reads
echo -e "Filter\t$(($(wc -l <"$proname_dir/HQ/HQ_simplex_seqs.fastq") / 4))" >>"$stats/nb_reads.tsv"
