#!/bin/bash
#SBATCH -J run_proname_import
#SBATCH -o /home/%u/work/job_logs/proname/output_%j.out
#SBATCH -e /home/%u/work/job_logs/proname/error_%j.out
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH -c 32

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
source "$SLURM_SUBMIT_DIR/config_nanopore.cfg"

# Checkings (particularly import to do this because Cutadapt doesn't handle well missing directories)
#
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

# Check MPLCONFIGDIR

apptainer exec --bind "$proname_dir":/data --pwd /data "$proname" proname_import \
	--inputpath "$demux" \
	--threads 32 \
	--duplex no \
	--trimadapters no \
	--trimprimers yes \
	--fwdprimer GTACACACCGCCCGTCG \
	--revprimer CGCCTSCSCTTANTDATATGC \
	--primertrimmingmode hard \
	--plotformat html
