#!/bin/bash
#SBATCH -J run_proname
#SBATCH -o /home/%u/work/job_logs/dorado/output_%j.out
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH -c 32

IMPORT=FALSE
FILTER=FALSE
REFINE=FALSE

usage() {
	echo "Usage: $0 -ifr"
	echo "  -i  proname_import is executed"
	echo "  -f  proname_filter is executed"
	echo "  -r  proname_refine is executed"
	echo "  -h  Show this help message"
	exit 0
}

while getopts "hifr" option; do
	case $option in
	h)
		usage
		;;
	i)
		IMPORT=TRUE
		;;
	f)
		FILTER=TRUE
		;;
	r)
		REFINE=TRUE
		;;
	\?)
		echo "Invalid option" >&2
		usage
		;;

	esac
done

if [ $IMPORT = FALSE ] && [ $FILTER = FALSE ] && [ $REFINE = FALSE ]; then
	IMPORT=TRUE
	FILTER=TRUE
	REFINE=TRUE
fi

#Load modules
module purge

module load containers/Apptainer/1.4.1

#default
data="$HOME/work/Nanopore/run1/"
stats="$HOME/work/Nanopore/run1/stats/"
demux="$HOME/work/Nanopore/run1/demux/"
proname_dir="$HOME/work/Nanopore/proname"
proname="$HOME/work/Nanopore/proname/proname_v2.3.0-amd64.sif"
#Charge config file (a litle trick to make sure it's form the same directory as the script)
source "$SLURM_SUBMIT_DIR/config_nanopore.cfg"

# Checkings (particularly import to do this because Cutadapt doesn't handle well missing directories)

if [ ! -d "$stats" ]; then
	mkdir "$stats"
fi

if [ ! -d "$demux" ]; then
	mkdir "$demux"
fi
if [ ! -d "$data" ]; then
	mkdir "$stats"
fi
if [ ! -d "$proname" ]; then
	mkdir "$stats"
fi

if [ ! -f "$proname" ]; then
	cd "$proname_dir" || exit 1
	apptainer pull docker://benn888/proname:v2.3.0-amd64
fi

apptainer run --bind "$data":/data "$proname"

# Check MPLCONFIGDIR

if [ $IMPORT = TRUE ]; then
	proname_import \
		--inputpath "$demux" \
		--threads 32 \
		--duplex no \
		--trimadapters no \
		--trimprimers yes \
		--fwdprimer GTACACACCGCCCGTCG \
		--revprimer CGCCTSCSCTTANTDATATGC \
		--primertrimmingmode hard
fi

if [ $FILTER = TRUE ]; then
	proname_filter \
		--datatype simplex \
		--filtminlen 500 \
		--filtmaxlen 2000 \
		--filtminqual 15 \
		--threads 32 \
		--inputpath "$demux"
fi

if [ $REFINE = TRUE ]; then
	proname_refine \
		--clusterid 0.97 \
		--clusterthreads 32 \
		--inputpath "$demux" \
		--polisher medaka \
		--minreadspercluster 2 \
		--polisherthreads 32 \
		--chimeramethod denovo \
		--polishermodel r1041_e82_400bps_sup_v5.2.0 \
		--qiime2import no
fi
