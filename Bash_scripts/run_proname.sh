#!/bin/bash
#SBATCH -J run_proname
#SBATCH -o /home/%u/work/job_logs/proname/output_%j.out
#SBATCH -e /home/%u/work/job_logs/proname/error_%j.out
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
working_dir="$HOME/work/Nanopore/"
run="$working_dir/run1/"
demux="$run/demux/"
stats="$run/stats/"
proname_dir="$working_dir/proname/run1/"
proname="$proname_dir/../proname_v2.3.0-amd64.sif"
#Charge config file (a litle trick to make sure it's form the same directory as the script)
source "$SLURM_SUBMIT_DIR/config_nanopore.cfg"

# Checkings (particularly import to do this because Cutadapt doesn't handle well missing directories)

if [ ! -d "$stats" ]; then
	mkdir -p "$stats"
fi

if [ ! -d "$demux" ]; then
	mkdir -p "$demux"
fi
if [ ! -d "$proname_dir" ]; then
	mkdir -p "$proname_dir"
fi

if [ ! -f "$proname" ]; then
	cd "$proname_dir" || exit 1
	apptainer pull docker://benn888/proname:v2.3.0-amd64
fi

apptainer run --bind "$proname_dir":/data "$proname"
echo "--- Test binding data file ---"
echo "Current directory after running docker image:"
pwd
cd "$proname_dir" || exit 1

echo "Directoty after cd:"
pwd
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
		--primertrimmingmode hard \
		--plotformat html
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
		--clusteringmethod vsearch \
		--inputpath "$demux" \
		--polisher dorado \
		--minreadspercluster 2 \
		--polisherthreads 32 \
		--chimeramethod denovo \
		--polishermodel r1041_e82_400bps_sup_v5.2.0 \
		--qiime2import no
fi
