#!/bin/bash
#SBATCH -J run_taxinfo
#SBATCH -o /home/%u/work/job_logs/output_%j.out
#SBATCH --mem=8G
#SBATCH -c 4

#Load modules
module purge

module load compilers/gcc/15.1.0
module load statistics/R/4.5.0

#default
input_clean="$HOME/save/cleaned_OTU_tables"
pq_verify="$HOME/work/data_taxinfo/verified_pq"
pq_occur="$HOME/work/data_taxinfo/occur_pq"
pq_range="$HOME/work/data_taxinfo/range_pq"
pq_traits="$HOME/work/data_taxinfo/traits_pq"
plot_traits="$HOME/work/plots/traits"
plot_maps="$HOME/work/plots/maps"
plot_range="$HOME/work/plots/range"

#Charge config file (a liitle trick to make sure it's form the same directory as the script)
R_SCRIPTS="$SLURM_SUBMIT_DIR/../R_scripts"
source "$SLURM_SUBMIT_DIR/config.cfg"

if [ "$#" -lt 1 ]; then
	Rscript "$R_SCRIPTS"/verify_gna_from_clean_data.R "$input_clean" "$pq_verify"
else
	case "$1" in
	"verify")
		Rscript "$R_SCRIPTS"/verify_gna_from_clean_data.R "$input_clean" "$pq_verify"
		;;
	"occur")
		Rscript "$R_SCRIPTS"/add_occur.R "$pq_verify" "$pq_occur"
		;;
	*)
		echo "Unknown arg : $1" >&2
		exit 1
		;;
	esac
fi
