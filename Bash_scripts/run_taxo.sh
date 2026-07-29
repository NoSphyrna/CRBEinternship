#!/bin/bash
#SBATCH -J run_taxo
#SBATCH -o /home/%u/work/job_logs/taxo/output_%j.out
#SBATCH -e /home/%u/work/job_logs/taxo/error_%j.out
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH -c 32

#Load modules
module purge

module load bioinfo/VSEARCH/2.29.3

#default
databases="$HOME/work/database"
original="$databases/original"
trimmed="$databases/trimmed"
unite_utax="$trimmed/utax_reference_dataset_all_19.02.2025_full_ITS_ITS5_ITS4_trimmed.fasta"
euk_utax="$trimmed/SINTAX_EUK_ITS_v2.0_full_ITS_ITS5_ITS4_trimmed.fasta"
kraken="$databases/kraken"
dnabarcoder="$databases/dnabarcoder"

working_dir="$HOME/work/Nanopore/"
run="$working_dir/run1/"

proname_dir="$working_dir/proname/run1/"
reads_OTU="$proname_dir/rep_seqs.fasta"
taxo="$run/taxo/"
#Charge config file (a liitle trick to make sure it's form the same directory as the script)
source "$SLURM_SUBMIT_DIR/config_databases.cfg"
source "$SLURM_SUBMIT_DIR/config_nanopore.cfg"

# Checkings (particularly import to do this because Cutadapt doesn't handle well missing directories)

METHOD=all
usage() {
	echo "Usage: $0 [-h][-m method]"
	echo "	-m method used for the taxonomic assignation (default all) available :
			vsearch, sintax, dnabarcoder"
	exit 1
}

while getopts "m:h" opt; do
	case $opt in
	m)
		METHOD="$OPTARG"
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

if [ ! -d "$original" ]; then
	echo "Invalid argument: $original doesn't exist"
	exit 1
fi

if [ ! -d "$taxo" ]; then
	mkdir -p "$taxo"
fi

if [ "$METHOD" = "all" ] || [ "$METHOD" = "vsearch" ]; then
	echo "Running vsearch assignation with vsearch on: $reads_OTU, with unite"
	vsearch --usearch_global "$reads_OTU" \
		--dbmask none \
		--qmask none \
		--rowlen 0 \
		--notrunclabels \
		--userfields query+id+target \
		--maxaccepts 10 \
		--maxrejects 32 \
		--db "$unite_utax" \
		--id 0.7 \
		--iddef 2 \
		--userout "$taxo/taxonomy_OTU_vsearch_unite.txt"

	# --threads $NB_CORES \ by default it uses all threads available so it's not necessary

	echo "Running vsearch assignation with vsearch on: $reads_OTU, with eukaryome"
	vsearch --usearch_global "$reads_OTU" \
		--dbmask none \
		--qmask none \
		--rowlen 0 \
		--notrunclabels \
		--userfields query+id+target \
		--maxaccepts 10 \
		--maxrejects 32 \
		--db "$euk_utax" \
		--id 0.7 \
		--iddef 2 \
		--userout "$taxo/taxonomy_OTU_vsearch_euk.txt"

fi
if [ "$METHOD" = "all" ] || [ "$METHOD" = "sintax" ]; then
	echo "Running vsearch assignation with sintax on: $reads_OTU, with unite"
	vsearch --sintax "$reads_OTU" \
		--db "$unite_utax" \
		--sintax_cutoff 0.5 \
		--tabbedout "$taxo/taxonomy_OTU_sintax_unite.txt" \
		--strand plus

	echo "Running vsearch assignation with sintax on: $reads_OTU, with eukaryome"
	vsearch --sintax "$reads_OTU" \
		--db "$euk_utax" \
		--sintax_cutoff 0.5 \
		--tabbedout "$taxo/taxonomy_OTU_sintax_euk.txt" \
		--strand plus

fi
