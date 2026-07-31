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

if [ "$DNABAR" = true ]; then

	# Dependencies
	module load bioinfo/NCBI_Blast+/2.15.0+ bioinfo/Krona/2.8.1 bioinfo/IQ-TREE/2.4.0 bioinfo/MAFFT/7.505 bioinfo/ClustalOmega/1.2.4

	module load bioinfo/dnabarcoder/1.0.7

	dnabarcoder_dir=$(dirname "$(which dnabarcoder.py)")

	# First transform the database according to
	"$dnabarcoder_dir"/aidscripts/filterClassificationFromSequenceHeaders.py -i "$original/sh_general_release_dynamic_s_all_19.02.2025_dev.fasta" -prefix "$dnabarcoder/unite2025ITS"

	#Compute cutoffs estimations (adapatation of the script unite2024.sh)
	cd "$dnabarcoder" || exit 1

	#select unique sequences
	"$dnabarcoder_dir"/aidscripts/selectsequences.py -i unite2025ITS.fasta -unique yes -c unite2025ITS.classification -o unite2025ITS.unique.fasta
	dnabarcoder.py length -i unite2025ITS.unique.fasta -l 100
	dnabarcoder.py overview -i unite2025ITS.unique.fasta -c unite2025ITS.unique.classification

	#select sequences having taxonomic information at the species level
	"$dnabarcoder_dir"/aidscripts/selectsequences.py -i unite2025ITS.unique.fasta -c unite2025ITS.unique.classification -rank species -o unite2025ITS.unique.species.fasta
	#predict similarity cutoffs for all the genera. For this big dataset, we do not compute the species similarity cutoffs for higher taxonomic levels
	dnabarcoder.py predict -i unite2025ITS.unique.species.fasta -c unite2025ITS.unique.species.classification -st 0.7 -et 1 -s 0.001 -rank species -prefix unite2025ITS.unique -higherrank genus -maxproportion 0.9 -removecomplexes yes
	#predict a global similarity cutoff
	dnabarcoder.py predict -i unite2025ITS.unique.species.fasta -c unite2025ITS.unique.species.classification -st 0.9 -et 1 -s 0.001 -rank species -prefix unite2025ITS.unique -removecomplexes yes
	#remove all the created files
	rm unite2025ITS.unique.species.*

	#genus
	"$dnabarcoder_dir"/aidscripts/selectsequences.py -i unite2025ITS.unique.fasta -c unite2025ITS.unique.classification -rank genus -maxseqnopergroup 1000 -o unite2025ITS.unique.genus.fasta
	dnabarcoder.py predict -i unite2025ITS.unique.genus.fasta -c unite2025ITS.unique.genus.classification -st 0.7 -et 1 -s 0.001 -rank genus -higherrank family -prefix unite2025ITS.unique -maxproportion 0.75
	dnabarcoder.py predict -i unite2025ITS.unique.genus.fasta -c unite2025ITS.unique.genus.classification -st 0.7 -et 1 -s 0.001 -rank genus -prefix unite2025ITS.unique
	rm unite2025ITS.unique.genus.*

	#family
	"$dnabarcoder_dir"/aidscripts/selectsequences.py -i unite2025ITS.unique.fasta -c unite2025ITS.unique.classification -rank family -maxseqnopergroup 1000 -o unite2025ITS.unique.family.fasta
	dnabarcoder.py predict -i unite2025ITS.unique.family.fasta -c unite2025ITS.unique.family.classification -st 0.5 -et 1 -s 0.001 -rank family -higherrank order -prefix unite2025ITS.unique -maxproportion 0.75
	dnabarcoder.py predict -i unite2025ITS.unique.family.fasta -c unite2025ITS.unique.family.classification -st 0.5 -et 1 -s 0.001 -rank family -prefix unite2025ITS.unique
	rm unite2025ITS.unique.family.*

	#order
	"$dnabarcoder_dir"/aidscripts/selectsequences.py -i unite2025ITS.unique.fasta -c unite2025ITS.unique.classification -rank order -maxseqnopergroup 1000 -o unite2025ITS.unique.order.fasta
	dnabarcoder.py predict -i unite2025ITS.unique.order.fasta -c unite2025ITS.unique.order.classification -st 0.5 -et 1 -s 0.001 -rank order -higherrank class -prefix unite2025ITS.unique -maxproportion 0.75
	dnabarcoder.py predict -i unite2025ITS.unique.order.fasta -c unite2025ITS.unique.order.classification -st 0.5 -et 1 -s 0.001 -rank order -prefix unite2025ITS.unique
	rm unite2025ITS.unique.order.*

	#class
	"$dnabarcoder_dir"/aidscripts/selectsequences.py -i unite2025ITS.unique.fasta -c unite2025ITS.unique.classification -rank class -maxseqnopergroup 1000 -o unite2025ITS.unique.class.fasta
	dnabarcoder.py predict -i unite2025ITS.unique.class.fasta -c unite2025ITS.unique.class.classification -st 0.5 -et 1 -s 0.001 -rank class -higherrank phylum -prefix unite2025ITS.unique -maxproportion 0.75
	dnabarcoder.py predict -i unite2025ITS.unique.class.fasta -c unite2025ITS.unique.class.classification -st 0.5 -et 1 -s 0.001 -rank class -prefix unite2025ITS.unique
	rm unite2025ITS.unique.class.*

	#visualization
	dnabarcoder.py predict -i unite2025ITS1.unique.fasta -c unite2025ITS1.unique.classification -rank species,genus,family,order,class

	#compute best cutoffs.
	#If at a taxonomic level, the confidence measure obtained for of a clade is lower
	#than the confidence measure obtained for all sequences then the similarity cutoff predicted for all will be taken.
	#This is to avoid to the problem that sequences being classified wrongly due to the fact that some groups are in need for reclassification.

	dnabarcoder.py best -i dnabarcoder/unite2025ITS.unique.cutoffs.json -c unite2025ITS.unique.classification -mincutoff 0.71

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
