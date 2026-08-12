#!/bin/bash
#SBATCH -J verify_dnabar
#SBATCH -o /home/%u/work/job_logs/verify/output_%j.out
#SBATCH -e /home/%u/work/job_logs/verify/error_%j.out
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH -c 128

# the 128 are for the dnabarcoder threads
# needs hitac conda env

#Load modules
module purge

module load devel/Miniconda/Miniconda3

#default parameters
databases="$HOME/work/database"
dnabarcoder="$databases/dnabarcoder"
working_dir="$HOME/work/Nanopore/"
run="$working_dir/run1/"

proname_dir="$working_dir/proname/run1/"
reads_OTU="$proname_dir/rep_seqs.fasta"
taxo="$run/taxo/"

dnabar_ref="$dnabarcoder/unite2025ITS.fasta"
dnabar_class="$dnabarcoder/unite2025ITS.classification"
dnabar_cutoffs="$dnabarcoder/dnabarcoder/unite2025ITS.unique.cutoffs.json"
dnabar_classified="$taxo/$(basename "$reads_OTU" .fasta).unite2025ITS_BLAST.classified"

#Charge config file (a liitle trick to make sure it's form the same directory as the script)
source "$SLURM_SUBMIT_DIR/config_databases.cfg"
source "$SLURM_SUBMIT_DIR/config_nanopore.cfg"

# Checkings (particularly import to do this because Cutadapt doesn't handle well missing directories)

if [ ! -d "$dnabarcoder" ]; then
	echo "Invalid argument: $dnabarcoder doesn't exist"
	exit 1
fi

if [ ! -f "$dnabar_classified" ]; then
	echo "Classified file '$dnabar_classified' doesn't exist"
	exit 1
fi

if [ ! -f "$dnabar_ref" ]; then
	echo "Reference fasta file '$dnabar_ref' doesn't exist"
	exit 1
fi

if [ ! -f "$dnabar_class" ]; then
	echo "Classification file '$dnabar_class' doesn't exist"
	exit 1
fi

if [ ! -f "$dnabar_cutoffs" ]; then
	echo "Cutoffs file '$dnabar_cutoffs' doesn't exist"
	exit 1
fi

module load bioinfo/NCBI_Blast+/2.15.0+ bioinfo/Krona/2.8.1 bioinfo/IQ-TREE/2.4.0 bioinfo/MAFFT/7.505 bioinfo/ClustalOmega/1.2.4

module load bioinfo/dnabarcoder/1.0.7

# We go to the dnbarcoder folder because the out puts will be created there
cd "$dnabarcoder" || exit 1

dnabarcoder.py verify -i "$dnabar_classified" -c "$dnabar_class" -r unite2025ITS.fasta -f "$reads_OTU" -rank -cutoffs "$dnabar_cutoffs" -method cutoff -o "$taxo"
