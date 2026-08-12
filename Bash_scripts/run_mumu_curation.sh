#!/bin/bash
#SBATCH -J run_mumu
#SBATCH -o /home/%u/work/job_logs/mumu/output_%j.out
#SBATCH -e /home/%u/work/job_logs/mumu/error_%j.out
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH -c 32

#Load modules
module purge

module load bioinfo/VSEARCH/2.29.3

# For mumu
module load compilers/gcc/15.1.0

#default parameters
working_dir="$HOME/work/Nanopore/"
run="$working_dir/run1/"

proname_dir="$working_dir/proname/run1/"
reads_OTU="$proname_dir/rep_seqs.fasta"
OTU_table="$proname_dir/rep_table.tsv"

mumu_output="$run/final_mumu/"
OTU_TABLE_form="$mumu_output/OTU_table.tsv"
OTU_SEQ="$mumu_output/OTU_seqs.fasta"
OTU_TABLE_MUMU="$mumu_output/OTU_table_mumu.tsv"
MATCH_LIST="$mumu_output/matches.list"
#Charge config file (a liitle trick to make sure it's form the same directory as the script)
source "$SLURM_SUBMIT_DIR/config_nanopore.cfg"

# Checkings (particularly import to do this because Cutadapt doesn't handle well missing directories)

if [ ! -f "$OTU_table" ]; then
	echo "OTU table '$OTU_table' doesn't exist, maybe proname refine wasn't correctly launched"
	exit 1
fi
if [ ! -f "$reads_OTU" ]; then
	echo "OTU seq file '$reads_OTU' doesn't exist, maybe proname refine wasn't correctly launched"
	exit 1
fi

if [ ! -d "$mumu_output" ]; then
	mkdir -p "$mumu_output"
fi

#first we just format the table with OTU as first column

ID=$(head -n1 "$OTU_table" | cut -f1)
sed "s/$ID/OTU/" "$OTU_table" >"$OTU_TABLE_form"

cp "$reads_OTU" "$OTU_SEQ"

# We then create the match list
# previous --id 0.84
vsearch \
	--usearch_global "$OTU_SEQ" \
	--db "$OTU_SEQ" \
	--self \
	--id 0.95 \
	--iddef 1 \
	--userfields query+target+id \
	--maxaccepts 0 \
	--query_cov 0.9 \
	--maxhits 10 \
	--threads 32 \
	--userout - |
	sed -r 's/;size=[0-9]+;//g' >"$MATCH_LIST"

# And then create the curated OTU_table with mumu
mumu \
	--otu_table "$OTU_TABLE_form" \
	--match_list "$MATCH_LIST" \
	--log "$mumu_output/MUMU_merging_stats.txt" \
	--new_otu_table "$OTU_TABLE_MUMU" \
	--minimum_match 95 \
	--minimum_ratio 1 \
	--minimum_ratio_type "min" \
	--minimum_relative_cooccurence 0.90
