#!/bin/bash
#SBATCH -J run_taxo
#SBATCH -o /home/%u/work/job_logs/taxo/output_%j.out
#SBATCH -e /home/%u/work/job_logs/taxo/error_%j.out
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH -c 128

# the 128 are for the dnabarcoder threads
# needs hitac conda env

#Load modules
module purge

module load bioinfo/VSEARCH/2.29.3
module load devel/Miniconda/Miniconda3

#default parameters
databases="$HOME/work/database"
original="$databases/original"
trimmed="$databases/trimmed"
unite_utax="$trimmed/utax_reference_dataset_all_19.02.2025_full_ITS_ITS5_ITS4_trimmed.fasta"
euk_utax="$trimmed/SINTAX_EUK_ITS_v2.0_full_ITS_ITS5_ITS4_trimmed.fasta"
kraken="$databases/kraken"
dnabarcoder="$databases/dnabarcoder"

dnabar_ref="$dnabarcoder/unite2025ITS.fasta"
dnabar_class="$dnabarcoder/unite2025ITS.classification"
dnabar_cutoffs="$dnabarcoder/dnabarcoder/unite2025ITS.unique.cutoffs.json"

working_dir="$HOME/work/Nanopore/"
run="$working_dir/run1/"

proname_dir="$working_dir/proname/run1/"
reads_OTU="$proname_dir/rep_seqs.fasta"
taxo="$run/taxo/"
taxo_stats="$taxo/stats/"

hitac_models="$HOME/work/HiTaC_models/"
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

if [ "$METHOD" = "all" ] || [ "$METHOD" = "dnabarcoder" ]; then
	if [ ! -d "$dnabarcoder" ]; then
		echo "Invalid argument : $dnabarcoder doesn't exist"
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
fi

if [ ! -d "$taxo" ]; then
	mkdir -p "$taxo"
fi

if [ ! -d "$taxo_stats" ]; then
	mkdir -p "$taxo_stats"
fi

# The output file for stats of the different classifiers
STATS="$taxo_stats/classifiers_ressource_usage.tsv"

# First we write the header of the stat file : (with -e to translat \t as tab)
echo -e "id\tcommand\ttime(s)\tuser(s)\tsys(s)\tCPU_time\tavgMEM(kB)\tmaxMEM(kB)\tswaps\tdata(kB)\tjob_cores\tjob_mem" >"$STATS"

# Format for time :
TIMEFMT="%C\t%e\t%U\t%S\t%P\t%K\t%M\t%W\t%D"

if [ "$METHOD" = "all" ] || [ "$METHOD" = "vsearch" ]; then
	echo "Running vsearch assignation with vsearch on: $reads_OTU, with unite"
	/usr/bin/time -o "$STATS" -a -f "vsearch_unite\t$TIMEFMT\t$SLURM_CPUS_PER_TASK\t$SLURM_MEM_PER_NODE" vsearch --usearch_global "$reads_OTU" \
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
		--userout "$taxo/taxonomy_OTU_vsearch_unite.tsv"

	# --threads $NB_CORES \ by default it uses all threads available so it's not necessary

	echo "Running vsearch assignation with vsearch on: $reads_OTU, with eukaryome"
	/usr/bin/time -o "$STATS" -a -f "vsearch_euk\t$TIMEFMT\t$SLURM_CPUS_PER_TASK\t$SLURM_MEM_PER_NODE" vsearch --usearch_global "$reads_OTU" \
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
		--userout "$taxo/taxonomy_OTU_vsearch_euk.tsv"

fi
if [ "$METHOD" = "all" ] || [ "$METHOD" = "sintax" ]; then
	echo "Running vsearch assignation with sintax on: $reads_OTU, with unite"
	/usr/bin/time -o "$STATS" -a -f "sintax_unite\t$TIMEFMT\t$SLURM_CPUS_PER_TASK\t$SLURM_MEM_PER_NODE" vsearch --sintax "$reads_OTU" \
		--db "$unite_utax" \
		--sintax_cutoff 0.5 \
		--tabbedout "$taxo/taxonomy_OTU_sintax_unite.tsv" \
		--strand plus

	echo "Running vsearch assignation with sintax on: $reads_OTU, with eukaryome"
	/usr/bin/time -o "$STATS" -a -f "vsearch_euk\t$TIMEFMT" vsearch --sintax "$reads_OTU" \
		--db "$euk_utax" \
		--sintax_cutoff 0.5 \
		--tabbedout "$taxo/taxonomy_OTU_sintax_euk.tsv" \
		--strand plus

fi

if [ "$METHOD" = "all" ] || [ "$METHOD" = "hitac" ]; then

	source activate ~/work/conda/envs/hitac #TODO: check if env exists else create it
	models=("$hitac_models"/*)
	for model in "${models[@]}"; do
		model_name=$(basename "$model")

		echo "Running hitac assignation with model $model_name on: $reads_OTU"
		/usr/bin/time -o "$STATS" -a -f "hitac_$model_name\t$TIMEFMT\t$SLURM_CPUS_PER_TASK\t$SLURM_MEM_PER_NODE" hitac-classify \
			--classifier "$model" \
			--reads "$reads_OTU" \
			--classification "$taxo/taxonomy_OTU_hitac_$model_name.tsv"
	done
	source deactivate
fi

if [ "$METHOD" = "all" ] || [ "$METHOD" = "dnabarcoder" ]; then

	module load bioinfo/NCBI_Blast+/2.15.0+ bioinfo/Krona/2.8.1 bioinfo/IQ-TREE/2.4.0 bioinfo/MAFFT/7.505 bioinfo/ClustalOmega/1.2.4

	module load bioinfo/dnabarcoder/1.0.7

	# We go to the dnbarcoder folder because the out puts will be created there
	cd "$dnabarcoder" || exit 1

	# First serach for best match sequences
	/usr/bin/time -o "$STATS" -a -f "dnabarcoder_search\t$TIMEFMT\t$SLURM_CPUS_PER_TASK\t$SLURM_MEM_PER_NODE" dnabarcoder.py search -i "$reads_OTU" -r unite2025ITS.fasta -o "$taxo"

	# Then we classify using the cutoffs

	/usr/bin/time -o "$STATS" -a -f "dnabarcoder_classify\t$TIMEFMT\t$SLURM_CPUS_PER_TASK\t$SLURM_MEM_PER_NODE" dnabarcoder.py classify -i "$taxo"/"$(basename "$reads_OTU" .fasta)".unite2025ITS_BLAST.bestmatch -c "$dnabar_class" -cutoffs "$dnabar_cutoffs" -o "$taxo"

fi
