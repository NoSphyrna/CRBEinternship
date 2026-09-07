# ==================================================================== #
# ==================== ONITS Default Config file ===================== #
# ==================================================================== #

# This is the default config file for the ONITS pipeline if you need to
# modify it, a good practice is to create a copy and modify it.
# Tips : Create multiple config files for each of your tests/runs

# ============================= Options ============================== #
# This part contains the options of the options for each part of the
# pipeline that you might most likely want to change

# Working directory
WORKING_DIR="$HOME/work/ONITS/"
RUN_NAME="run1"

#Inputs
POD5="/home/bperez/save/MITI/Nanopore_run1/"

#Databases
databases="$HOME/work/database"
trimmed="$databases/trimmed"
unite_utax="$trimmed/utax_reference_dataset_all_19.02.2025_full_ITS_ITS5_ITS4_trimmed.fasta"
euk_utax="$trimmed/SINTAX_EUK_ITS_v2.0_full_ITS_ITS5_ITS4_trimmed.fasta"

# Basecalling
kit_name="SQK-NBD114-24"
model_name="dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
# demux ?
DEMUX_b="True"
# ITS extraction ?
ITSXRUST_b="True"
# Primers
primer_fwd="GTACACACCGCCCGTCG"
primer_rev="CGCCTSCSCTTANTDATATGC"
# Filters
min_quality=15
min_length=300
max_length=2000
#Mumu
id_match=0.95
# Taxonomic assignations
sintax_prob_cutoff=0.5
vsearch_pct_cutoff=0.97

# =========================== Installations ========================== #
# This part contains places for the installations of dependeces

CONDA_ENVS="$HOME/work/conda/envs/"
BIOPY="$CONDA_ENVS/biopy"
ITSXRUST="$CONDA_ENVS/itsxrust"
DORADO_MODELS="$WORKING_DIR/dorado_models/"
PRONAME_SIF="$WORKING_DIR/proname/"
APPS="$HOME/work/app/" # For mumu

# ========================== Generic output folders ================== #

# =========================== Older config file ====================== #
working_dir="$HOME/work/Nanopore/"
run="$working_dir/run1/"

# Basecalling
pod5="/home/bperez/save/MITI/Nanopore_run1/"
basecalled="$run/basecalled_sup/"
model_dir="$working_dir/dorado_models/"
model="$model_dir/$model_name"

# Merging
merge_fastq="$basecalled/merged.fastq"

# Demultiplexing
linked_adapters="$run/adapters/linked_adapters.fasta"
pre_demux="$run/pre_demux/"
demux="$run/demux/"
stats="$run/stats/"

# Proname
proname_dir="$working_dir/proname/run1/"
proname="$proname_dir/../proname_v2.3.0-amd64.sif"

# MUMU
reads_OTU="$proname_dir/rep_seqs.fasta"
OTU_table="$proname_dir/rep_table.tsv"

mumu_output="$run/final_mumu/"
OTU_TABLE_form="$mumu_output/OTU_table.tsv"
OTU_SEQ="$mumu_output/OTU_seqs.fasta"
OTU_TABLE_MUMU="$mumu_output/OTU_table_mumu.tsv"
MATCH_LIST="$mumu_output/matches.list"

# Taxo
taxo="$run/taxo/"
taxo_stats="$taxo/stats/"

dnabar_ref="$dnabarcoder/unite2025ITS.fasta"
dnabar_class="$dnabarcoder/unite2025ITS.classification"
dnabar_cutoffs="$dnabarcoder/dnabarcoder/unite2025ITS.unique.cutoffs.json"
dnabar_classified="$taxo/$(basename "$reads_OTU" .fasta).unite2025ITS_BLAST.classified"

hitac_models="$HOME/work/HiTaC_models/"
