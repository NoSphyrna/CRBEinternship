#!/bin/bash
#SBATCH -J run_basecalling
#SBATCH -o /home/%u/work/job_logs/demux/output_%j.out
#SBATCH -e /home/%u/work/job_logs/demux/error_%j.out
#SBATCH -t 24:00:00
#SBATCH --mem=8G
#SBATCH -c 16

# For Cutadapt limiting factor is CPUs not ram

#Load modules
module purge

module load bioinfo/Cutadapt/5.0
module load devel/Miniconda/Miniconda3

#default
merge_fastq="$HOME/work/Nanopore/run1/basecalled_sup/merged.fastq"
stats="$HOME/work/Nanopore/run1/stats/"
pre_demux="$HOME/work/Nanopore/run1/pre_demux/"
demux="$HOME/work/Nanopore/run1/demux/"
linked_adapters="$HOME/work/Nanopore/run1/adapters/linked_adapters.fasta"
biopy="$HOME/work/conda/envs/biopy"
#Charge config file (a liitle trick to make sure it's form the same directory as the script)
source "$SLURM_SUBMIT_DIR/config_nanopore.cfg"

# we activate the conda environment
# TODO :Check if conda env is here and create it if not
source activate "$biopy"
# Checkings (particularly import to do this because Cutadapt doesn't handle well missing directories)
if [ ! -f "$merge_fastq" ]; then
	echo "Invalid argument: $merge_fastq isn't a file"
	exit 1
fi

if [ ! -f "$linked_adapters" ]; then
	echo "Invalid argument: $linked_adapters isn't a file"
	exit 1
fi

if [ ! -d "$stats" ]; then
	mkdir -p "$stats"
fi

if [ ! -d "$pre_demux" ]; then
	mkdir -p "$pre_demux"
fi

if [ ! -d "$demux" ]; then
	mkdir -p "$demux"
fi

# -g file:"$linked_adapt" <- Here we use linked adapters like ADPATFWD...ADAPTREV and the -g option allows to force presence of both adapters for a read to be trimmed
# --revcomp \ # allows to check reverse complement of ther read (in that case write the reversercomplent in th ourput file)
# --rename '{header}\tCT:r:{rc}' \ # This allows to place the rc of the --revcomp after a tab and a prefix to find it with the python script and parse the info-file.tsv
# --cores=0 \ # automatically detects the number of cpu cores available in the job
# -e 0.1 --no-indels \ # Error max at 10% of the read length ~ 3 for the nanopore seq and no indels allowed
# --discard-untrimmed \ # When an adapter has not been found, read is not trimmed then the read is discarded
# -o "demux/{name}.fastq" \ # write the output fastq to the name of the sample (given with the adapter in the linked_adapt.fasta file)
# "$input" # the input fastq file
cutadapt \
	-g file:"$linked_adapters" \
	--revcomp \
	--rename '{header}\tCT:r:{rc}' \
	--cores=0 \
	-e 0.1 --no-indels \
	--discard-untrimmed \
	--info-file "$stats/info.tsv" \
	--json="$stats/demux.cutadapt.json" \
	-o "$pre_demux/{name}.fastq" \
	"$merge_fastq" >"$stats/demux_cutadapt.txt"

# After we need to trimm move tables if they exist :
python "$SLURM_SUBMIT_DIR/../Python_scripts/demux_moves.py" "$stats/info.tsv" "$pre_demux" "$demux"
