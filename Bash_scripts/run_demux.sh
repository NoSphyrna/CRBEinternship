#!/bin/bash
#SBATCH -J run_basecalling
#SBATCH -o /home/%u/work/job_logs/dorado/output_%j.out
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH -c 16

#Load modules
module purge

module load bioinfo/Cutadapt/5.0

#default
merge_fastq="$HOME/work/Nanopore/run1/basecalled_sup/merged.fastq"
stats="$HOME/work/Nanopore/run1/stats/"
demux="$HOME/work/Nanopore/run1/demux/"
linked_adapters="$HOME/work/Nanopore/Adapters/linked_adapters.fasta"
#Charge config file (a liitle trick to make sure it's form the same directory as the script)
source "$SLURM_SUBMIT_DIR/config_nanopore.cfg"

# Checkings (particularly import to do this because Cutadapt doesn't handle well missing directories)
if [ ! -f "$merge_fastq" ]; then
	echo "Invalid argument: $merge_fastq isn't a file"
fi

if [ ! -f "$linked_adapters" ]; then
	echo "Invalid argument: $linked_adapters isn't a file"
fi

if [ ! -d "$stats" ]; then
	mkdir "$stats"
fi

if [ ! -d "$demux" ]; then
	mkdir "$demux"
fi

# -g file:"$linked_adapt" <- Here we use linked adapters like ADPATFWD...ADAPTREV and the -g option allows to force prensence of both adapters for a read to be trimmed
# --revcomp \ # allows to check reverse complement of ther read (in that case write the reversercomplent in th ourput file)
# --cores=0 \ # automatically detects the number of cpu cores available in the job
# -e 0.1 --no-indels \ # Error max at 10% of the read length ~ 3 for the nanopore seq and no indels allowed
# --discard-untrimmed \ # When an adapter has not been found, read is not trimmed then the read is discarded
# -o "demux/{name}.fastq" \ # write the output fastq to the name of the sample (given with the adapter in the linked_adapt.fasta file)
# "$input" # the input fastq file
cutadapt \
	-g file:"$linked_adapters" \
	--revcomp \
	--cores=0 \
	-e 0.1 --no-indels \
	--discard-untrimmed \
	--json="$stats/demux.cutadapt.json" \
	-o "$demux/{name}.fastq" \
	"$merge_fastq" >"$stats/demux_cutadapt.txt"
