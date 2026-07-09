# Basecall for ONT

## Module load

```bash
module load bioinfo/Dorado/2.0.1
module load nvidia/CudaToolkit/12.4
```

## Downloading model

First check the availability of models :

```bash
dorado download --list
```

```bash
[2026-07-03 11:35:01.472] [info] > simplex models
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_fast@v4.2.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_hac@v4.2.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_sup@v4.2.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_fast@v4.3.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_hac@v4.3.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_sup@v4.3.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_fast@v5.0.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_hac@v5.0.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_sup@v5.0.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_apk_sup@v5.0.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_fast@v5.2.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_hac@v5.2.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_sup@v5.2.0
[2026-07-03 11:35:01.472] [info]  - dna_r10.4.1_e8.2_400bps_hac@v6.0.0
[2026-07-03 11:35:01.472] [info]  - rna004_130bps_fast@v3.0.1
[2026-07-03 11:35:01.472] [info]  - rna004_130bps_hac@v3.0.1
[2026-07-03 11:35:01.472] [info]  - rna004_130bps_sup@v3.0.1
```

Here we are in simplex mode with dna and want the higheqt accuracy model
"sup" so we check for the most recent one :

Make a directory and then download the model of interest with the following

```bash
mkdir dorado_models
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.2.0 --directory dorado_models/
```

Basecalling :

- --recursive allows to treat all pod5 files given in the input folder even when it's in different folders
- --device cuda:all allows to make sure we use all available gpus when runnong the basecall
- --emit-fastq change the output from bam files to fastq files
- --kit-name the kit used to sequence the samples
- --emit-summary to have a summary for pycoQC
- if you want to treat separatly pass and fail add a --min-qscore for the treshold between pass and fail
- -o output folder where the arborescence of fastq files whlie be placed
- model : the model used for the basecalling (by default you can choose fast, hac or sup and it will download the corresponding one)
- but here we downloaded first the model with dorado download then used the downloaded model.
- input : the input folder containg all pod5 files we want to treat (as we use recursive the pod5 can be in lower folder
  e.g. input/flow1/pod5/file.pod5 input/flow2/pod5/file2/pod5)

```bash
dorado basecaller --recursive --device auto --emit-fastq --kit-name SQK-NBD114-24 --emit-summary -o output_dorado_sup/ dorado_models/dna_r10.4.1_e8.2_400bps_sup\@v5.2.0/ /home/bperez/save/MITI/Nanopore_run1
```
