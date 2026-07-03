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

```bash
dorado basecaller --recursive --device auto --emit-fastq --kit-name SQK-NBD114-24 --emit-summary -o output_dorado_sup/ dorado_models/dna_r10.4.1_e8.2_400bps_sup\@v5.2.0/ /home/bperez/save/MITI/Nanopore_run1
```
