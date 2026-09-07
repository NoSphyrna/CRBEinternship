# ONITS

A pipeline for metabarcoding with Oxford Nanopore Technologies over ITS (but not only)

## Description

This pipeline allows to treat and filter reads from raw [pod5](https://software-docs.nanoporetech.com/output-specifications/latest/read_formats/pod5/)
datas to filtered and curated OTUs that are then assigned to a taxonomy (optional part)

## Steps

- **Basecalling** : Takes pod5 files with raw datas and returns sequences in fastQ files (demultiplexed from Nanopore barcodes not personnalised ones)
- **Quality stats** : Uses pycoQC to get quality infos over the sequencing
- **Demultiplexing/Trimming** :
