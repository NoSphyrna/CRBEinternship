# Download and format databases for each classifiers

## Download databases

- [UNITE](https://doi.plutof.ut.ee/doi/10.15156/BIO/3301246)
  (choose the desired version)

- [Eukaryome](https://eukaryome.org/sintax/)
  (again choose the desired version)

Put the compressed verison in a folder in the cluster.

To decompress the unite database it's quite straightforward :

example with the all eukaryotes version :

```bash
gzip -d utax_reference_dataset_all_19.02.2025.fasta.gz
```

For Eukaryome you 7zip but it's made for windows so on the cluster, you need to
use p7zip :

For this you can install it with bioconda : (see conda.md)

```bash
module load devel/Miniconda/Miniconda3

conda create -p ~/work/conda/envs/p7zip -c bioconda -c conda-forge p7zip --override-channels -y

conda activate ~/work/conda/envs/biopy
```

Then you can use 7zip as following:

```bash
7z e SINTAX_EUK_ITS_v2.0.fasta
```

Then you can clean the folder and keep only the wanted databases.

## Format databases

Set the folder configuration for the databases in the file **config_databases.cfg**

Then run the script **prepare_databases.sh** with

```bash
sbatch prepare_databases.sh
```
