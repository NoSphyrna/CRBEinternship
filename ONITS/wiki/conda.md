# Create a conda envirronment for running python scripts

## first load Miniconda

```bash
module load devel/Miniconda/Miniconda3
```

Then create a conda envirronment with biopython :

-p is to choose the path where we store it and the name (here biopy for exemple)

We then install python 3.11 and biopython
in this new conda envirronment

-y allows to automatically say yes when install ask for it

-c conda-forge --override-channels : fixes an issue with the cluster's conda
config by sepcifying explicitely the channel "conda-forge" for the downloading
and installation

```bash
conda create -p ~/work/conda/envs/biopy python=3.11 biopython -y -c conda-forge --override-channels
```

and if you are in interactive mode you can activate it with

```bash
conda activate ~/work/conda/envs/biopy
```

and in a sbatch job :

```bash
source activate ~/work/conda/envs/biopy
```

## ITSxRust

With bioconda, we need to put conda-forge first so it looks like :

```bash
conda create -p ~/work/conda/envs/itsxrust -c conda-forge -c bioconda itsxrust --override-channels -y
```

## p7zip env

```bash
conda create -p ~/work/conda/envs/p7zip -c conda-forge -c bioconda p7zip --override-channels -y
```

## HiTac

Here we use the solver libmamba to handle conflicts of differents versions

```bash
conda create -p ~/work/conda/envs/hitac -c conda-forge -c bioconda hitac --solver=libmamba -y
```
