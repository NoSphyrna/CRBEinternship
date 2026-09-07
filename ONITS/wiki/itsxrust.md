# Usage of ITSxRust

First install it with conda
to handle alls installs and dependancies
Then load ITSx
And find the HMM model

```bash
conda activate ~/work/conda/envs/itsxrust
module load bioinfo/ITSx/1.1.3
find $(dirname $(which ITSx))/../ -name "F.hmm" 2>/dev/null
```

On the cluster :

```bash
/usr/local/bioinfo/src/ITSx/ITSx_1.1.3/../ITSx_1.1.3/ITSx_db/HMMs/F.hmm
```

Then we can either give this path to ITSxRust or copy that into the desired folder
