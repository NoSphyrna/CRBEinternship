# TaxInfo usages

Here you will find scripts of different usages of TaxInfo.

## Genotoul

In the folder run_on_genotoul, you can find an sbatch script "run_job.sh" that
can run the other scripts that are specific to a function of Taxinfo.

These scripts where designed to run on every file of the input file (# needs to be extended to all types of file not just)

Usage:

First set the config file according to your file organisation with your prefered text editor

```bash
vim config.sh
```

Then be sure that you are in the folder "run_on_genotoul"

And you can run the sbtach script

Example :

```bash
sbatch runjob.sh verify
```

This will launch the script gna_verifier_pq on all files of your input folder (this is also the default behaviour)
