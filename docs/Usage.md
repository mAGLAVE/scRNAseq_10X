# Usage on Flamingo, the GR's computing cluster

1. make the parameters file according to your needs (configfile: see Parameters files section)
2. indicate the path to this file in the path_to_configfile variable
3. run the snakemake command
```
#parameters
path_to_configfile="<path/to/your_configfile.yaml>"
path_to_pipeline="/mnt/beegfs/pipelines/single-cell"

#launch
snakemake --profile ${path_to_pipeline}/profiles/slurm -s ${path_to_pipeline}/Snakefile --configfile ${path_to_configfile}
```
NB: The first utilisation can take some time because of the installation of conda sub-environment (automatic).

# Usage on your local machine

1. make the parameters file according to your needs (configfile: see Parameters files section)
2. indicate the path to this file in the path_to_configfile variable
3. indicate the path to this pipeline in the path_to_pipeline variable
4. run the snakemake command
```
#parameters
path_to_configfile="<path/to/your_configfile.yaml>"
path_to_pipeline="<path/to/single-cell>"

#launch
snakemake --profile ${path_to_pipeline}/profiles/local -s ${path_to_pipeline}/Snakefile --configfile ${path_to_configfile}
```
NB: The first utilisation can take some time because of the installation of conda sub-environment (automatic).