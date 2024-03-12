# Installation on Flamingo, the GR's computing cluster

## Installation of pipeline
No installation of pipeline is required. Everything is downloaded, unzip and installed for you.

## Installation of working environment
You need a working environment with:
- python (version 3.8.5)
- snakemake (version 5.32.2)
- singularity (version 3.6.3)

It is possible to create this environment with conda by:
```
conda env create -f /mnt/beegfs/pipelines/single-cell/envs/conda/single-cell_user.yaml --prefix="<path/MyNameEnvir>"
```
NB: singularity is not in this conda environment because of a namespace issue for slurm, but singularity is available by `module load`.

To activate complete environnement:
```
source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate <path/MyNameEnvir>
module load singularity
```
# Manual installation on your local machine

## Installation of pipeline
To download this workflow you can use 2 methods:
1. install git (see: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
2. `git clone https://github.com/mAGLAVE/single-cell.git`
3. go in single-cell/resources/WHITELISTS/ and unzip the 3M-february-2018.txt.gz file.
4. built singularity environment (singularity can be installed thanks to the "Installation of working environment" next section):
```
        singularity build single_cell.simg single_cell.def
        singularity build single_cell_TCR_BCR.simg single_cell_TCR_BCR.def
        singularity build single_cell_oldcerebro.simg single_cell_oldcerebro.def
```

Or:
1. download this workflow
![clone the pipeline](https://github.com/mAGLAVE/single-cell/blob/master/images/clone_pipeline.png).
2. unzip the dowloaded file
3. go in single-cell/resources/WHITELISTS/ and unzip the 3M-february-2018.txt.gz file.
4. built singularity environment (singularity can be installed thanks to the "Installation of working environment" next section):
```
        singularity build single_cell.simg single_cell.def
        singularity build single_cell_TCR_BCR.simg single_cell_TCR_BCR.def
        singularity build single_cell_oldcerebro.simg single_cell_oldcerebro.def
```

## Installation of working environment
You need a working environment with:
- python (version 3.8.5)
- snakemake (version 5.32.2)
- singularity (version 3.6.3)

It is possible to create this environment with conda:
1. install conda (see: https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
2. install the environment by: `conda env create -f /mnt/beegfs/pipelines/single-cell/envs/conda/single-cell_user_local.yaml -n <MyNameEnvir>`

To activate complete environnement:
```
conda activate <MyNameEnvir>
```