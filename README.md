# prepro
Snakemake pipeline for intitial processing (trimming, errorcorrection, merging) of Illumina data

## Typical usage

### Regular run on single node (assumes all data are in working directory)
```bash
snakemake --use-singularity -s Snakefile --singularity-args "-B $(pwd)"
```

### Distribute on cluster (VSC4)
```bash
snakemake -s Snakefile \
	--jobs 1000 --latency-wait 300 \
	--use-singularity --singularity-args "-B $BINFL -B /gpfs/data/fs71312/hahnc" \
	--cluster-config data/vsc4_config.yaml \
	--cluster '$(pwd)/bin/immediate_submit.py {dependencies} slurm' \
	--immediate-submit --notemp -pr
```

### Distribute on cluster (Sauron)
```bash
snakemake -s Snakefile \
	--jobs 1000 --latency-wait 300 \
	--use-singularity --singularity-args "-B /cl_tmp/hahnc" \
	--cluster-config data/sauron_config.yaml \
	--cluster "$(pwd)/bin/immediate_submit.py '{dependencies}' sge" \
	--immediate-submit --notemp -pr
```

## Prerequisites
### Set up conda and install Snakemake

To set us up with Snakemake we are going to use a package management system called `conda` (see [Documentation](https://docs.conda.io/en/latest/)).

First we need to setup `conda`. I will create a new directory called `conda`, then move to it, then download the installation script.
```bash
(user@host)-$ mkdir conda
(user@host)-$ cd conda
(user@host)-$ wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh
```

The installation script we just downloaded contains all instructions to set up Miniconda. Execute it (the `bash` part is telling the system that it should interpret the content in the Bash language). 
```bash
(user@host)-$ bash ./Miniconda3-py37_4.8.2-Linux-x86_64.sh
```

Now you will be prompted to agree to the licence - press enter, then space until you reach the end of the licence agreement (reading it along the way). At the bottom, you need to confirm with 'yes'. Then the installer asks you where you want to put the installation. For now you can just stick with the default and press enter. Then the installer will start working, downloading stuff, etc.. Might take a minute or two. Then you're asked if you want Miniconda to be initialized. Answer 'yes', and that should be almost it.

If you answered yes and you want to activate `conda` immediately in this session you can then do the following. If you're done for today you can just exit the session and `conda` will be available automatically for you once you start the next session.

Assuming you want to activate immediately.
```bash
(user@host)-$ . ~/.bashrc
```

If you install `conda` in this way it will by default be configured so that whenever you start a new shell you'll be setup in a virtual (base) conda environment - you may see your prompt change somewhat after you executed the above command. Nothing wrong with that in principle, but I don't want that so I disable this behaviour by executing the following command.
```bash
(base) (user@host)-$ conda config --set auto_activate_base false
```

If you are in the base environment at this point (note the `(base)` at the beginning of the prompt line), I suggest you deactivate the virtual env for now.
```bash
(base) (user@host)-$ conda deactivate
```

Now we want to set us up a virtual environment that has Snakemake installed using the `conda` package managing system. For reproducibility it's a good idea to put all the instructions for the setup in a configuration file. Conda expects this to come in the so-called `YAML` format. Let's use `nano` to make a new file.

```bash
(user@host)-$ nano snakemake_config.yml
```
Enter (or copy/paste) the following text to our config file. This will define the name of the environment to be created and the snakemake version you want to be installed. Conda will automatically fetch all packages that this particular snakemake version depends on. Within the conda universe software packages are more or less losely grouped in so-called channels and we'll ask conda to look for any software dependies of snakemake in these two channels. It's not so obvious why, but the order in which you list the channels may matter and in this case it's important to stick to this order (hats off to [@HannesOberreiter](https://github.com/HannesOberreiter) for calling that!)
```bash
name: snakemake
channels:
  - conda-forge
  - bioconda
dependencies:
  - snakemake==5.9.1
```
Save and close. 

Now, let's create the environment (this will take about 10 minutes in the test run).
```bash
(user@host)-$ conda env create -f snakemake_config.yml
```

After that, if you want to run Snakemake, you first need to enter the environment, then give it a whirl.
```bash
(user@host)-$ conda activate snakemake
(snakemake) (user@host)-$ snakemake -h
```
