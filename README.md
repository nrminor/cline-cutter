# Variant-calling Passerina WGS data with nf-core/sarek

## Overview

Bioinformatics has come a long way in the past few years. Containerized software and efficient, open-source workflow managers are now the norm, making it easier than ever to do science rigorously and reproducibly. Biologists have coalesced most widely around the workflow manager [Nextflow](https://nextflow.io/), which offers [many advantages for high-throughput biological data](https://www.nextflow.io/blog/2022/learn-nextflow-in-2022.html). These include that Nextflow workflows follow [FAIR Guiding Principles for scientific data management and stewardship](https://www.nature.com/articles/sdata201618). They are also portable, meaning they can be used in the same way across a wide variety of compute infrastructures (e.g., desktop personal computers, HPC Clusters, Amazon and Google cloud services, and Kubernetes). Nextflow handles all parallel execution itself, making it much more efficient than linear pipeline scripts, and also scalable to any amount of data.

Perhaps best of all, Nextflow boasts an active and growing support community, which, for many oft-used workflow applications, includes [nf-core](https://nf-co.re/). Pipelines in nf-core are among the most sophisticated and widely used pipelines to date, and are also completely open-source, well-documented, deeply configurable, and available to all via GitHub.

One such pipeline is [nf-core/sarek](https://nf-co.re/sarek), which is nf-core's standard variant-calling pipeline for diploid organisms (plus some add-ons for comparing tumor and normal cells for clinical bioinformatics). Sarek automates the process of breaking large sequence read files into small chunks and processing them in parallel, which massively reduces the compute time required to analyze large datasets. It also generates a wide variety of data QA/QC visualizations, allowing you to fine-tune your analysis to strengths and weaknesses of your data that you may otherwise be blind to.

For my thesis project, a key file format is the VCF. With the existence of nf-core/sarek, there truly is no reason to reinvent the bioinformatics wheel. Instead, I made a simple configuration file for running Sarek on the University of Wyoming Advanced Research Computing Center's HPC Cluster _Beartooth_ ([read more here](https://arccwiki.atlassian.net/wiki/spaces/DOCUMENTAT/pages/1683587073)). Below I detail the steps I took to configure and run Sarek on Beartooth. Later, I will add my post-processing analysis codebase to a second repo.

## Set-Up on the ARCC Beartooth Cluster

To proceed through the following three steps and reproduce my results, you will need:

1. Access to the Beartooth HPC cluster.
2. Access to resources via a project on Beartooth, such as _passerinagenome_.

Next, to get all your files in place, open a Terminal and `ssh` into beartooth, like so:

```
ssh <username>@beartooth.arcc.uwyo.edu
```

Then, change to your user folder within the associated project folder, create a directory where you will run this pipeline, and run:

```
module load arcc/1.0 gcc/12.2.0 git/2.37.0
https://github.com/nrminor/uwyo-thesis-project.git .
```

This will pull my configuration of files into the working directory.

### Step 1: Pre-run shell setup

Next, we need to set-up our shell so that it has the software it needs to run:

```
module load arcc/1.0 gcc/12.2.0 nextflow/22.10.4 singularity/3.10.3
currentdate=$(date +'%Y%m%d')
```

Also, be sure to download your reference genome of choice--mine was [GenBank Assembly GCA_014549065.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_014549065.1/)--and place it in the `resources/` folder.

### Step 2: Double check configuration

#### The Configuration File

This is arguably the most important step. I have specified most settings specific to Beartooth and my data in the critical `beartooth.config` file. Settings you should review for use with your own project include:

- The reference genome in FASTA format to be used. You specify the name of the genome to be used as shorthand as well as the file path to that FASTA.
- Beartooth run settings with the parameters `cluster_account`, `email`, and `date`.
- Where the input samplesheet, which specifies file paths and other important information related to your raw data, has been placed. By default I keep this in the `resources/` folder in the launch directory.
- Where to place the results (parameter `outdir`)
- Sequencing platform settings with parameter `seq_platform`
- Whether to variant call with paired reads only, in the case that you used Illumina paired-end sequencing (parameter `only_paired_variant_calling`).

If you are using Beartooth yourself, this config file also includes a number of important settings that you should _not_ change, including:

- To use slurm to execute tasks with `process.executor = 'slurm'`. Slurm is the workload manager and queueing infrastructure used at ARCC.
- To ignore certain settings that will confuse Beartooth with the parameter `schema_ignore_params`.

#### The Sample Sheet

The input file required by Sarek is a CSV spreadsheet with the following columns:

```
patient,sample,lane,fastq_1,fastq_2
```

The most important information in my case are the columns `sample`, `fastq_1`, and `fastq_2`, the last two of which are the absolute file paths to my paired read FASTQ files. This is how the pipeline finds your input data and pulls them into the workflow. My samplesheet, as an example, is in `resources/`.

### Step 3: Launch the pipeline

Finally, to run the pipeline with my configuration, I set up my working directory like so:

```
working_directory/
├── config/
│   └── beartooth.config
├── REAME.md
├── resources/
│   ├── GCA_014549065.1_CCar_1.0_genomic.fna
│   └── samplesheet.csv
```

This can be set up in any way if you alter the file paths in `beartooth.config`. Then, change into the working directory, and you're ready to proceed.

#### Running in the background

Usually, I use the Nextflow flag `-bg` to run the pipeline in the background, which is akin to using `screen` on Linux. This ensures that the pipeline will run to finish, even if your computer goes to sleep, your multi-factor authentication expires (this happens all the time on ARCC resources), the cluster has an outage, or some other disruption takes place.

To run the workflow in the background, use:

```
nextflow -bg run nf-core/sarek -r 3.0.1 -profile singularity -c config/beartooth.config

```

Once the pipeline gets going, Nextflow handles all resource requesting with Slurm, meaning you do not need to allocate resources with `sbatch` yourself.

If the workflow crashes, hangs up, or errors out, run the same command as above but with the `-resume` flag, like so:

```
nextflow -bg run nf-core/sarek -r 3.0.1 -profile singularity -c config/beartooth.config -resume
```

Thanks to Nextflow's detailed logging, the pipeline will be able to pick up where it left off instead of running all the way from the beginning.

#### Running interactively

If you keep running into errors or other issues, I recommend you run the pipeline interactively on a subset of your data (say, the first two rows only of the sample sheet). To do so, simply create a truncated copy of the sample sheet, run the same command as above with the new config file and without the `-bg` flag:

```
nextflow run nf-core/sarek -r 3.0.1 -profile singularity -c config/new.config
```

Resume the run after fixing any errors with:

```
nextflow run nf-core/sarek -r 3.0.1 -profile singularity -c config/new.config -resume
```
