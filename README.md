# Variant-calling Passerina WGS data with nf-core/sarek

### Overview

Bioinformatics has come a long way in the past few years. Containerized software and efficient, open-source workflow managers are now the norm, making it easier than ever to do science rigorously and reproducibly. Biologists have coalesced most widely around the workflow manager [Nextflow](https://nextflow.io/), which offers [many advantages for high-throughput biological data](https://www.nextflow.io/blog/2022/learn-nextflow-in-2022.html). These include that Nextflow workflows follow [FAIR Guiding Principles for scientific data management and stewardship](https://www.nature.com/articles/sdata201618). They are also portable, meaning they can be used in the same way across a wide variety of compute infrastructures (e.g., desktop personal computers, HPC Clusters, Amazon and Google cloud services, and Kubernetes). Nextflow handles all parallel execution itself, making it much more efficient than linear pipeline scripts, and also scalable to any amount of data.

Perhaps best of all, Nextflow boasts an active and growing support community, which, for many oft-used workflow applications, include [nf-core](https://nf-co.re/). Pipelines in nf-core are among the most sophisticated and widely used pipelines to date, and are also completely open-source, excellently documented, and available via GitHub.

One such pipeline is [nf-core/sarek](nf-core/sarek), which is nf-core's standard variant-calling pipeline for diploid organisms (plus some add-ons for comparing tumor and normal cells for clinical bioinformatics). Sarek automates the process of breaking large sequencing read files into small chunks and processing them in parallel, which massively reduces the amount of compute time required to analyze large datasets. It also generates a wide variety of data QA/QC visualizations, allowing you to fine-tune your analysis to strengths and weaknesses of your data that you may otherwise be blind to.

For my thesis project, a key file format is the VCF. With the existence of nf-core/sarek, there truly is no reason to reinvent the bioinformatics wheel. Instead, I made a simple configuration file for running SAREK on the University of Wyoming Advanced Research Computing Center's HPC Cluster _Beartooth_. Below I detail the steps I took to configure and run Sarek on Beartooth. Later, I will add my post-processing analysis codebase to a second repo.

### Step 1: Pre-run shell setup

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

The input file required by SAREK is a CSV spreadsheet with the following columns:

```
patient,sample,lane,fastq_1,fastq_2
```

The most important information in my case are the columns `sample`, `fastq_1`, and `fastq_2`, the last two of which are the absolute file paths to my paired read FASTQ files. This is how the pipeline finds your input data and pulls them into the workflow. My samplesheet, as an example, is in `resources/`.

### Step 3: Launch the pipeline

To run the pipeline with my configuration, I set up my working directory like so:

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
nextflow -bg run nf-core/sarek -profile singularity -c config/beartooth.config

```

If the workflow crashes, hangs up, or errors out, run the same command as above but with the `-resume` flag, like so:

```
nextflow -bg run nf-core/sarek -profile singularity -c config/beartooth.config -resume
```

#### Running interactively

If you keep running into errors or other issues, I recommend you run the pipeline interactively on a subset of your data. To do so, simply run the same command as above but without the `-bg` flag:

```
nextflow run nf-core/sarek -profile singularity -c config/beartooth.config
```

Resume the run after fixing any errors with:

```
nextflow run nf-core/sarek -profile singularity -c config/beartooth.config -resume
```
