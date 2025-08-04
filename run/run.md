## Introduction

Welcome to today's workshop about the Snakemake workflow management software. This is day 1 of 2, where we'll be focusing on how to get pre-made workflows up and running. We'll touch on the basics of Snakemake's syntax and terminology, and learn how to debug some common problems. We'll also learn how to use Snakemake in conjunction with SLURM to really scale up your analyses.

But first, you'll need to download the workshop repository! You can do this in a couple of ways.

## Getting started

### Downloading the workshop

First, login to the cluster and navigate to the location in which you want to download the repo.

#### 1. If you are familiar with *git*:

If you'd like to use the `git` command line tool, you can clone the repo as follows:

```bash
git clone https://github.com/harvardinformatics/snakemake-workshop.git
```

#### 2. If you'd rather not use *git*:

You can also just download the repository archive directly and extract it:

```bash
wget https://github.com/harvardinformatics/snakemake-workshop/archive/refs/heads/main.zip
unzip snakemake-workshop
```

Either way, after you have the repo downloaded, enter the workshop directory:

```bash
cd snakemake-workshop
```

### Installing Snakemake

Snakemake is the workflow management software we'll be using today. We recommend installing it through the conda/mamba software managers (For more info about conda and mamba, see [our tutorial]()).

With `conda`, you should be able to setup an environment for Snakemake with:

```bash
conda create -n snakemake-env -c conda-forge snakemake-minimal
```

Confirm the prompts on the screen. Then when the environment is created, activate it:

```bash
conda activate snakemake-env
```

You should be able to type `snakemake` and see the message:

```
Error: No Snakefile found, tried Snakefile, snakefile, workflow/Snakefile, workflow/snakefile
```

This means snakemake is installed and ready to run! If you see any variation of the error message `command not found`, then something has gone wrong with your installation. Make sure your envrionment is activated (you should see `(snakemake-env)` at the front of your prompt), and if it still doesn't work let us know!

#### snakemake vs. snakemake-minimal

If you search for the snakemake package on conda-forge, you'll see several come up, with the two most popular being `snakemake` and `snakemake-minimal`. These will both install snakemake into your environment, however `snakemake` includes many extra dependencies, including those for remote execution and storage. Since our most likely use case is running on the cluster, we won't need those dependencies so we can stick with `snakemake-minimal`. But if you do ever need to store data or run workflows on places like AWS or Google, remember to use `snakemake`.

## Workflows

What do we mean when we say **workflow**? Well, let's picture some typical steps for a phylogenomic analysis:

1. Download genomes of interest and their annotations
2. Extract annotated regions of each genome (*e.g.* longest transcript from of each gene)
3. Predict orthology between transcripts of different species
4. Align orthologous sequences
5. Filter poorly aligned sequences
6. Create a phylogeny for each aligned gene
7. Assess rates of evolution across gene trees

Or:

![Example workflow diagram](../img/workflow-demo-01.drawio.png)

You can see that each of these steps requires input files, produces output files, and likely has an associate piece of software or a custom script written by the researcher. These steps and the tools associated with them are your analytical workflow - **the output of one step becomes the input of the next step**.

Importantly, each step needs to be done for each sample or gene. For steps 1 and 2, you may have to run whatever tool you use for every sample in your data set, maybe dozens or even hundreds of times. For steps 3-7, you may have to run whatever tool you use for every gene in your analysis, possibly thousands of times.

Let's take a look at a single step, (4) align orthologous sequences. If we use a common program, such as [MAFFT](), which aligns a single sequence at a time, what are our possible solutions to automate this for our potentially thousands of sequences?

### A bash script for every step

One of the most common way to automate analyses prior to workflow languages was the use of custom shell scripts. A shell is the program behind the text interface with which you interface with your computer, with `bash` being one of the most common shells. Shells themselves are like mini-programming languages, giving users the power to write code and commands to customize their environment. With that being said, one could write a **bash script** that loops over every unaligned sequence file in a directory, runs an alignment program on it, and saves the output:

```bash
#!/bin/bash

for locus in /path/to/loci/*.fasta
do
    locus_base=$(basename "$locus" .fasta)
    mafft "$locus" > /path/to/alignments/"$locus_base".aln
done
```

One might (and should) save this set of commands as a script for reproducibility, possibly calling it something like `04_run_mafft.sh`. That way when you look back at your analysis, you'll be able to easily remember what commands you ran and be able to run them again if needed.

Of course, this is just one step in our workflow. We would likely need a script for each step: `01_download_samples.sh`, `02_extract_genes.sh`, ... and so on. Then we'd have to run them each individually and deal with individual errors, or perhaps write a meta-script that runs each of them sequentially.

There are two problems with this approach:

1. High maintenance
2. Difficult to parallelize

### Job arrays

**Job arrays** can help with the parallelization problem. A job array is a feature of most modern **job schedulers** that are installed on institutional clusters. A job scheduler is a a program that handles user requests for resources and allocates compute nodes based on resources available. On the Cannon cluster at Harvard, we use the [SLURM]() job scheduler.

In a typical SLURM job, you would write a SLURM script to submit a request to run some command on a compute node:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00

mafft locus1.fa > alignments/locus1.aln
```

Here, the `#SBATCH` comments give information to the SLURM scheduler. **Job arrays** allow us to scale this up, running multiple inputs on the same command, in parallel if resources are available.

For a job array, in the context of aligning sequences, you would first create a text file that just lists the paths to the input sequences. This is sometimes called a **manifest** file:

```
loci/locus1.fasta
loci/locus2.fasta
loci/locus3.fasta
...
loci/locus1453.fasta
```

In this example, let's call this file `loci.txt`. Then, in your SLURM script, you specify the number of jobs you'd like to submit and edit the script so it uses


```bash
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --array=1-1453

# Get locus name for this array index as the line number in the file
locus=$(sed -n "${SLURM_ARRAY_TASK_ID}p" loci.txt)

locus_base=$(basename "$locus" .fasta)
mafft "$locus" > alignments/"$locus_base".aln
```

In this case, we specify that this job array will create 1453 tasks by using `#SBATCH --array=1-1453`. Then, in our script, we use a command called `sed` to pull the locus ID as the line number in the manifest file and use that ID as the basis for input and ouput from our alignment program. Then in the alignment program itself, instead of specifying an exact file name, we use the file name constructed from the locus ID. This will be constructed and submitted as a separate task for all 1453 loci.

Job arrays effectively solve the parallelization problem: as many of the tasks that there are resources available for will be submitted at once. However, you would still need to maintain a separate job array script for each step of the workflow.

### Workflow languages

Workflow languages try to lower the development and maintenance costs of compiling steps into a workflow while easily paralellizing tasks and even integrating with job schedulers. They also add infrastracture that allows one to start and stop your workflow, resume workflows if stopped or an error is encountered, and add samples to the workflow potentially without having to re-run all steps for all samples.

Today, of course, we'll be talking about [**Snakemake**](), but there are other workflow languages out there, including [Nextflow]() and [Common Workflow Language]().

## Snakemake

Snakemake is a scripting language built from Python (with a bit of YAML formatting mixed in). Its philosophy is based on the standard `make` tool used commonly when building software: there should be an end result of the script which should be a file or files. These are the **targets**, and the workflow builds from the target backwards. It is essentially asking itself, "What do I need to do to produce these target files?"

In our example of a phylogenomic analysis, the target would be the final `.csv` file with the results of the rate analysis. Snakemake would see that as the target and look backwards: Are all the gene trees that it was expecting to be there actually there? If not, then it will look backwards: Are all the alignments that it was expecting to be there actually there? And so on until it encounters a step where all the outputs exist, from which it would start moving forward through the pipeline.

### Rules

Snakemake and other workflow languages are made up of **rules**. A rule encompasses one step of your workflow. For example, a rule for our alignment step may look something like this:

```yaml
rule mafft_align:
    input:
        "data/{locus}.fasta"
    output:
        "alignments/{locus}.aln"
    shell:
        "mafft {input} > {output}"
```



## Preparing inputs

### wildcards

## the config file


## the dryrun


## rulegraphs and dags


## running a workflow


## running on a cluster

### resources

### profiles
