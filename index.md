---
title: "[Workshop] Introduction to Snakemake"
authors:
    - Gregg Thomas
    - Lei Ma
author_header: Workshop Developers
---

# Introduction to Snakemake

{{ author_row(page) }}

This is a two session workshop intended to introduce analysts to the workflow scripting language [Snakemake :octicons-link-external-24:](https://snakemake.readthedocs.io/en/stable/){ target="_blank" }. Workflow management systems, such as Snakemake, allow one to create reproducible and scalable data analysis pipelines. In the first session, we will cover the basics of **running** a Snakemake workflow. In the second session, we'll learn the basics of **writing and developing** a Snakemake workflow from scratch. This will aid those who wish to convert their analyses to workflows.

## Prerequisites

While we do not require attendees to have previous experience writing Snakemake scripts, this is an intermediate level workshop, so some computational experience is required:

- Basic knowledge of navigating a file system from the command line
- Previous experience with command line software and ability to read basic bash scripts
- Basic understanding software environments (*e.g.* conda) and containers (*e.g.* singularity)

!!! warning "Cluster account required"

    Additionally, since this workshop will involve exercises on the Cannon cluster, a [Cannon account through FASRC is required :octicons-link-external-24:](https://docs.rc.fas.harvard.edu/kb/quickstart-guide/){ target="_blank" }. Be sure you can login to the cluster **BEFORE CLASS**.

## Getting Started

This workshop will use files hosted on a [github repository :octicons-link-external-24:](https://github.com/harvardinformatics/snakemake-workshop){ target="_blank" } that you will need to download to the Cannon cluster.

### Downloading the workshop

First, login to the cluster and navigate to the location in which you want to download the repo. Since this is a temporary workshop, feel free to use the SCRATCH space (`/netscratch/YOUR_LAB/Lab/USERNAME`). 

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

Snakemake is the workflow management software we'll be using during the workshop. We recommend installing it through the conda/mamba software managers (For more info about conda and mamba, see [our tutorial](../../resources/tutorials/installing-command-line-software-conda-mamba.md)).

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

<!-- I think we should just have them install snakemake regular. Are we trying to save space here? -->

### Working with text files on the cluster

During this workshop we will be editing text files on the cluster, which means you must either be comfortable using a text editor on the command line, or you must be able to connect to the cluster with a remote desktop client such as [VSCode](), or edit the files on an OpenOnDemand instance.

If you want to use a text editor on the command line, we recommend using `nano`, which is a great beginner text editor. You open a file using the command `nano <filename>`. This will change your terminal to a text editor interface. You will need to use your arrow keys to navigate the text (not your mouse/clicking). You can type normally to add/change text. When you're done editing, you can save the file by pressing (on a PC) <!-- not sure -->, then `Enter`. On a Mac your keyboard shortcut for saving is `Ctrl + O`, then `Enter`. To exit, press <!-- something --> on a PC or `Ctrl + X` on a Mac, then `Enter`.

If you want to use a the text editor VSCode to edit files on the cluster, you can follow the instructions on the [FASRC Docs page](https://docs.rc.fas.harvard.edu/kb/vscode-remote-development-via-ssh-or-tunnel/) on how to set that up. You can either run VSCode in your browser, or you can install it locally and connect to the cluster. Give this a try, but if it doesn't work, just use `nano` for now and we can try to troubleshoot it later.

You can also open the files by using remote desktop or RStudio server. To start an OpenOnDemand session, you follow the directions on this [FASRC Docs page](https://docs.rc.fas.harvard.edu/kb/virtual-desktop/). <!-- more instructions here -->

## Workshop content

Part 1: Running a Snakemake workflow

[Snakemake Run](run/run.md){ .md-button }

Part 2: Developing a Snakemake workflow

[Snakemake Develop](run/run.md){ .md-button }