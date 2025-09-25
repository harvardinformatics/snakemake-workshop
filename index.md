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
- Experience editing text files
- Basic understanding software environments (*e.g.* conda) and containers (*e.g.* singularity)

!!! warning "Cluster account required"

    Additionally, since this workshop will involve exercises on the Cannon cluster, a [Cannon account through FASRC is required :octicons-link-external-24:](https://docs.rc.fas.harvard.edu/kb/quickstart-guide/){ target="_blank" }. Be sure you can login to the cluster **BEFORE CLASS**.

## Getting Started

This workshop will use files hosted on a [github repository :octicons-link-external-24:](https://github.com/harvardinformatics/snakemake-workshop){ target="_blank" } that you will need to download to the Cannon cluster.

### Downloading the workshop

First, login to the cluster and navigate to the location in which you want to download the repo. 

#### 1. If you are familiar with `git`

If you'd like to use the `git` command line tool, you can clone the repo as follows:

```bash
git clone https://github.com/harvardinformatics/snakemake-workshop.git
cd snakemake-workshop
```

??? example "Code breakdown"

    | Code                                                                     | Description |
    | ------------------------------------------------------------------------ | ----------- |
    | `git clone https://github.com/harvardinformatics/snakemake-workshop.git` | This uses `git` to make a copy (clone) of our workshop directory, specified by the URL. |
    | `cd snakemake-workshop`                                                  | This changes your working directory to the newly cloned directory. |

#### 2. If you'd rather not use `git`

You can also just download the repository archive directly and extract it:

```bash
wget https://github.com/harvardinformatics/snakemake-workshop/archive/refs/heads/main.zip
unzip main.zip
cd snakemake-workshop-main
```

??? example "Code breakdown"

    | Code                                                                                        | Description |
    | ------------------------------------------------------------------------------------------- | ----------- |
    | `wget https://github.com/harvardinformatics/snakemake-workshop/archive/refs/heads/main.zip` | This uses the `wget` tool to download a compressed copy of our workshop repo, specified by the URL. |
    | `unzip main.zip`                                                                            | This decompresses the workshop directory.
    | `cd snakemake-workshop-main`                                                                | This changes your working directory to the newly downloaded workshop directory. |

### Installing Snakemake

Snakemake is the workflow management software we'll be using during the workshop. We recommend installing it through the conda/mamba software managers (For more info about conda and mamba, see [our tutorial :material-arrow-top-right:](../../resources/tutorials/installing-command-line-software-conda-mamba.md){ target="_blank" }).

#### Option 1: Loading conda if you don't already have it installed

If you don't already have conda or mamba installed, you will have to load the cluster modules. First, check if they are installed by typing either `conda` or `mamba` into the command line. If you see a help menu, you are good to continue. 

If you see any variation of the error `command not found`, you will have to load the modules with the following:

```bash
module load python
```

??? example "Code breakdown"

    | Code     | Description |
    | -------- | ----------- |
    | `module` | Run the cluster's module tool for pre-installed software. |
    | `load`   | The sub-command to run to load a given module. | 
    | `python` | The name of the module to load. |     

Then, if you check which modules are loaded with `module list` you should see:

```
Currently Loaded Modules:
  1) Miniforge3/25.3.1-fasrc01   2) python/3.12.11-fasrc01
```

Miniforge3 is the module that contains `conda`

#### Option 2: Using your own conda installation

If you already have conda/mamba installed and don't use the cluster's module, feel free to use that! Just follow the appropriate instructions below for activating your environments (or do it how you normally do it).

#### 2. Creating an environment with Snakemake

With `conda` loaded, you should be able to setup an environment for Snakemake with:

```bash
conda create -n snakemake-env -c bioconda snakemake-minimal -y 
```

??? example "Code breakdown"

    | Code                | Description |
    | ------------------- | ----------- |
    | `conda`             | Run the conda package and environment manager command line tool. |
    | `create`            | The sub-command to run from conda to create a new software environment. | 
    | `-n snakemake-env`  | Specifies the name of the new environment |
    | `-c bioconda`       | Specifies the conda channel from which to install packages into the environment while it is created. |
    | `snakemake-minimal` | The name of the packages to install while the environment is created. |
    | `-y`                | Automatically say 'yes' to any prompts during environment creation. |

With `-y` all prompts should be confirmed automatically. Then when the environment is created, activate it:

#### 3. Activating the environment

##### Option 1: Using the cluster module conda install

If you are using the cluster's conda module's, do the following to activate your environment:

```bash
source activate snakemake-env
```

??? example "Code breakdown"

    | Code            | Description |
    | --------------- | ----------- |
    | `source`        | Run the subsequent command as a bash script. |
    | `create`        | The sub-command to run to activate a previously created conda environment. | 
    | `snakemake-env` | The name of the environment to activate |

##### Option 2: Using your own conda install:

If you have your own conda installation, feel free to activate the environment as you normally do, but the following will likely work:

```bash
conda activate snakemake-env
```

??? example "Code breakdown"

    | Code            | Description |
    | --------------- | ----------- |
    | `conda`         | Run the conda package and environment manager command line tool. |
    | `create`        | The sub-command to run to activate a previously created conda environment. | 
    | `snakemake-env` | The name of the environment to activate |

##### Confirming that snakemake is installed

You should be able to type `snakemake` and see the message:

```
Error: No Snakefile found, tried Snakefile, snakefile, workflow/Snakefile, workflow/snakefile
```

This means snakemake is installed and ready to run! If you see any variation of the error message `command not found`, then something has gone wrong with your installation. Make sure your envrionment is activated (you should see `(snakemake-env)` at the front of your prompt), and if it still doesn't work let us know!

!!! note "snakemake vs. snakemake-minimal"

    If you search for the snakemake package on conda-forge, you'll see several come up, with the two most popular being `snakemake` and `snakemake-minimal`. These will both install snakemake into your environment, however `snakemake` includes many extra dependencies, including those for remote execution and storage. Since our most likely use case is running on the cluster, we won't need those dependencies so we can stick with `snakemake-minimal`. But if you do ever need to store data or run workflows on places like AWS or Google, remember to use `snakemake`.

### Working with text files on the cluster

During this workshop we will be writing scripts and editing text files on the cluster, which means you must either be comfortable using a text editor on the command line, or you must be able to connect to the cluster with a remote desktop client such as [VSCode :octicons-link-external-24:](https://code.visualstudio.com/){ target="_blank" }, or edit the files on an OpenOnDemand instance.

If you want to use a text editor on the command line, we recommend using `nano`, which is a great beginner text editor. You open a file using the command `nano <filename>`. This will change your terminal to a text editor interface. You will need to use your arrow keys to navigate the text (not your mouse/clicking). You can type normally to add/change text. When you're done editing, you can save the file by pressing `Ctrl + O`, then `Enter`. To exit, press `Ctrl + X`, then `Enter`. Look at the bottom of the interface for other options and messages.

If you want to use a the text editor VSCode to edit files on the cluster, you can follow the instructions on [this FASRC Docs page :octicons-link-external-24:](https://docs.rc.fas.harvard.edu/kb/vscode-remote-development-via-ssh-or-tunnel/){ target="_blank" } on how to set that up. You can either run VSCode in your browser, or you can install it locally and connect to the cluster. Give this a try, but if it doesn't work, just use `nano` for now and we can try to troubleshoot it later.

You can also open the files by using remote desktop or RStudio server. To start an OpenOnDemand (remote desktop) session, follow the directions on [this FASRC Docs page :octicons-link-external-24:](https://docs.rc.fas.harvard.edu/kb/virtual-desktop/){ target="_blank" }. Briefly, the steps are:

1. Connect to VPN using your FASRC credentials
2. Go to [https://rcood.rc.fas.harvard.edu :octicons-link-external-24:](https://rcood.rc.fas.harvard.edu){ target="_blank" }
3. Start a jupyterlab or R session with 4 GB memory, 1 CPU, and 4 hours 
4. Navigate to the directory where you want to put your repo
5. Use the terminal to install the snakemake conda environment & download the workshop repo as described above

## Workshop content

:material-calendar-clock: Links to content will appear as the date of each session approaches!

<!-- Click the buttons below to access the workshop content.

Part 1: Running a Snakemake workflow

[Snakemake Run](run/run.md){ .md-button .md-button--primary }

Part 2: Developing a Snakemake workflow

:construction_site: Under construction, check back soon! :construction: -->

<!-- [Snakemake Develop](run/run.md){ .md-button .md-button--primary } -->

---
