# Snakemake Workshop: Develop 

## Introduction

Welcome to the second part of our Snakemake workshop! In this section, we will learn how to begin creating a Snakemake workflow from a series of shell scripts. We will learn about what types of workflows might benefit from being converted to Snakemake, how to go from a series of shell scripts to a pipeline, and how to make a pipeline configurable so that it can be used with different datasets. 

## What workflows are suitable for Snakemake?

Although we love workflow managers like Snakemake, that doesn't mean that every since time we're running something on our computers we are wrapping them in a Snakemakw script. A workflow manager has the following powerful features:

* **Resumability**: If a step in a workflow fails, you can rerun from that step.
* **Reproducibility**: You can rerun the same workflow with the same input data and get the same results.
* **Parallelization**: You can run multiple steps in parallel, and subsequent steps run without waiting for unrelated steps to finish.
* **Scalability**: Running the same workflow for 1 or 100 files is the same

However, there is a good amount of overhead involved in writing a Snakemake pipeline, so in order to take advantge of these features, your existing workflow should have the following characteristics so that it is worth the effort:

* **Multiple, interrelated steps**: If you have multiple file transformations that depend on each other, this is a good candidate for Snakemake. Don't use Snakemake if you have just one or two simple steps that are easily launched in a shell script.
* **Need for reruns**: If you need to rerun the same steps with different parameters (e.g. for benchmarking), or if you need to compare the results of different runs, Snakemake could help make this easier and more reproducible.
* **Can benefit from simple parallelization**: You have many independent steps that can be run in parallel, such as multiple input files to be processed with the same steps
* **Complex dependencies across steps**: Each rule in Snakemake can run in a different environment, so this would be useful if you want to compartmentalize your software environments

A classic example of a workflow that would NOT be worth converting into a Snakemake workflow is some one-off exploratory or data cleaning script that you aren't sure you will ever run again. Another example is a super convoluted workflow that you don't fully understand or didn't write yourself, and don't have time to debug. In that case it might be better to just start fresh rather than trying to convert. 

## Converting shell scripts to Snakemake

### Drawing the input-output of your workflow

Read through the shell script files in the `shell-scripts` directory. These scripts take a set of text files in the `data` directory, count the number of lines and words in each file, and then combine those counts into a summary file for each input file.

Draw a diagram of the input-output relationships of these scripts. For example, you might draw something like this:

```
data/sample1.txt  -->  count_lines.sh  -->  results/sample1.lines
```

### Writing the first Snakemake rule

Create a `dev.smk` file in the `develop` directory. Write a rule that implements the `count_lines.sh` script. For now, you can hardcode the input and output filename (i.e. `dev_data/sample1.txt` and `results/sample1.lines`). You should have something like this:

```python
rule count_lines:
    input: "data/sample1.txt"
    output: "results/sample1.lines"
    shell: 
        "bash shell-scripts/01_count_lines.sh {input}"
```

What if we didn't want to rely on the shell script? We could directly put the bash command in the `shell` directive:

```python
rule count_lines:
    input: "data/sample1.txt"
    output: "results/sample1.lines"
    shell: 
        "wc -l {input} | awk '{{print $1}}' > {output}"
```

!! Note
    The `shell` directive allows you to run any shell command, and you can use `{input}` and `{output}` to refer to the input and output files of the rule. The double curly braces `{{ }}` are used to escape the curly braces in the `awk` command, so that Snakemake does not interpret them as placeholders.

In practice, whether you want to directly write the command in your Snakemake file or call an external script is up to you. In this workshop, we will directly write the commands in the Snakemake file because they are short and easy to read. 

### Rule all

Having a rule does not do anything by itself. You will nee a special rule called `all` that specifies the final output of your workflow. For now, our one rule workflow only has one output, so we can write that as our desired final output:

```python
rule all:
    input: "results/sample1.lines"
```

Your full `dev.smk` file should now look like this:

```python
rule all:
    input: "results/sample1.lines"

rule count_lines:
    input: "data/sample1.txt"
    output: "results/sample1.lines"
    shell:
        "wc -l {input} | awk '{{print $1}}' > {output}"
``` 

### Testing your Snakemake workflow pt 1

The first thing you should do to test your Snakemake workflow is to run it with the `--dry-run` option. This will show you what Snakemake would do without actually running any commands. You can run this command in the terminal:

```bash
snakemake --snakefile develop/dev.smk --dry-run
```

### Generalizing with wildcards

### Adding more rules

### Decorations

## What is the proper scope of a rule?

