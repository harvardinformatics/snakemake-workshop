# Snakemake Workshop: Develop 

## Introduction

Welcome to the second part of our Snakemake workshop! In this section, we will learn how to begin creating a Snakemake workflow from a series of shell scripts. We will learn about what types of workflows might benefit from being converted to Snakemake, how to go from a series of shell scripts to a pipeline, and how to make a pipeline configurable so that it can be used with different datasets. 

## What workflows are suitable for Snakemake?

Although we love workflow managers like Snakemake, that doesn't mean that every single time we're running something on our computers we are wrapping them in a Snakemake script. A workflow manager has the following powerful features:

* **Resumability**: If a step in a workflow fails, you can rerun from that step.
* **Reproducibility**: You can rerun the same workflow with the same input data and get the same results.
* **Parallelization**: You can run multiple steps in parallel, and subsequent steps run without waiting for unrelated steps to finish.
* **Scalability**: Running the same workflow for 1 or 100 files is the same

However, there is a good amount of overhead involved in writing a Snakemake pipeline, so in order to take advantage of these features, your existing workflow should have the following characteristics so that it is worth the effort:

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

In practice, whether you want to directly write the command in your Snakemake file or call an external script is up to you. In this workshop, we will directly write the commands in the Snakemake file because they are short and easy to read. Directly writing the command in the `shell:` section has the benefit of making the workflow more self-contained and easier to understand. On the other hand, calling an external script is useful for running a more complex command. Note that if you do find yourself making a Snakemake workflow with few rules that all refer to complex scripts, it may be worth considering breaking up that complex script into multiple steps/rules. We will discuss this later in the section "What is the proper scope of a rule?". 

### Rule all

Having a rule does not do anything by itself. You will need a special rule called `all` that specifies the final output of your workflow. This is because Snakemake needs a way to know what the final product of your workflow is, so it knows what to work towards. Remember: snakemake is based on GNU Make, which is typically used to build software, and in that context, the `all` rule is like the final completed executable. 

When you begin to write a workflow, however, you build it up in steps. Rather than specifying the final For now, our workflow only has one output, so we can write that as our desired final output:

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

Let's leave everything hard coded for now so that we can focus on getting the workflow working. Soon, we will learn about how to use wildcards so that these rules can be generalized and take multiple files. 

### Testing your Snakemake workflow pt 1

The first thing you should do to test your Snakemake workflow is to run it with the `--dry-run` option. This will show you what Snakemake would do without actually running any commands. You can run this command in the terminal:

```bash
snakemake --snakefile develop/dev.smk --dry-run
```

A dry run is for checking the logic of your workflow. It does not read or modify any files. Its main function is to check each rule and find out if the required input files are available or will be generated by another rule. Make sure to read the output of your dry run to check that the logic matches your intention. 

Importantly, this dry run WILL NOT check whether the commands in the `shell:` directives are valid or will succeed. It just assumes that they will and only looks at the `input` and `output` sections. The dry run should look like this:

```
Building DAG of jobs...
Job stats:
job            count
-----------  -------
all                1
count_lines        1
total              2


[Tue Aug 19 15:59:10 2025]
rule count_lines:
    input: data/sample1.txt
    output: results/sample1.lines
    jobid: 1
    reason: Missing output files: results/sample1.lines
    resources: tmpdir=<TBD>
[Tue Aug 19 15:59:10 2025]
rule all:
    input: results/sample1.lines
    jobid: 0
    reason: Input files updated by another job: results/sample1.lines
    resources: tmpdir=<TBD>
Job stats:
job            count
-----------  -------
all                1
count_lines        1
total              2

Reasons:
    (check individual jobs above for details)
    input files updated by another job:
        all
    output files have to be generated:
        count_lines
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

<!-- explanation of the components of the dry run here -->

### Generalizing with wildcards

Now that we have a working rule, let's generalize it so that it can work with multiple input files. We will use wildcards to do this. Wildcards are placeholders that Snakemake uses to match file names. For example, we can change our `count_lines` rule to use a wildcard for the input file:

```python
rule count_lines:
    input: "data/{sample}.txt"
    output: "results/{sample}.lines"
    shell:
        "wc -l {input} | awk '{{print $1}}' > {output}"
```

The part of the file path that we are changing for the wildcards is the part where the filenames go `sample1`, `sample2`, etc. This is the variable part of the filename. However, just because we created a wildcard doesn't mean that Snakemake will scan your filesystem to look for files that match that pattern. The wildcards function only as a way to define the structure of the input and output file paths. In our case, we will need to provide a list of samples for Snakemake to process. 

```python
SAMPLES = ["sample1", "sample2"]
```

Now that we've listed the samples, how do we actually substitute the wildcards with those sample names? In Snakemake, there is a special function called `expand()` that we can use to fill in the wildcards with the actual sample names. For example, 

```python
expand("data/{sample}.txt", sample=['sample1', 'sample2'])
```

becomes

```python
["data/sample1.txt", "data/sample2.txt"]
```

Where will we put this `expand()` function? We can put it in the `input` section of rules. `expand` is rarely used in the output section of rules. In our simple Snakemake workflow, we only need to put the `expand()` function in the `input` section of final `all` rule. Our `dev-02.smk` should now look like this

```python
SAMPLES = ["sample1", "sample2"]

rule all:
    input: expand("results/{sample}.lines", sample=SAMPLES)

rule count_lines:
    input: "data/{sample}.txt"
    output: "results/{sample}.lines"
    shell:
        "wc -l {input} | awk '{{print $1}}' > {output}"
```

<!-- test again with dry-run and then run it for real, then delete the output folder-->

#### More advanced `expand`ing

You can have multiple wildcards in a single expand call. For example, you can expand both the sample name and the file type:

```python
expand("data/{sample}.{ext}", sample=['sample1', 'sample2'], ext=['txt', 'csv'])
```

becomes

```python
["data/sample1.txt", "data/sample2.txt", "data/sample1.csv", "data/sample2.csv"]
```

### Exercises 

In the next few sections, we will gradually build new rules for our Snakemake workflow based on the bash scripts in the `develop/shell-scripts` directory. We will name them `dev-01.smk`, `dev-02.smk`, etc. As we test and run our Snakemake workflows, we will be removing the `results` directory and its contents. This is because when we do a dry run, Snakemake checks to see if any intermediate or final files have already been created by the rules, and if so, excludes that rule from the dry run logic. This could cause unexpected behavior if we are intending to test the entire logic of the workflow. 

> **Exercise**: At this point, try to add another rule `count_words` that counts the number of words in each input file. Test your code with `--dry-run` and then run it for real to make sure it worked. Remember to remove the `results` directory before running the workflow again.

> **Solution**: Only showing new ormodified rules

```python
rule all:
    input: 
        expand("results/{sample}.{ext}", sample=SAMPLES, ext=["lines", "words"])

rule count_words:
    input: "data/{sample}.txt"
    output: "results/{sample}.words"
    shell:
        "wc -w {input} | awk '{{print $1}}' > {output}"
```

The next shell script we want to add to our Snakemake workflow is the `03_combine_counts.sh` script, which will combine the line and word counts for each sample into a single output file for each sample. This rule will take as an input the line and word count files for each sample. The final output will now be a `.summary` file for each sample. This means we will have to edit both `rule all` as well as add the new `rule combine_counts`. Additionally, we will have two inputs for the `combine_counts` rule: the line count file and the word count file. 

To add two inputs to a rule in Snakemake, you can specify them on separate lines and name them accordingly. Then, when you need to refer to those inputs, you can use the `input` keyword followed by the name of that input. This also applies for the output block. For example:

```python
rule example:
    input:
        A="{input}.A",
        B="{input}.B"
    output:
        C="{input}.C"
    shell:
        "cat {input.A} {input.B} > {output.C}"
```

Not every input needs to be a wildcard. You may have a rule that takes a fixed input file name and then a series of variable files. For example:

```python
rule example:
    input:
        A="fixed_input.A",
        B="{input}.B"
    output:
        C="{input}.C"
    shell:
        "cat {input.A} {input.B} > {output.C}"
```

> **Exercise**: Try making the necessary edits now in a new file called `dev-04.smk`. You will need to replace the input for `rule all` with the summary files and add a rule `combine_counts`.

**Solution**: Only showing modified rules.

```python
rule all:
    input: expand("results/{sample}.summary", sample=SAMPLES)

rule combine_counts:
    input:
        lines="results/{sample}.lines",
        words="results/{sample}.words"
    output:
        summary="results/{sample}.summary"
    shell:
        """
        echo -n "lines\t" > {output.summary}
        cat {input.lines} >> {output.summary}
        echo -n "words\t" >> {output.summary}
        cat {input.words} >> {output.summary}
        """
```

<!-- We should probably talk about multi-line shell scripts, tabs, and whitespace stuff here -->

### Aggregation/reduction rule

We will now move on to writing the final rule of our Snakemake workflow, which will adapt the script `04_aggregate.sh`. This script takes all the `.summary` files and combines them into a nice tab-separated file (TSV) with a header. Until now, we have only made rules that have a one-to-one correspondence between input and outputs, i.e. each input creates a unique output. This time, the `rule aggregate` will take many files and output a single file. Because this is a more complicated rule, let's break up the exercise into steps.

>**Exercise**: First step: what is the input to `rule aggregate` and how would you specify it?

**Solution**: The input to `rule aggregate` is all the `.summary` files. We know that all the summary files will be named following the pattern of `sample1.summary`, `sample2.summary`, etc. So we need to create an `expand()` function that will get all of these. 

```python
rule aggregate:
    input: expand("results/{sample}.summary", sample = SAMPLES)
```

>**Exercise**: Based on the shell script, what should the output be?

**Solution**: The output should be a single named file. Because it does not vary based on the input, it should not contain any wildcards. 

```python 
rule aggregate:
    input: expand("results/{sample}.summary", sample = SAMPLES)
    output: "results/combined.summary"
```

Now we will get to the difficult part of converting this script. First, notice how the input differs from previous rules' inputs. In our `rule combine_counts`, the input files were specified as separate lines for each count type (lines and words). So the input is a list of dictionaries like this:

```python
{
    "lines": "results/sample1.lines",
    "words": "results/sample1.words"
}
```

Then each of those files will get referenced in the shell command when we do `{input.lines}` and `{input.words}`. The shell interpreter then sees the input as `"results/sample1.lines"` and `"results/sample1.words"`, respectively. 

In contrast, when we use the `expand()` function in the `rule aggregate`, we create one big list of files that looks like this:

```python
["results/sample1.summary", "results/sample2.summary"]
```

When we refer to `{input}` in the shell command in this case, the shell interpreter will see the input as a single string with all the file names separated by spaces:

So the below shell directive 

```python
shell:
    cat {input} 
```

becomes

```bash
cat results/sample1.summary results/sample2.summary
```

Therfore, in order to work with individual files, we will need to use a loop to iterate over the list of input files. If you look at the `04_aggregate.sh` script, you will see that it already has the for loop written. Now we need to figure out where to put the `{input}`. 

>**Exercise**: Complete the `rule aggregate` and adapt the shell command. Save this has `dev-05.smk`

```python
rule aggregate:
    input:
        expand("results/{sample}.summary", sample = SAMPLES)
    output:
        "results/aggregate-summary.tsv"
    shell:
        """
        echo -e "sample\tlines\twords" > {output}
        for summary_file in {input}; do
            SAMPLE_NAME=$(basename "$summary_file" .summary)
            LINES=$(grep -e "^lines\t" "$summary_file" | cut -f2)
            WORDS=$(grep -e "^words\t" "$summary_file" | cut -f2)
            echo -e "$SAMPLE_NAME\t$LINES\t$WORDS" >> {output}
        done
        """
```

>**Exercise:** Modify the `rule all` to the correct final file and then test the workflow using dry-run. 

```python
rule all:
    input:
        "results/aggregate-summary.tsv"
```

### Other ways to specify inputs

### Decorations

## What is the proper scope of a rule?

