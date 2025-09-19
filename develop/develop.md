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

Create a `dev-01.smk` file in the `develop` directory. Write a rule that implements the `count_lines.sh` script. For now, you can hardcode the input and output filename (i.e. `dev_data/sample1.txt` and `results/sample1.lines`). You should have something like this:

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

Your full `dev-01.smk` file should now look like this:

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
snakemake --snakefile develop/dev-01.smk --dry-run
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

### Writing a `clean` rule to reset

Now that we'll be actively developing the workflow, it's a good idea to have a way to reset the state between runs. We can do this by writing a `clean` rule that removes all the output files. This will ensure that we start with a clean slate each time we run the workflow. Add the following rule to your `dev-01.smk` Snakefile:

```python
rule clean:
    shell: "rm -rf results/*"
```

Next time we want to run the workflow, we can first run `snakemake -s <snakefile> -R clean --cores 1`. This will remove all the output files and allow us to start fresh.

### Generalizing with wildcards

Now that we have a working rule, let's generalize it so that it can work with multiple input files. Create a `dev-02.smk` file in the `develop` directory and copy over the contents of `dev-01.smk`. Now, we will modify the rules to incorporate wildcards. Wildcards are placeholders that Snakemake uses to match file names. For example, we can change our `count_lines` rule to use a wildcard for the input file:

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

Where will we put this `expand()` function? We can put it in the `input` section of rules. `expand` is rarely used in the output section of rules. In our simple Snakemake workflow, we only need to put the `expand()` function in the `input` section of final `all` rule. Our `dev-02.smk` should now look like this:

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

Run `snakemake -s dev-02.smk -R clean --cores 1` to reset your project directory. Then do a dry run of `snakemake -s dev-02.smk --dry-run`. If that looks good, run it for real with `snakemake -s dev-02.smk --cores 1`.

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

In the next few sections, we will gradually build new rules for our Snakemake workflow based on the bash scripts in the `develop/shell-scripts` directory. We will name them `dev-03.smk`, `dev-04.smk`, etc. As we test and run our Snakemake workflows, we will be removing the `results` directory and its contents using the `clean` rule. This is because when we do a dry run, Snakemake checks to see if any intermediate or final files have already been created by the rules, and if so, excludes that rule from the dry run logic. This could cause unexpected behavior if we are intending to test the entire logic of the workflow.

> **Exercise**: Try to add another rule `count_words` that counts the number of words in each input file. This should do the same thing as the shell script `shell-scripts/02_count_words.sh`. When you add the new rule, make sure to also modify the `all` rule to reflect the new workflow logic and the final product. If you are lost, consult your workflow diagram that you made earlier. Test your code with `--dry-run` and then run it for real to make sure it worked. Remember to remove the `results` directory before running the workflow again.

> **Solution**: Only showing new or modified rules

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

The next shell script we want to add to our Snakemake workflow is the `03_combine_counts.sh` script, which will combine the line and word counts for each sample into a single output file for each sample. This rule will take as an input the line and word count files for each sample. The final output will now be a `.summary` file for each sample. As before, we will need to edit both `rule all` as well as add the new `rule combine_counts`. Additionally, we will have two inputs for the `combine_counts` rule: the line count file and the word count file. 

To add two inputs to a rule in Snakemake, you can specify them on separate lines, separated by a comma, and name them accordingly. Then, when you need to refer to those inputs, you can use the `input` keyword followed by the name of that input. This also applies for the output block. For example:

```python
rule example:
    input:
        A="{sample}.A",
        B="{sample}.B"
    output:
        C="{sample}.C"
    shell:
        "cat {input.A} {input.B} > {output.C}"
```

Not every input needs to be a wildcard. You may have a rule that takes a fixed input file name and then a series of variable files. For example:

```python
rule example:
    input:
        A="fixed_input.A",
        B="{sample}.B"
    output:
        C="{sample}.C"
    shell:
        "cat {input.A} {input.B} > {output.C}"
```

Another concept that is important in Snakemake formatting is multi-line shell commands. If you have a long shell command, you can split it across multiple lines for better readability. You can do this by using triple quotes `"""` to enclose the entire command. We will be using this feature in the next exercise. For example:

```python
rule example:
    input:
        A="{sample}.A",
        B="{sample}.B"
    output:
        C="{sample}.C"
    shell:
        """
        cat {input.A} {input.B} > {output.C}
        some_other_command here
        and so on
        """
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

### Aggregation/reduction rule

We will now move on to writing the final rule of our Snakemake workflow, which will adapt the script `04_aggregate.sh`. This script takes all the `.summary` files and combines them into a nice tab-separated file (TSV) with a header. Until now, we have only made rules that have a one-to-one correspondence between input and outputs, i.e. each input creates a unique output. This time, the `rule aggregate` will take many files and output a single file. Because this is a more complicated rule, let's break up the exercise into steps.

>**Exercise**: First step: what is the input to `rule aggregate` and how would you specify it using wildcards?

**Solution**: The input to `rule aggregate` is all the `.summary` files. We know that all the summary files will be named following the pattern of `sample1.summary`, `sample2.summary`, etc. So we need to create an `expand()` function that will get all of these. 

```python
rule aggregate:
    input: expand("results/{sample}.summary", sample = SAMPLES)
```

>**Exercise**: Second step: based on the shell script, what should the output be?

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

That means instead of one file being passed to the rule, we want a list of files being passed at the same time. When we refer to `{input}` in the shell command in this case, the shell interpreter will see the input as a single string with all the file names separated by spaces:

So the below shell directive 

```python
shell:
    cat {input} 
```

becomes

```bash
cat results/sample1.summary results/sample2.summary
```

Therfore, in order to work with individual files, we will need to use a loop to iterate over the list of input files. If you look at the `04_aggregate.sh` script, you will see that it already has the for loop written. Now we need to figure out where to put the `{input}` and `{output}`. 

To recap:

* One-to-one rules (like combine_counts) use wildcards to process one set of files per sample.
* Many-to-one rules (like aggregate) use an explicit input list (from expand()), combining all files in a single rule execution.

>**Exercise**: Third step: Complete the `rule aggregate` and adapt the shell command. Save this has `dev-05.smk`

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

>**Exercise:** Fourth step: Modify the `rule all` to the correct final file and then test the workflow using dry-run. 

```python
rule all:
    input:
        "results/aggregate-summary.tsv"
```

In the following sections, we will introduce some more advanced topics related to Snakemake that allow you to better customize how your workflow runs. It will be lighter on coding exercises because there is a lot to cover.

## Other ways to specify inputs

### Automatic file discovery using glob

So far, we've used wildcards and hardcoded sample names to specify input files.

Here's how we did the hard-coded sample name:

```python
rule count_lines:
    input: "data/sample1.txt"
```

Here's how we did the wildcard sample names:

```python
SAMPLES = ["sample1", "sample2"]
rule all:
    input: expand("results/{sample}.lines", sample=SAMPLES)
```

Real world workflows need more flexible ways to specify inputs. Let's do a more advanced form of finding files using `glob()` and `glob_wildcards()`. `glob` is a python function that uses pattern matching to find files and directories. It can use regular expressions or python wildcards (aka `*`). Below is an example of using glob to find all files that match the pattern `*.txt` in the data folder, and then extracting the sample name. This can be plugged directly into the `dev-05.smk` file in place of the first line.

```python
import glob
import os

SAMPLES = []
for txt_file in glob.glob("data/*.txt"):
    sample_name = os.path.basename(txt_file).replace(".txt", "")
    SAMPLES.append(sample_name)
```

Another way to glob files is to use the snakemake-specific glob function called `glob_wildcards()`. This one lets you do in one line what we did in a for loop above.

```python
SAMPLES = glob_wildcards("data/{sample}.txt").sample
```

Breaking it down:

1. glob_wildcards("data/{sample}.txt") searches for all files matching the pattern in the data directory.
2. the {sample} attribute extracts the sample names from the matched file paths.
3. glob_wildcards("data/{sample}.txt") returns a list of all file paths, while the .sample attribute returns a list of all parts of the file paths without the extension.

You can use this method to generate multiple sample names from a directory of files without having to list them all explicitly. You can glob more than one thing in your filename. For example:

```python
# Extract both sample name and condition from filenames
SAMPLES, CONDITIONS = glob_wildcards("data/{sample}_{condition}.txt")

# Create all combinations
rule all:
    input: 
        expand("results/{sample}_{condition}.lines", 
               zip, sample=SAMPLES, condition=CONDITIONS)

rule count_lines:
    input: "data/{sample}_{condition}.txt"
    output: "results/{sample}_{condition}.lines"
    shell:
        "wc -l {input} | awk '{{print $1}}' > {output}"
```

When you are writing your own glob patterns, it's important to double check what is happening by having print statements at each step. You can even create a python file just to practice globbing. For example, let's create a `glob_test.py` file with the following code to make sure that what we're globbing is correct:

```python
import glob
import os
from snakemake.io import glob_wildcards

# Test globbing
SAMPLES = glob_wildcards("data/{sample}.txt").sample
print(f"Found samples using glob_wildcards: {SAMPLES}")

SAMPLES = []
for txt_file in glob.glob("data/*.txt"):
    sample_name = os.path.basename(txt_file).replace(".txt", "")
    SAMPLES.append(sample_name)

print(f"Found samples using regular python glob: {SAMPLES}")
```

Then in your terminal, run `python glob_test.py` to double check that you are getting the right wild cards. If you have a more complicated file pattern that you want to match, LLMs are decent at generating regular expressions to capture the parts of paths that you are interested in. A combination of feeding the LLMs some sample file paths, the part you want to extract, and trial and error with a test script, is a good way to do more advanced globbing.

### Using samplesheets

While globbing for your files is effective, it does have some drawbacks.

1. It requires your files to be in an organized pattern that can be easily matched with wildcards, so it can be less flexible when dealing with complex directory structures or varying file names.
2. There is no record of the original file paths, which can make debugging more difficult. If you accidentally add a file that you didn't want, but that matches the glob, it'll have downstream effects.  
3. Files that don't match the glob are silently skipped over, which can lead to incomplete analysis if you're not careful. 

Another way to manage your inputs is by listing your samples in a samplesheet file and reading that in. Because snakemake is written in Python, you can use any Python code to read in your samplesheet file. To follow along, create a samplesheet called `samplesheet.txt` in the `develop` directory with the following contents:

```
sample1
sample2
```

Now, we can use base python to read that text file at the top of our snakefile. Create a new snakefile based off of `dev-05.smk` and call it `dev-samplesheet.smk`. Replace the first line `SAMPLES = ["sample1", "sample2"]` with:

```python
SAMPLES = []
with open("samplesheet.txt", "r") as f:
    SAMPLES = [line.strip() for line in f.readlines()]
```

What this does is read the `samplesheet.txt` file line by line, strips away whitespace (`line.strip()`), and stores the result in the `SAMPLES` list. This `SAMPLES` list is identical to `["sample1", "sample2"]`, except that it was generated dynamically from the contents of the `samplesheet.txt` file.

!!! Question
    What do you think would happen if we had sample3 in our `samplesheet.txt` file? Try it out!

For more complicated samplesheets, you can use the `pandas` library to read in a CSV or other tab-delimited file. This is useful if you have a collection of metadata associated with each sample, or a column of sample names plus a column of file paths. `pandas` comes pre-installed when you install snakemake. Below is an example CSV of forward and reverse reads and how it might be parsed into a workflow:

```
sample,R1,R2,condition
sample1,/path/to/data/sample1_R1.fastq,/path/to/data/sample1_R2.fastq,sample
sample2,/path/to/data/sample2_R1.fastq,/path/to/data/sample2_R2.fastq,sample
```

Now we read it using pandas and then parse out the sample column to a list. For the file paths, we create a dictionary where the key is the sample name, and the value is the R1 or R2 file path. Then, in the input section, we use a lambda function to dynamically retrieve the correct file paths for each sample.

```python
import pandas as pd

df = pd.read_csv("samplesheet.txt", sep="\t")
SAMPLES = df["sample"].tolist()
R1_FILES = dict(zip(df["sample"], df["R1"]))
R2_FILES = dict(zip(df["sample"], df["R2"]))

rule all:
    input:
        expand("results/{sample}.lines", sample=SAMPLES)

rule process_sample:
    input:
        R1=lambda wildcards: R1_FILES[wildcards.sample],
        R2=lambda wildcards: R2_FILES[wildcards.sample]
    output:
        "results/{sample}.lines"
    shell:
        "cat {input.R1} {input.R2} | wc -l > {output}"
```

### Using input functions

Sometimes you may want to use more complex logic to determine the input files for a rule. When you might have a multi-line specification to define a set of files, it can be useful to encapsulate that logic in a function. Snakemake allows you to define input functions that can take wildcards as arguments and return the appropriate file paths. For example, you might want to use logic to check that a file exists before adding it to the list. In the below code, we only want to get the fasta files if both the reverse and forward reads exist:

```python
# R1_FILES & R2_FILES are dictionaries generated elsewhere, perhaps by reading a samplesheet

def get_input_files(wildcards):
    sample = wildcards.sample
    r1 = R1_FILES.get(sample)
    r2 = R2_FILES.get(sample)
    if not os.path.exists(r1) or not os.path.exists(r2):
        raise ValueError(f"No input files found for sample {sample} in R1_FILES or R2_FILES.")
    return {"R1": r1, "R2": r2}

rule process_sample:
    input: get_input_files
    output: "results/{sample}.lines"
    shell:
        "cat {input.R1} {input.R2} | wc -l > {output}"
```

## Passing information to external scripts

So far we have just used shell commands in our rules. Sometimes, we might want to use external scripts in python, R, or other languages. To do this, we can use the `script:` directive in our rules in place of `shell:`. Within python or R scripts, you can access Snakemake input, output, and other objects using the `snakemake` object that will be passed to that script's environment. For example:

Here is what it might look like for writing a rule based on a python script

```python
rule run_python_script:
    input: "data.csv"
    output: "results.txt"
    script: "scripts/my_script.py"
```

```python
import pandas as pd

df = pd.read_csv(snakemake.input[0])
# Do some processing
df.to_csv(snakemake.output[0])
```

Similarly, this is how it might look inside an R script:

```R
library(dplyr)

df <- read.csv(snakemake@input[["data"]])
# Do some processing
write.csv(df, snakemake@output[["results"]])
```

In the following sections, we'll learn about more customization Snakemake can do. Things like params, configs, etc, are also passed to the external scripts through the `snakemake` object and can be accessed in the same way.

## Customizing how rules are run

In this section we will learn important concepts for controlling the execution of rules in Snakemake. Most bioinformatics workflows are complex and require parameters, compute resources, and software dependencies. Here is where we'll go over the ways to specify that in a rule. 

### Params

A `params` directive can be used to specify parameters for a rule. These parameters can be accessed within the shell command or script and can be used to customize the behavior of the rule. Params are typically used to define options for tools or scripts that are being called within the rule.

For example, if I wanted to change my shell script to take a parameter, it might look something like this:

```python
rule param_shell:
    input: "data.csv"
    output: "results.txt"
    params: param1="value"
    shell:
        "some_command --param1 {params.param1} {input} > {output}"
```

If I had an external R script that was running a function, I could pass parameters to it like this:

```python
rule run_r_script:
    input: "data.csv"
    output: "results.txt"
    params: param1="value"
    script: "scripts/my_script.R"
```

And my R script might look like this:

```R
param1 <- snakemake@params[["param1"]]
some_function(snakemake@input[["data"]], param1)
```

In the above, I gave my parameter the name `param1` and now I can access it by that name in the R script. You can have as many parameters as you need, separated by commas, and they can all be accessed via keyword within the script. 

Often, parameters are the `-` or `--` options in command-line tools. For example, in `bwa mem`, you have the option of specifying `-t` for number of threads and `-v` for verbosity, among others. The rule might look something like this:

```python
rule align_sequences:
    input: "raw_data/{sample}.fastq"
    output: "aligned/{sample}.bam"
    params:
        verbosity=4,
        num_threads=4,
    shell:
        "bwa mem -t {params.num_threads} -v {params.verbosity} reference.fasta {input} > {output}"

```

Because of the way we named these parameters in the rule, it actually makes the command line more readable. 

Importantly, just because we put the number of threads as a parameter, that doesn't mean that the rule will actually run with that many threads. Parameters are just keywords that you tell the tool to use; they don't enforce any behavior on their own. In the next section, we will see how to specify the actual compute resources that a rule should use. If we try to run this rule on a machine that doesn't have the required number of CPUs, it may run more slowly because bwa is expecting 4 CPUs. 

Params can also be global, meaning they can be defined outside of a rule and apply to the entire workflow. For example, you might want to have a param called `workdir` which points to where you want your working directory and all your files to be written. If your workflow has many parameters, it is a good idea to collect them all into a **config file**, which we will discuss later. 

### Resources

Snakemake has a built-in directive called `threads` that works for both resource allocation and passing that parameter to the shell. You can specify the number of threads a rule should use like this:

```python
rule my_rule:
    input: "data.csv"
    output: "results.txt"
    threads: 4
    shell:
        "some_command {input} -t {threads} > {output}"
```

This directive automatically sets the `cpus_per_task` to the number 4 and also passes it to the shell command. We can also reserve arbitrary resources like memory using the `resources` directive:

```python
rule my_rule:
    input: "data.csv"
    output: "results.txt"
    threads: 4
    resources:
        mem = 24 GB,
        runtime = 1 h
    shell:
        "some_command {input} -t {threads} > {output}"
```

To learn more about the resources directive, see the [snakemake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources).

### Software

Another important component of a workflow is software dependencies. So far, we have been using basic bash or python code (or pseudocode), but in practice, we will be using specific software tools and libraries. One of the benefits of workflow managers is that you can specify a totally separate software environment for each rule, including different versions of software or even different programming languages. This allows you to keep your rules modular and contained so that they can be easily reused or modified without affecting other parts of the workflow.

In Snakemake, you can specify what software to use for the rule using either the `conda` directive or the `container` directive.

#### Conda

You can use the `conda` directive in a few ways. The first way is to specify an existing conda environment that has already been created. You can reference the environment by name or by the path to the environment. 

```python
rule some_rule:
    input: "data.csv"
    output: "output.txt"
    conda: "my_env_name"
    script: "scripts/my_script.R"
```

```python
rule some_rule:
    input: "data.csv"
    output: "output.txt"
    conda: "/path/to/my_env_name"
    script: "scripts/my_script.R"
```

However, doing this means that you need to manage the conda environment yourself. This is an issue if you want to run this script on a different computer where the environment doesn't yet exist. Another way to specify a conda environment is to use a `yaml` file that contains the requirements for the environment you want the rule to run in. You can get the `yaml` file by writing it yourself, like:

```yaml
name: my_env_name
channels:
  - conda-forge
dependencies:
  - pandas=1.3.0
  - numpy=1.21.0
  - scikit-learn=0.24.2
```

You would then save that yaml file in your project directory. Perhaps under a folder named `envs`. And then you would use it in your Snakemake rule like this:

```python
rule some_rule:
    input: "data.csv"
    output: "output.txt"
    conda: "envs/my_env_name.yaml"
    script: "scripts/my_script.py"
```

What this does is causes Snakemake to create a temporary conda environment based on the specifications in the `yaml` file. This environment will be activated whenever the rule is executed, ensuring that all the required dependencies are available. This environment will be cached so that any other rule that uses it or any subsequent execution of the same rule can reuse it without having to recreate it. 

The major benefit of using this yaml file is that it will always travel with your project so that it essentially makes the environment documented and portable. By specifying the version numbers of each software package, you can better ensure that your workflow will have consistent behavior. 

#### Containers

Containers are a way of packaging all the requirements of a software into one image file. Containers are run by software called Docker or Singularity. You can usually find container images on dockerhub or the biocontainers registry. Here is an example of a rule using the software `mafft` from the biocontainers registry.

```python
rule mafft:
    input: "data.fasta"
    output: "data_aligned.fasta"
    container: "biocontainers/mafft:7.475--hdfd78af_0"
    shell:
        "mafft {input} > {output}"
```

<!-- can the full galaxyproject url be used for singularity?  e.g. https://depot.galaxyproject.org/singularity/mafft:7.525--h031d066_1 -->

Once you have defined your rule with the container directive, you can then run snakemake with the option `--sdm conda singularity`. (sdm stands for "software deployment method") Then, depending on whether you've used the `conda` or `container` directive in the rule, Snakemake will automatically create and manage the necessary environments for you.

!!! Tip
    In case you don't want to remember to run snakemake with the `--sdm conda singularity` option every time, you can create a configuration file (e.g., `config.yaml`) and specify the default software deployment method there. We will cover config files in more detail later.


### Exercise

**Exercise:** Let's go back to our `combine_counts` rule. Imagine that your collaborator decided that they wanted to combine counts in a different way, using the program `pandas`. Also, your collaborator wants to add an option to return a CSV instead of a TSV. Thankfully, your collaborator has provided you with a python file. The python file is below (you can also find it in the `scripts` folder). Additionally, you've been given a `yaml` file for a conda environment `combine_counts.yml`. And following that is the skeleton of the `combine_counts` rule. 

```python
# scripts/combine_counts.py
import pandas as pd

# Access Snakemake variables directly (no command line args needed!)
lines_file = snakemake.input.lines
words_file = snakemake.input.words
output_file = snakemake.output[0]
output_format = snakemake.params.output_format

# Read the count values
with open(lines_file, 'r') as f:
    line_count = int(f.read().strip())

with open(words_file, 'r') as f:
    word_count = int(f.read().strip())

# Create DataFrame
data = {
    'metric': ['lines', 'words'],
    'count': [line_count, word_count]
}
df = pd.DataFrame(data)

# Write output
if output_format == 'csv':
    df.to_csv(output_file, index=False)
else:  # tsv format
    df.to_csv(output_file, sep='\t', index=False)
```

As a refresher, this is the original combine_counts rule:

```python
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

Your task is to modify the `combine_counts` rule to do the following:

0.  create a new snakefile called `dev-customization.smk` and copy over the contents of `dev-05.smk`. Then we will work with the `combine_counts` rule for the next items. 
1.  modify the `shell` command to call the `combine_counts.py` script instead of using shell commands.
2.  add a parameter called `output_format` with a default value of `tsv`.
3.  add a conda environment specification to the `envs/combine_counts.yml` file.
4.  add a default resource of 1 thread and 1024 MB of memory

Solution

```python
rule combine_counts:
    input:
        lines="results/{sample}.lines",
        words="results/{sample}.words"
    output:
        "results/{sample}.summary"
    params: 
        output_format="tsv"
        threads: 1
        resources:
        mem_mb=1024
    conda:
        "envs/combine_counts.yml"
    script:
        "scripts/combine_counts.py"
```

### Logs

Our final topic on this section is the `logs` directive. This allows you to record log file for each individual rule. In fact, as we will see soon, it will allow you to record log files for each instance of a rule (aka each time a rule runs on a different piece of data). But first, let's discuss how to implement logging in our Snakemake workflow.

The first thing we need to do is to add the directive `log` to the rule. Then we need to specify a place for the log to go. Typically, this is a folder called "logs" that you create in your project directory. Importantly, this folder needs to be manually created before you run your snakemake workflow or else the program won't be able to find the path to save the log files. 

```python
rule aggregate:
    input:
        expand("results/{sample}.summary", sample = SAMPLES)
    output:
        "results/aggregate-summary.tsv"
    log:
        "logs/aggregate.log"
    shell:
        """
        echo -e "sample\tlines\twords" > {output}
        for summary_file in {input}; do
            SAMPLE_NAME=$(basename "$summary_file" .summary)
            LINES=$(grep -e "^lines\t" "$summary_file" | cut -f2)
            WORDS=$(grep -e "^words\t" "$summary_file" | cut -f2)
            echo -e "$SAMPLE_NAME\t$LINES\t$WORDS" >> {output}
        done &> {log}
        """
```

We have also modified the shell command. This is because just writing the log file is not enough. You need to make sure that the std out and std err streams of your command are properly directed to the file. std out is what is printed to the screen if you were to run this in a terminal window. std err is what is printed if there is an error. By adding `&> {log}` to the end of the command, you ensure that both streams are captured in the log file. 

Other ways of incorporating a log into your run is to use the `{log}` variable if your command line tool already has an option for writing log. For example:

```python
rule abc:
    input: "input.txt"
    output: "output.txt"
    log: "logs/abc.log"
    shell: "somecommand --log {log} {input} {output}"
```

If you are running a script, you can use the same redirect to write any print statements to the log:

```python
rule print_log:
    input: "input.txt"
    output: "output.txt"
    log: "logs/print_log.log"
    shell: "python scripts/print_log.py {input} {output} &> {log}"
```

Another way to do this is to use the built-in logging library in python or the logger package in R to write to the log files. So within your python or R script, you would have lines that are like print statements except they will be time-stamped and written to your choice of std out or std err. This gives you more fine-grained control over how your script records itself. 

In the above rule, `aggregate` only runs once, so we can just name the log file based on the rule. However, the other rules run multiple times, based on the number of samples. We can add wildcards to the log file name so that each instance of the rule gets its own log file.

```python
rule count_lines:
    input: "data/{sample}.txt"
    output: "results/{sample}.lines"
    log: "logs/{sample}_count_lines.log"
    shell:
        "wc -l {input} | awk '{{print $1}}' > {output}"
```

It is good practice to name your logs based on the rule, so that you know where it is coming from. In a more complex workflow, you can instead have directories for each rule's logs. Just make sure to create those directories first (perhaps write a rule or add python code to your snakefile to create all log directories!).

## Checkpointing

## Writing profiles & Config files
