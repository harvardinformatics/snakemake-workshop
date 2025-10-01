SAMPLES = ["sample1", "sample2"]

rule all:
    input: 
        "results/aggregate-summary.tsv"

rule count_lines:
    input: "data/{sample}.txt"
    output: "results/{sample}.lines"
    shell:
        "wc -l {input} | awk '{{print $1}}' > {output}"

rule count_words:
    input: "data/{sample}.txt"
    output: "results/{sample}.words"
    shell:
        "wc -w {input} | awk '{{print $1}}' > {output}"

def get_combine_counts_input(wildcards):
    current_sample = wildcards.sample
    lines_file = f"results/{current_sample}.lines"
    words_file = f"results/{current_sample}.words"
    return [lines_file, words_file]

# def get_combine_counts_input(wildcards):
#     current_sample = wildcards.sample
#     lines_file = f"results/{current_sample}.lines"
#     words_file = f"results/{current_sample}.words"
#     return { "lines": lines_file, "words": words_file }
## alternate solution using dictionaries and unpack

rule combine_counts:
    input:
        get_combine_counts_input
        # unpack(get_combine_counts_input)  ## alternate solution using dictionaries and unpack
    output:
        summary="results/{sample}.summary"
    shell:
        """
        echo -n "lines\t" > {output.summary}
        cat {input[0]} >> {output.summary}
        echo -n "words\t" >> {output.summary}
        cat {input[1]} >> {output.summary}
        """

    # shell:
    #     """
    #     echo -n "lines\t" > {output.summary}
    #     cat {input.lines} >> {output.summary}
    #     echo -n "words\t" >> {output.summary}
    #     cat {input.words} >> {output.summary}
    #     """
    ## alternate solution using dictionaries and unpack

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

rule clean:
    shell: "rm -r results/*"
