    SAMPLES = ["sample1", "sample2"]

rule all:
    input: 
        "results/aggregate-summary.tsv"

def get_count_lines_input(wildcards):
    current_sample = wildcards.sample
    current_file = f"data/{current_sample}.txt"
    return current_file

rule count_lines:
    input: get_count_lines_input
    output: "results/{sample}.lines"
    shell:
        "wc -l {input} | awk '{{print $1}}' > {output}"

rule count_words:
    input: "data/{sample}.txt"
    output: "results/{sample}.words"
    shell:
        "wc -w {input} | awk '{{print $1}}' > {output}"

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

def get_aggregate_input(wildcards):
    my_samples = expand("results/{sample}.summary", sample = SAMPLES)
    return my_samples

rule aggregate:
    input:
        get_aggregate_input
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
