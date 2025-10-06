SAMPLES = ["sample1", "sample2"]

rule all:
    input:
        "results/aggregate-summary.tsv"

rule count_lines:
    input:
        "data/{sample}.txt"
    output:
        data = "results/{sample}.lines"
    shell:
        "wc -l {input} | awk '{{print $1}}' > {output}"

rule count_words:
    input:
        "data/{sample}.txt"
    output:
        "results/{sample}.words"
    shell:
        "wc -w {input} | awk '{{print $1}}' > {output}"

checkpoint perform_qc:
    input:
        lines = expand("results/{sample}.lines", sample = SAMPLES),
        words = expand("results/{sample}.words", sample = SAMPLES)
    output:
        "results/qc-manifest.txt"
    run:
        with open(output[0], "w") as manifest:
            for sample in SAMPLES:
                lines_file = f"results/{sample}.lines"
                words_file = f"results/{sample}.words"

                num_words = int(open(words_file).read().strip())
                if num_words <= 7:
                    manifest.write(f"{sample}\n")

rule combine_counts:
    input:
        lines="results/{sample}.lines",
        words="results/{sample}.words",
        manifest="results/qc-manifest.txt"
    output:
        "results/{sample}.summary"
    shell:
        """
        echo -n "lines\t" > {output}
        cat {input.lines} >> {output}
        echo -n "words\t" >> {output}
        cat {input.words} >> {output}
        """

def get_aggregate_input(wildcards):
    manifest = checkpoints.perform_qc.get().output[0] # Wait for checkpoint and get output path
    with open(manifest) as mf:
        samples = [line.strip() for line in mf if line.strip()]
    return [f"results/{sample}.summary" for sample in samples]

rule aggregate:
    input: get_aggregate_input
    output: "results/aggregate-summary.tsv"
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
