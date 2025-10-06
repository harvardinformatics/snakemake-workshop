import os

SAMPLES = config["sample_list"]
INPUT_DIR = config["input_dir"]
OUTPUT_DIR = config["output_dir"]

print(OUTPUT_DIR)

rule all:
    input: 
        os.path.join(OUTPUT_DIR, "aggregate-summary.tsv")

rule count_lines:
    input: os.path.join(INPUT_DIR, "{sample}.txt")
    output: os.path.join(OUTPUT_DIR, "{sample}.lines")
    shell:
        "wc -l {input} | awk '{{print $1}}' > {output}"

rule count_words:
    input: os.path.join(INPUT_DIR, "{sample}.txt")
    output: os.path.join(OUTPUT_DIR, "{sample}.words")
    shell:
        "wc -w {input} | awk '{{print $1}}' > {output}"

rule combine_counts:
    input:
        lines=os.path.join(OUTPUT_DIR, "{sample}.lines"),
        words=os.path.join(OUTPUT_DIR, "{sample}.words")
    output:
        summary=os.path.join(OUTPUT_DIR, "{sample}.summary")
    shell:
        """
        echo -n "lines\t" > {output.summary}
        cat {input.lines} >> {output.summary}
        echo -n "words\t" >> {output.summary}
        cat {input.words} >> {output.summary}
        """

rule aggregate:
    input:
        expand(os.path.join(OUTPUT_DIR, "{sample}.summary"), sample = SAMPLES)
    output:
        os.path.join(OUTPUT_DIR, "aggregate-summary.tsv")
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
