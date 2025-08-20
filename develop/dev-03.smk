SAMPLES = ["sample1", "sample2"]

rule all:
    input: 
        expand("results/{sample}.{ext}", sample=SAMPLES, ext=["lines", "words"])

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