SAMPLES = ["sample1", "sample2"]

rule all:
    input: expand("results/{sample}.lines", sample=SAMPLES)

rule count_lines:
    input: "data/{sample}.txt"
    output: "results/{sample}.lines"
    shell:
        "wc -l {input} | awk '{{print $1}}' > {output}"

rule clean:
    shell: "rm -r results/*"
