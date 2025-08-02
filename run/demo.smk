import yaml

configfile: "demo-config.yml"
samples = config["samples"]

rule all:
    input:
        "results/aggregate_summary.tsv"

rule count_lines:
    input:
        "demo-data/{sample}.txt"
    output:
        "results/{sample}.lines"
    shell:
        "wc -l {input} > {output}"

rule count_words:
    input:
        "demo-data/{sample}.txt"
    output:
        "results/{sample}.words"
    shell:
        "wc -w {input} > {output}"

rule combine_counts:
    input:
        lines="results/{sample}.lines",
        words="results/{sample}.words"
    output:
        "results/{sample}.summary"
    run:
        with open(input.lines) as lin, open(input.words) as wor, open(output[0], "w") as out:
            num_lines = lin.read().split()[0]
            num_words = wor.read().split()[0]
            out.write(f"lines\t{num_lines}\nwords\t{num_words}\n")

rule aggregate:
    input:
        expand("results/{sample}.summary", sample=samples)
    output:
        "results/aggregate_summary.tsv"
    run:
        with open(output[0], "w") as out:
            out.write("sample\tlines\twords\n")
            for sample, summary_file in zip(samples, input):
                with open(summary_file) as f:
                    lines = f.readlines()
                    num_lines = lines[0].strip().split("\t")[1]
                    num_words = lines[1].strip().split("\t")[1]
                    out.write(f"{sample}\t{num_lines}\t{num_words}\n")