import os

sample_sheet = config["sample_sheet"]
with open(sample_sheet) as f:
    samples = [line.strip() for line in f if line.strip()]

INDIR = config["input_directory"]
OUTDIR = config.get("output_directory", "results")

rule all:
    input:
        os.path.join(OUTDIR, "aggregate-summary.tsv")

rule count_lines:
    input:
        os.path.join(INDIR, "{sample}.txt")
    output:
        os.path.join(OUTDIR, "{sample}.lines")
    shell:
        "wc -l {input} > {output}"

rule count_words:
    input:
        os.path.join(INDIR, "{sample}.txt")
    output:
        os.path.join(OUTDIR, "{sample}.words")
    log:
        os.path.join(OUTDIR, "{sample}.words.log")
    shell:
        "wc-demo-error -w {input} > {output} 2> {log}"

rule combine_counts:
    input:
        lines = os.path.join(OUTDIR, "{sample}.lines"),
        words = os.path.join(OUTDIR, "{sample}.words")
    output:
        os.path.join(OUTDIR, "{sample}.summary")
    run:
        with open(input.lines) as lin, open(input.words) as wor, open(output[0], "w") as out:
            num_lines = lin.read().split()[0]
            num_words = wor.read().split()[0]
            out.write(f"lines\t{num_lines}\nwords\t{num_words}\n")

rule aggregate:
    input:
        expand(os.path.join(OUTDIR, "{sample}.summary"), sample=samples)
    output:
        os.path.join(OUTDIR, "aggregate-summary.tsv")
    run:
        with open(output[0], "w") as out:
            out.write("sample\tlines\twords\n")
            for sample, summary_file in zip(samples, input):
                with open(summary_file) as f:
                    lines = f.readlines()
                    num_lines = lines[0].strip().split("\t")[1]
                    num_words = lines[1].strip().split("\t")[1]
                    out.write(f"{sample}\t{num_lines}\t{num_words}\n")