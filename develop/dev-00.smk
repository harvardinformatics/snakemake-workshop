rule count_lines:
    input: "data/sample1.txt"
    output: "results/sample1.lines"
    shell:
        "wc -l {input} | awk '{{print $1}}' > {output}"
