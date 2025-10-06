# Input function to discover splits.
def get_split_outputs(wildcards):
    manifest = checkpoints.split.get().output[1]
    with open(manifest) as f:
        splits = [line.strip() for line in f if line.strip()]
    return [f"splits/{split}.txt" for split in splits]

rule all:
    input: get_split_outputs

# Checkpoint splits the file and writes a manifest.
checkpoint split:
    input:
        "complete/split-data.txt"
    output:
        directory("splits"),
        "manifest/splits.txt"
    run:
        import os

        os.makedirs("splits", exist_ok=True)
        splits = []
        with open(input[0]) as infile:
            lines = [line.strip() for line in infile if line.strip()]
            for l in lines:
                outfile = f"splits/{l}.txt"
                with open(outfile, "w") as out:
                    out.write(l + "\n")
                splits.append(l)
        os.makedirs("manifest", exist_ok=True)
        with open(output[1], "w") as mf:
            for s in splits:
                mf.write(s + "\n")

