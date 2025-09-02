import glob
import os
from snakemake.io import glob_wildcards

# Test globbing
SAMPLES = glob_wildcards("data/{sample}.txt").sample
print(f"Found samples using glob_wildcards: {SAMPLES}")

SAMPLES = []
for txt_file in glob.glob("data/*.txt"):
    sample_name = os.path.basename(txt_file).replace(".txt", "")
    SAMPLES.append(sample_name)

print(f"Found samples using regular python glob: {SAMPLES}")
