# scripts/combine_counts.py
import pandas as pd

# Access Snakemake variables directly (no command line args needed!)
lines_file = snakemake.input.lines
words_file = snakemake.input.words
output_file = snakemake.output[0]
output_format = snakemake.params.output_format

# Read the count values
with open(lines_file, 'r') as f:
    line_count = int(f.read().strip())

with open(words_file, 'r') as f:
    word_count = int(f.read().strip())

# Create DataFrame
data = {
    'metric': ['lines', 'words'],
    'count': [line_count, word_count]
}
df = pd.DataFrame(data)

# Write output
if output_format == 'csv':
    df.to_csv(output_file, index=False)
else:  # tsv format
    df.to_csv(output_file, sep='\t', index=False)
