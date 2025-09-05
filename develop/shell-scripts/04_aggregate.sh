#!/bin/bash

FILENAME="results/aggregate-summary.tsv"
echo -e "sample\tlines\twords" > "$FILENAME"

for sample in results/*.summary; do
    SAMPLE_NAME=$(basename "$sample" .summary)
    LINES=$(cat "$sample" | grep -e "^lines\t" | cut -f2)
    WORDS=$(cat "$sample" | grep -e "^words\t" | cut -f2)
    echo -e "$SAMPLE_NAME\t$LINES\t$WORDS" >> "$FILENAME"
done
