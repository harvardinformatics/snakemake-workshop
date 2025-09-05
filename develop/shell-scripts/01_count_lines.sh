#!/bin/bash

SAMPLE=$1
SAMPLE_NAME=$(basename "$SAMPLE" .txt)

wc -l "$SAMPLE" | awk '{print $1}' > results/"$SAMPLE_NAME".lines
