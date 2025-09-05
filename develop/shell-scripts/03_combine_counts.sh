#!/bin/bash

SAMPLE_NAME=$1
SAMPLE_LINES="results/${SAMPLE_NAME}.lines"
SAMPLE_WORDS="results/${SAMPLE_NAME}.words"

echo -n "lines	" > results/$SAMPLE_NAME.summary
cat $SAMPLE_LINES >> results/$SAMPLE_NAME.summary
echo -n "words	" >> results/$SAMPLE_NAME.summary
cat $SAMPLE_WORDS >> results/$SAMPLE_NAME.summary
