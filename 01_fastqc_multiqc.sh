#!/usr/bin/env bash

# Input and output directories
INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Run FastQC
fastqc \
    -t 12 \
    --nogroup \
    -o "$OUTPUT_DIR" \
    "$INPUT_DIR"/*.fastq.gz

# Summarize FastQC reports
multiqc \
    "$OUTPUT_DIR" \
    -o "$OUTPUT_DIR"/multiqc \
    --filename fastqc_multiqc