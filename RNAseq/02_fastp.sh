#!/usr/bin/env bash

# Input and output directories
INPUT_DIR="$1"
TRIMMED_DIR="$2"
REPORT_DIR="$3"

# Run fastp for each FASTQ file
for file in "$INPUT_DIR"/*.fastq.gz
do
    filename=$(basename "$file" .fastq.gz)

    fastp \
        -i "$file" \
        -o "$TRIMMED_DIR/${filename}_trim.fastq.gz" \
        -w 16 \
        -j "$REPORT_DIR/${filename}_fastp.json" \
        -h "$REPORT_DIR/${filename}_fastp.html"
done

# Summarize fastp reports
multiqc \
    "$REPORT_DIR" \
    -o "$REPORT_DIR/multiqc" \
    --filename fastp_multiqc