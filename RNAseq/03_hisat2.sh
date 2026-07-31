#!/usr/bin/env bash

# Input and output paths
GENOME_FASTA="$1"
INDEX_PREFIX="$2"
INPUT_DIR="$3"
SAM_DIR="$4"
BAM_DIR="$5"

THREADS=8

# Build the HISAT2 genome index
hisat2-build \
    -p "$THREADS" \
    "$GENOME_FASTA" \
    "$INDEX_PREFIX"

# Align trimmed single-end reads to the reference genome
for FILE in "$INPUT_DIR"/*.fastq.gz
do
    SAMPLE=$(basename "$FILE" .fastq.gz)

    hisat2 \
        -p "$THREADS" \
        --dta \
        -x "$INDEX_PREFIX" \
        -U "$FILE" \
        -S "$SAM_DIR/${SAMPLE}.sam" \
        --summary-file "$SAM_DIR/${SAMPLE}_summary.txt"
done

# Convert SAM files to coordinate-sorted BAM files
for SAM_FILE in "$SAM_DIR"/*.sam
do
    SAMPLE=$(basename "$SAM_FILE" .sam)

    samtools sort \
        -@ "$THREADS" \
        -o "$BAM_DIR/${SAMPLE}.sorted.bam" \
        "$SAM_FILE"
done