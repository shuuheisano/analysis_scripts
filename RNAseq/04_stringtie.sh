#!/usr/bin/env bash

# Input and output paths
GTF_FILE="$1"
BAM_DIR="$2"
OUTPUT_DIR="$3"
PREPDE_PY3="$4"
READ_LENGTH="$5"

THREADS=8

# Quantify expression for each sample
for BAM_FILE in "$BAM_DIR"/*.sorted.bam
do
    SAMPLE=$(basename "$BAM_FILE" .sorted.bam)

    stringtie \
        -e \
        -p "$THREADS" \
        -G "$GTF_FILE" \
        -o "$OUTPUT_DIR/${SAMPLE}_stringtie.gtf" \
        -A "$OUTPUT_DIR/${SAMPLE}_gene_abundance.txt" \
        "$BAM_FILE"
done

# Create a list of sample IDs and StringTie GTF files
GTF_LIST="$OUTPUT_DIR/gtf_list.txt"
> "$GTF_LIST"

for SAMPLE_GTF in "$OUTPUT_DIR"/*_stringtie.gtf
do
    SAMPLE=$(basename "$SAMPLE_GTF" _stringtie.gtf)
    printf "%s\t%s\n" "$SAMPLE" "$SAMPLE_GTF" >> "$GTF_LIST"
done

# Generate gene- and transcript-level count matrices
python3 "$PREPDE_PY3" \
    -i "$GTF_LIST" \
    -l 100 \
    -g "$OUTPUT_DIR/gene_count_matrix.csv" \
    -t "$OUTPUT_DIR/transcript_count_matrix.csv"