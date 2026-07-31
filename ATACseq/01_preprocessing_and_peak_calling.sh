#!/bin/bash

############################################################
# 01_preprocessing_and_peak_calling.sh
#
# ATAC-seq preprocessing and peak-calling pipeline
#
# Workflow:
# 1. Raw-read quality control with FastQC and MultiQC
# 2. Nextera adapter trimming with Trimmomatic
# 3. Quality control of trimmed reads
# 4. Alignment to the hg38 reference genome with BWA-MEM
# 5. PCR duplicate removal and chromosome filtering with SAMtools
# 6. Peak calling with MACS2
#
# Input FASTQ naming convention:
#   <sample_id>_1.fastq.gz
#   <sample_id>_2.fastq.gz
############################################################

set -euo pipefail


############################################################
# User-defined paths and parameters
############################################################

# Root directory of the ATAC-seq project
PROJECT_DIR="/path/to/ATACseq_project"

# Input FASTQ directory
INPUT_DIR_FASTQ="${PROJECT_DIR}/data/raw/fastq"

# Metadata and scripts directories
METADATA_DIR="${PROJECT_DIR}/data/metadata"
SCRIPT_DIR="${PROJECT_DIR}/scripts"

# Output directories
FASTQC_DIR="${PROJECT_DIR}/data/processed/fastqc"
TRIM_DIR="${PROJECT_DIR}/data/processed/trimmomatic"
FASTQC_TRIM_DIR="${PROJECT_DIR}/data/processed/fastqc_trimmomatic"
BWA_DIR="${PROJECT_DIR}/data/processed/bwa"
MACS2_DIR="${PROJECT_DIR}/data/processed/macs2"

# Sample ID list
SAMPLE_LIST="${METADATA_DIR}/sample_ids.txt"

# hg38 reference genome FASTA file
REFERENCE_FASTA="/path/to/reference/hg38/genome.fa"

# BWA index prefix
BWA_INDEX_DIR="${METADATA_DIR}/bwa_index"
BWA_INDEX_PREFIX="${BWA_INDEX_DIR}/UCSC_hg38_index"

# Nextera adapter sequence file provided with Trimmomatic
NEXTERA_ADAPTERS="/path/to/trimmomatic/adapters/NexteraPE-PE.fa"

# Number of threads
THREADS=12


############################################################
# Create directories
############################################################

mkdir -p \
    "${METADATA_DIR}" \
    "${SCRIPT_DIR}" \
    "${FASTQC_DIR}" \
    "${TRIM_DIR}" \
    "${FASTQC_TRIM_DIR}" \
    "${BWA_DIR}" \
    "${MACS2_DIR}" \
    "${BWA_INDEX_DIR}"


############################################################
# 1. Quality control of raw FASTQ files
############################################################

# Check FastQC version
fastqc --version

# Run FastQC for all raw FASTQ files
fastqc \
    -t "${THREADS}" \
    -o "${FASTQC_DIR}" \
    --nogroup \
    "${INPUT_DIR_FASTQ}"/*.fastq.gz

# Check MultiQC version
multiqc --version

# Summarize raw-read FastQC results with MultiQC
mkdir -p "${FASTQC_DIR}/multiqc"

multiqc \
    "${FASTQC_DIR}" \
    --outdir "${FASTQC_DIR}/multiqc" \
    --filename fastqc_multiqc


############################################################
# 2. Generate sample ID list
############################################################

# Extract sample IDs from paired-end FASTQ filenames:
# <sample_id>_1.fastq.gz and <sample_id>_2.fastq.gz
ls "${INPUT_DIR_FASTQ}"/*_{1,2}.fastq.gz |
    sed -E 's/.*\/(.+)_([12])\.fastq\.gz/\1/' |
    sort -u > "${SAMPLE_LIST}"


############################################################
# 3. Nextera adapter trimming with Trimmomatic
############################################################

# Check Trimmomatic version
trimmomatic -version

# Create a separate trimming command script
TRIMMING_SCRIPT="${SCRIPT_DIR}/trimming.sh"

while read -r sample_id; do
    echo "trimmomatic PE -threads ${THREADS} \
${INPUT_DIR_FASTQ}/${sample_id}_1.fastq.gz \
${INPUT_DIR_FASTQ}/${sample_id}_2.fastq.gz \
${TRIM_DIR}/${sample_id}_tr_1P.fastq.gz \
${TRIM_DIR}/${sample_id}_tr_1U.fastq.gz \
${TRIM_DIR}/${sample_id}_tr_2P.fastq.gz \
${TRIM_DIR}/${sample_id}_tr_2U.fastq.gz \
ILLUMINACLIP:${NEXTERA_ADAPTERS}:2:30:10 \
MINLEN:30"
done < "${SAMPLE_LIST}" > "${TRIMMING_SCRIPT}"

# Run Trimmomatic
bash "${TRIMMING_SCRIPT}"


############################################################
# 4. Quality control of trimmed paired-end reads
############################################################

# Only paired reads are used for downstream alignment and
# peak-calling analyses. Therefore, FastQC is performed for:
#
#   *_tr_1P.fastq.gz
#   *_tr_2P.fastq.gz
#
# Unpaired reads are not included in downstream analysis.

fastqc \
    -t "${THREADS}" \
    --nogroup \
    "${TRIM_DIR}"/*_tr_1P.fastq.gz \
    "${TRIM_DIR}"/*_tr_2P.fastq.gz \
    -o "${FASTQC_TRIM_DIR}"

# Summarize trimmed-read FastQC results with MultiQC
mkdir -p "${FASTQC_TRIM_DIR}/multiqc"

multiqc \
    "${FASTQC_TRIM_DIR}" \
    --outdir "${FASTQC_TRIM_DIR}/multiqc" \
    --filename fastqc_trimmomatic_multiqc


############################################################
# 5. Build BWA index
############################################################

# Check BWA version
brew info bwa

# Build BWA index
bwa index \
    -p "${BWA_INDEX_PREFIX}" \
    "${REFERENCE_FASTA}"


############################################################
# 6. Alignment and BAM processing
############################################################

# Alignment of trimmed FASTQ files, PCR duplicate removal,
# and removal of noncanonical chromosomes
while read sample_id
do
    echo "Processing sample: $sample_id"

    # Define directories and file names
    BASE_DIR="/path/to/ATACseq_project/data/processed/bwa"
    TRIM_DIR="/path/to/ATACseq_project/data/processed/trimmomatic"
    INDEX="/path/to/ATACseq_project/data/metadata/bwa_index/UCSC_hg38_index"

    # Intermediate files
    NAME_SORT_BAM="${BASE_DIR}/${sample_id}_namesort.bam"
    FIXMATE_BAM="${BASE_DIR}/${sample_id}_fixmate.bam"
    SORTED_BAM="${BASE_DIR}/${sample_id}_sorted.bam"
    DEDUP_BAM="${BASE_DIR}/${sample_id}_dedup.bam"
    FINAL_BAM="${BASE_DIR}/${sample_id}.bam"
    STATS="${BASE_DIR}/${sample_id}_map_stats.txt"

    # Step 1: BWA alignment followed by name sorting
    bwa mem -t 12 ${INDEX} \
        ${TRIM_DIR}/${sample_id}_tr_1P.fastq.gz \
        ${TRIM_DIR}/${sample_id}_tr_2P.fastq.gz | \
        samtools sort -n -@12 -o ${NAME_SORT_BAM} -

    # Step 2: Add mate information using fixmate
    samtools fixmate -m -@12 ${NAME_SORT_BAM} ${FIXMATE_BAM}

    # Step 3: Sort by genomic coordinates
    samtools sort -@12 \
        -T "${BASE_DIR}/${sample_id}_tmp" \
        -o ${SORTED_BAM} \
        ${FIXMATE_BAM}

    # Step 4: Create BAM index
    samtools index ${SORTED_BAM}

    # Step 5: Remove PCR duplicates
    samtools markdup -@12 -r ${SORTED_BAM} ${DEDUP_BAM}

    # Step 6: Create BAM index
    samtools index ${DEDUP_BAM}

    # Step 7: Retain chr1–22, X, and Y only
    samtools view -@12 -b ${DEDUP_BAM} \
        $(samtools idxstats ${DEDUP_BAM} |
        cut -f 1 |
        grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)$') \
        > ${FINAL_BAM}

    # Step 8: Create index for the final BAM file
    samtools index ${FINAL_BAM}

    # Step 9: Output alignment statistics
    samtools flagstat ${FINAL_BAM} > ${STATS}

    # Remove intermediate files
    rm -f \
        ${NAME_SORT_BAM} \
        ${FIXMATE_BAM} \
        ${SORTED_BAM} \
        ${SORTED_BAM}.bai \
        ${DEDUP_BAM} \
        ${DEDUP_BAM}.bai

    echo "Finished processing: $sample_id"

done < /path/to/ATACseq_project/data/metadata/sample_ids.txt

############################################################
# 7. Peak calling with MACS2
############################################################

# Check MACS2 version
macs2 --version

while read -r sample
do
    echo "Processing peak calling for sample: ${sample}"

    # Final BAM:
    # duplicate-removed and restricted to canonical chromosomes
    INPUT_BAM="${BWA_DIR}/${sample}.bam"

    PREFIX="${sample}_peak"

    macs2 callpeak \
        -t "${INPUT_BAM}" \
        -f BAMPE \
        -g hs \
        -n "${PREFIX}" \
        -q 0.01 \
        --call-summits \
        -B \
        --SPMR \
        --outdir "${MACS2_DIR}"

    echo "Finished peak calling for sample: ${sample}"
    echo

done < "${SAMPLE_LIST}"

