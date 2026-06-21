#!/bin/bash
# =============================================================================
# qc_reads.sh — Quality Control for 16S rRNA Amplicon Reads
# =============================================================================
# Description: Runs FastQC and fastp on paired-end amplicon FASTQ files
# Usage: bash scripts/qc_reads.sh <R1.fastq.gz> <R2.fastq.gz> <sample_name>
# Output: results/qc/<sample_name>/
# Author: Taruna Gupta | github.com/tarunagupta20
# =============================================================================

set -euo pipefail

# ── Arguments ─────────────────────────────────────────────────────────────────
R1=$1
R2=$2
SAMPLE=$3
OUTDIR="results/qc/${SAMPLE}"

# ── Validate inputs ────────────────────────────────────────────────────────────
if [ $# -ne 3 ]; then
    echo "Usage: bash scripts/qc_reads.sh <R1.fastq.gz> <R2.fastq.gz> <sample_name>"
    exit 1
fi

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "ERROR: Input files not found: $R1 / $R2"
    exit 1
fi

# ── Setup ──────────────────────────────────────────────────────────────────────
mkdir -p "${OUTDIR}/fastqc_raw"
mkdir -p "${OUTDIR}/fastqc_trimmed"
mkdir -p "data/processed"

echo "============================================"
echo "  16S QC Pipeline"
echo "  Sample:  ${SAMPLE}"
echo "  R1:      ${R1}"
echo "  R2:      ${R2}"
echo "  Output:  ${OUTDIR}"
echo "============================================"

# ── STEP 1: FastQC on raw reads ────────────────────────────────────────────────
echo "[STEP 1] FastQC on raw reads..."
fastqc "${R1}" "${R2}" \
    --outdir "${OUTDIR}/fastqc_raw" \
    --threads 4

# ── STEP 2: Trim with fastp ────────────────────────────────────────────────────
echo "[STEP 2] Trimming with fastp..."

R1_TRIM="data/processed/${SAMPLE}_R1_trimmed.fastq.gz"
R2_TRIM="data/processed/${SAMPLE}_R2_trimmed.fastq.gz"

fastp \
    --in1 "${R1}" --in2 "${R2}" \
    --out1 "${R1_TRIM}" --out2 "${R2_TRIM}" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 200 \
    --json "${OUTDIR}/${SAMPLE}_fastp.json" \
    --html "${OUTDIR}/${SAMPLE}_fastp.html" \
    --thread 4

# ── STEP 3: FastQC on trimmed reads ───────────────────────────────────────────
echo "[STEP 3] FastQC on trimmed reads..."
fastqc "${R1_TRIM}" "${R2_TRIM}" \
    --outdir "${OUTDIR}/fastqc_trimmed" \
    --threads 4

# ── Done ───────────────────────────────────────────────────────────────────────
echo "============================================"
echo "  QC COMPLETE"
echo "  FastQC raw:     ${OUTDIR}/fastqc_raw/"
echo "  FastQC trimmed: ${OUTDIR}/fastqc_trimmed/"
echo "  fastp report:   ${OUTDIR}/${SAMPLE}_fastp.html"
echo "  Trimmed reads:  data/processed/"
echo "============================================"
