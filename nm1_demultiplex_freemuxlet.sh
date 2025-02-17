#!/bin/bash

## Define Sample Name ##
SAMPLE="GEX1_jeong"  # Change this to run for a different sample

## Paths to files ##
BASE_DIR="/data/scrna/scrna_snubh_batch20250106_counts"
VCF="$BASE_DIR/vcf_merged/${SAMPLE}_final_chr_fixed.vcf"
SAMPLE_DIR="$BASE_DIR/$SAMPLE/outs"
COUNTS="$SAMPLE_DIR/filtered_feature_bc_matrix"
BARCODES="$COUNTS/barcodes.tsv"
BAM="$SAMPLE_DIR/possorted_genome_bam.bam"
INDS="$BASE_DIR/$SAMPLE/individuals.txt"
CELL_TAG="CB"
UMI_TAG="UB"
N=8

## Output directories ##
OUTDIR="$BASE_DIR/demuxafy"
FREEMUXLET_OUTDIR="$OUTDIR/freemuxlet"
PILEUP_OUTDIR="$FREEMUXLET_OUTDIR/pileup"

# Define output prefixes using sample name
PILEUP_PREFIX="$PILEUP_OUTDIR/${SAMPLE}_pileup"
FREEMUXLET_PREFIX="$FREEMUXLET_OUTDIR/${SAMPLE}_freemuxlet"

# Ensure output directories exist
mkdir -p "$OUTDIR"
mkdir -p "$PILEUP_OUTDIR"

# Check if necessary files exist before running
if [[ ! -f "$VCF" ]]; then
    echo "Error: VCF file not found: $VCF" >&2
    exit 1
fi
if [[ ! -f "$BAM" ]]; then
    echo "Error: BAM file not found: $BAM" >&2
    exit 1
fi
if [[ ! -f "$BARCODES" ]]; then
    echo "Error: Barcodes file not found: $BARCODES" >&2
    exit 1
fi

echo "All input files found, proceeding with analysis."

echo "Running PILEUP for $SAMPLE..."
singularity exec -B /data:/data Demuxafy.sif popscle_pileup.py \
    --sam "$BAM" \
    --vcf "$VCF" \
    --group-list "$BARCODES" \
    --tag-group "$CELL_TAG" \
    --tag-UMI "$UMI_TAG" \
    --out "$PILEUP_PREFIX"

if [[ $? -ne 0 ]]; then
    echo "Error: popscle_pileup.py failed!" >&2
    exit 1
fi
echo "PILEUP completed successfully."

echo "Running FREEMUXLET for $SAMPLE..."
singularity exec -B /data:/data Demuxafy.sif popscle freemuxlet \
    --plp "$PILEUP_PREFIX" \
    --out "$FREEMUXLET_PREFIX" \
    --tag-group "$CELL_TAG" \
    --tag-UMI "$UMI_TAG" \
    --group-list "$BARCODES" \
    --nsample "$N"

if [[ $? -ne 0 ]]; then
    echo "Error: popscle freemuxlet failed!" >&2
    exit 1
fi
echo "FREEMUXLET completed successfully."

echo "Pipeline completed successfully for $SAMPLE!"
