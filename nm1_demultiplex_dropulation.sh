#!/bin/bash

## Define Sample Name ##
# SAMPLE=GEX1_jeong  # Change this to run for a different sample
# SAMPLE=GEX2_han  # Change this to run for a different sample
# SAMPLE=GEX3_gil  # Change this to run for a different sample
SAMPLE=GEX4_stroke  # Change this to run for a different sample

## Paths to files ##
BASE_DIR=/data/scrna/scrna_snubh_batch20250106_counts
VCF=$BASE_DIR/vcf_merged/${SAMPLE}_final_chr_fixed.vcf.gz
NORM_VCF=$BASE_DIR/vcf_normalized/${SAMPLE}_normalized_sorted.vcf.gz
SAMPLE_DIR=$BASE_DIR/$SAMPLE/outs
COUNTS=$SAMPLE_DIR/filtered_feature_bc_matrix
BARCODES=$COUNTS/barcodes.tsv.gz
BAM=$SAMPLE_DIR/possorted_genome_bam.bam
INDS=$BASE_DIR/$SAMPLE/individuals.txt
CELL_TAG=CB
UMI_TAG=UB
GTF=/data/refs/gtf/gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz

## Dropulation output directory ##
DROPULATION_OUTDIR=$BASE_DIR/demuxafy/dropulation/$SAMPLE
mkdir -p $DROPULATION_OUTDIR

# singularity exec -B /data:/data Demuxafy.sif TagReadWithGeneFunction \
#           ANNOTATIONS_FILE=$GTF \
#           INPUT=$BAM \
#           OUTPUT=$DROPULATION_OUTDIR/${SAMPLE}_possorted_genome_bam_dropulation_tag.bam

# Define BAM file in DROPULATION_OUTDIR
BAM_DROPULATION=$DROPULATION_OUTDIR/${SAMPLE}_possorted_genome_bam_dropulation_tag.bam

## Ensure required files exist before running ##
for file in $NORM_VCF $BAM $BARCODES $INDS; do
    if [[ ! -f $file ]]; then
        echo "Error: Required file missing: $file" >&2
        exit 1
    fi
done

echo "All input files found, proceeding with analysis."

# echo "Running AssignCellsToSamples for $SAMPLE..."
# singularity exec -B /data:/data Demuxafy.sif AssignCellsToSamples -- \
#     --CELL_BC_FILE $BARCODES \
#     --INPUT_BAM $BAM_DROPULATION \
#     --OUTPUT $DROPULATION_OUTDIR/assignments.tsv.gz \
#     --VCF $NORM_VCF \
#     --SAMPLE_FILE $INDS \
#     --CELL_BARCODE_TAG 'CB' \
#     --MOLECULAR_BARCODE_TAG 'UB' \
#     --VCF_OUTPUT $DROPULATION_OUTDIR/assignment.vcf \
#     --MAX_ERROR_RATE 0.05 \

# if [[ $? -ne 0 ]]; then
#     echo "Error: AssignCellsToSamples failed!" >&2
#     exit 1
# fi
# echo "AssignCellsToSamples completed successfully."

# echo "Running DetectDoublets for $SAMPLE..."
# singularity exec -B /data:/data Demuxafy.sif DetectDoublets -- \
#     --CELL_BC_FILE $BARCODES \
#     --INPUT_BAM $BAM_DROPULATION \
#     --OUTPUT $DROPULATION_OUTDIR/likelihoods.tsv.gz \
#     --VCF $NORM_VCF \
#     --CELL_BARCODE_TAG $CELL_TAG \
#     --MOLECULAR_BARCODE_TAG $UMI_TAG \
#     --SINGLE_DONOR_LIKELIHOOD_FILE $DROPULATION_OUTDIR/assignments.tsv.gz \
#     --SAMPLE_FILE $INDS \
#     --MAX_ERROR_RATE 0.05

# if [[ $? -ne 0 ]]; then
#     echo "Error: DetectDoublets failed!" >&2
#     exit 1
# fi
# echo "DetectDoublets completed successfully."

echo "Running dropulation_call.R for $SAMPLE..."
singularity exec -B /data:/data Demuxafy.sif dropulation_call.R \
    --assign $DROPULATION_OUTDIR/assignments.tsv.gz \
    --doublet $DROPULATION_OUTDIR/likelihoods.tsv.gz \
    --out $DROPULATION_OUTDIR/updated_assignments.tsv.gz

echo "dropulation_call.R completed successfully."

echo "Pipeline completed successfully for $SAMPLE!"

# RUN!

# chmod +x $HOME/projects/nm1/nm1_demultiplex_dropulation.sh
# $HOME/projects/nm1/nm1_demultiplex_dropulation.sh
# $HOME/projects/nm1/nm1_demultiplex_dropulation.sh > $HOME/dropul_1.log 2>&1 echo $! > $HOME/dropul_1.pid
# ps -fp $(cat $HOME/dropul_1.pid)
# tail -f $HOME/dropul_1.log
