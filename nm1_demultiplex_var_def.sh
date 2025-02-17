## Paths to files ##
VCF=/data/scrna/scrna_snubh_batch20250106_counts/vcf_merged/GEX1_jeong_final_chr_fixed.vcf
COUNTS=/data/scrna/scrna_snubh_batch20250106_counts/GEX1_jeong/outs/filtered_feature_bc_matrix/
BARCODES=$COUNTS/barcodes.tsv
BAM=/data/scrna/scrna_snubh_batch20250106_counts/GEX1_jeong/outs/possorted_genome_bam.bam
INDS=/data/scrna/scrna_snubh_batch20250106_counts/GEX1_jeong/1_individuals.txt
# GTF=/path/to/genes.gtf ## We do Not provide this - it should be the gtf file that you used to align your data. Otherwise you can download an appropriate gtf file from https://www.gencodegenes.org/human/

## Output directories ##
OUTDIR=/data/scrna/scrna_snubh_batch20250106_counts/demuxafy
DEMUXALOT_OUTDIR=$OUTDIR/demuxalot
DROPULATION_OUTDIR=$OUTDIR/dropulation
DOUBLETDETECTION_OUTDIR=$OUTDIR/DoubletDetection
SCDBLFINDER_OUTDIR=$OUTDIR/scDblFinder
SCDS_OUTDIR=$OUTDIR/scds

# Demuxafy: 
singularity exec --bind /data/scrna/ Demuxafy.sif python Demuxalot.py \
        -b $BARCODES \
        -a $BAM \
        -n $INDS \
        -v $VCF \
        -o $DEMUXALOT_OUTDIR \
        -r True