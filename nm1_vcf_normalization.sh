#!/bin/bash

# Variables
SAMPLE=GEX1_jeong
BASE_DIR=/data/scrna/scrna_snubh_batch20250106_counts
VCF=$BASE_DIR/vcf_merged
VCF_NEW=$BASE_DIR/vcf_normalized

mkdir -p $VCF_NEW

# Step 1: Check REF mismatches (to detect lowercase bases)
echo "Checking reference mismatches..."
bcftools norm -c s -f /data/refs/refdata-gex-GRCh38-2024-A/fasta/genome.fa $VCF/${SAMPLE}_final_chr_fixed.vcf.gz 2>&1 | grep -c "does not match reference"

# Step 2: Convert REF/ALT bases to uppercase and save to a new file
echo "Converting REF and ALT bases to uppercase..."
zcat $VCF/${SAMPLE}_final_chr_fixed.vcf.gz | awk 'BEGIN{OFS="\t"} {if ($0 ~ /^#/) print; else {$4=toupper($4); $5=toupper($5); print}}' | bgzip > $VCF_NEW/${SAMPLE}_uppercase.vcf.gz
tabix -p vcf $VCF_NEW/${SAMPLE}_uppercase.vcf.gz

# Step 3: Remove variants with alternative alleles (or remove multi-allelic variants)
echo "Filtering out multi-allelic variants..."
bcftools view -m2 -M2 -v snps $VCF_NEW/${SAMPLE}_uppercase.vcf.gz -Oz -o $VCF_NEW/${SAMPLE}_filtered.vcf.gz
tabix -p vcf $VCF_NEW/${SAMPLE}_filtered.vcf.gz

# Step 4: Final normalization of the filtered VCF
echo "Performing final normalization..."
bcftools norm -m -both -f /data/refs/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
-Oz -o $VCF_NEW/${SAMPLE}_normalized_tmp.vcf.gz $VCF_NEW/${SAMPLE}_filtered.vcf.gz
mv $VCF_NEW/${SAMPLE}_normalized_tmp.vcf.gz $VCF_NEW/${SAMPLE}_normalized.vcf.gz
tabix -p vcf $VCF_NEW/${SAMPLE}_normalized.vcf.gz

# Run script in the background
# chmod +x $HOME/projects/nm1/nm1_vcf_normalization.sh
# nohup $HOME/projects/nm1/nm1_vcf_normalization.sh > norm_1.log 2>&1 &
# tail -f norm_1.log
# bcftools view /data/scrna/scrna_snubh_batch20250106_counts/vcf_normalized/GEX1_jeong_normalized.vcf.gz | head -n 80
# bcftools view /data/scrna/scrna_snubh_batch20250106_counts/vcf_normalized/GEX1_jeong_normalized_GQ_fixed.vcf.gz | head -n 80
# bcftools view /data/scrna/scrna_snubh_batch20250106_counts/vcf_normalized/GEX1_jeong_normalized_with_GQ.vcf.gz | head -n 80

# After fixing the GQ values: 
bgzip $VCF_NEW/${SAMPLE}_normalized_mod.vcf
tabix -p vcf $VCF_NEW/${SAMPLE}_normalized_mod.vcf.gz


# Step 5: 
## Sorting chromosomes

bcftools sort -Oz -o $VCF_NEW/${SAMPLE}_normalized_sorted.vcf.gz $VCF_NEW/${SAMPLE}_normalized_mod.vcf.gz
tabix -p vcf $VCF_NEW/${SAMPLE}_normalized_sorted.vcf.gz

