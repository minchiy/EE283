#!/bin/bash
#SBATCH --job-name=Filter_SNPs
#SBATCH --output=filter_snps.out
#SBATCH --error=filter_snps.err
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

VCF_DIR="/pub/minchiy/EE283/Week5/SNPs"

echo "Filtering SNPs where two strains are 0/0 and two are 1/1..."

awk '{sum=0; for (i=2; i<=5; i++) sum+=$i; if (sum==4) print $0}' ${VCF_DIR}/X_1Mb.012 > ${VCF_DIR}/X_1Mb_filtered.012

echo "Filtered SNPs saved in X_1Mb_filtered.012"

