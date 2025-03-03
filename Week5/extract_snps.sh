#!/bin/bash
#SBATCH --job-name=Extract_SNPs
#SBATCH --output=extract_snps.out
#SBATCH --error=extract_snps.err
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32

module load bcftools/1.15.1
module load vcftools/0.1.16

VCF_DIR="/pub/minchiy/EE283/Week5/SNPs"
FILTERED_VCF="${VCF_DIR}/allsamples.final.vcf.gz"

echo "Extracting SNPs from X chromosome (first 1Mb)..."

# 建立索引
bcftools index ${FILTERED_VCF}

# 提取 X 染色體前 1Mb 的 SNPs
bcftools view -r X:1-1000000 -O z -o ${VCF_DIR}/X_1Mb.vcf.gz ${FILTERED_VCF}

# 轉換為 012 格式並移除 INDEL
vcftools --gzvcf ${VCF_DIR}/X_1Mb.vcf.gz --remove-indels --012 --out ${VCF_DIR}/X_1Mb

echo "SNP extraction completed!"
