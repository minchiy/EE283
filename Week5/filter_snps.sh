q#!/bin/bash
#SBATCH --job-name=SNP_filtering
#SBATCH --output=snp_filtering.out
#SBATCH --error=snp_filtering.err
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=35

# 加載模組
module load bcftools/1.15.1

VCF_DIR="/pub/minchiy/EE283/Week5/SNPs"
RAW_VCF="${VCF_DIR}/allsamples.raw.vcf.gz"
FILTERED_VCF="${VCF_DIR}/allsamples.filtered.vcf.gz"

# 過濾低質量 SNPs
bcftools filter -i 'FS<40.0 && SOR<3 && MQ>40.0 && MQRankSum>-5.0 && MQRankSum<5 && QD>2.0 && ReadPosRankSum>-4.0 && INFO/DP<16000' \
    -O z -o ${FILTERED_VCF} ${RAW_VCF}

# 移除低品質個體基因型
bcftools filter -S . -e 'FMT/DP<3 | FMT/GQ<20' -O z -o ${VCF_DIR}/allsamples.filtered2.vcf.gz ${FILTERED_VCF}

# 移除多等位 SNPs、INDELs、以及單調 SNPs
bcftools filter -e 'AC==0 || AC==AN' --SnpGap 10 ${VCF_DIR}/allsamples.filtered2.vcf.gz | \
    bcftools view -m2 -M2 -v snps -O z -o ${VCF_DIR}/allsamples.final.vcf.gz

echo "SNP filtering completed!"

