#!/bin/bash
#SBATCH --job-name=SNP_calling
#SBATCH --output=snp_calling.out
#SBATCH --error=snp_calling.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=36

# 加載模組
module load gatk/4.2.6.1
module load samtools/1.15.1

# 設定路徑
BAM_DIR="/pub/minchiy/ee283_work/Processed_data/BAM"
OUT_DIR="/pub/minchiy/EE283/Week5/SNPs"
REF_GENOME="/pub/minchiy/EE283/Week4/dmel-all-chromosome-r6.13.fasta"

mkdir -p ${OUT_DIR}

SAMPLES=("ADL06" "ADL09")

# SNP 呼叫
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Calling SNPs for sample: ${SAMPLE}"

    gatk --java-options "-Xmx8g" HaplotypeCaller \
        -R ${REF_GENOME} \
        -I ${BAM_DIR}/${SAMPLE}.dedup.bam \
        --minimum-mapping-quality 30 \
        -ERC GVCF \
        -O ${OUT_DIR}/${SAMPLE}.g.vcf.gz
done

# 合併所有樣本的 GVCF
echo "Merging all sample GVCFs..."
gatk CombineGVCFs -R ${REF_GENOME} \
    $(printf -- "-V ${OUT_DIR}/%s.g.vcf.gz " ${SAMPLES[@]}) \
    -O ${OUT_DIR}/allsamples.g.vcf.gz

# 呼叫變異
echo "Genotyping variants..."
gatk --java-options "-Xmx8g" GenotypeGVCFs \
    -R ${REF_GENOME} \
    -V ${OUT_DIR}/allsamples.g.vcf.gz \
    -O ${OUT_DIR}/allsamples.raw.vcf.gz

echo "SNP calling completed!"

