#!/bin/bash
#SBATCH --job-name=problem1
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=6:00:00

set -euo pipefail

cd /pub/minchiy/EE283/Week4
exec > >(tee problem1.log) 2>&1

# 設定路徑
DATADIR="/pub/minchiy/ee283_work/Raw_data/DNA_Seq"
WORKDIR="/pub/minchiy/EE283/Week4"
REF="$WORKDIR/dmel-all-chromosome-r6.13.fasta"
TMPDIR="$WORKDIR/tmp"

# 清理並創建臨時目錄
rm -rf "$TMPDIR"
mkdir -p "$TMPDIR"

# 載入模組
module load bwa/0.7.17
module load samtools/1.15.1

# 定義樣本對應關係
declare -A samples=(
    ["A4"]="ADL06"
    ["A5"]="ADL09"
)

# 處理每個樣本
for strain in "${!samples[@]}"; do
    sample="${samples[$strain]}"
    echo "Processing ${strain} (${sample})..."
    
    for rep in {1..3}; do
        echo "Processing replicate ${rep}..."
        
        R1="${DATADIR}/${sample}_${rep}_1.fq.gz"
        R2="${DATADIR}/${sample}_${rep}_2.fq.gz"
        temp_bam="${TMPDIR}/${sample}_${rep}.temp.bam"
        sorted_bam="${TMPDIR}/${sample}_${rep}.sorted.bam"
        region_bam="${TMPDIR}/${strain}_${rep}_chrX_region.bam"
        
        # BWA 比對
        echo "Aligning ${sample} replicate ${rep}..."
        bwa mem -t 16 "$REF" "$R1" "$R2" | \
            samtools view -b -h -o "$temp_bam" -

        # 排序並索引BAM
        echo "Sorting and indexing BAM for ${sample} replicate ${rep}..."
        samtools sort -@ 16 -o "$sorted_bam" "$temp_bam"
        samtools index "$sorted_bam"
        
        # 提取區域
        echo "Extracting region for ${sample} replicate ${rep}..."
        samtools view -b -h -q 30 -o "$region_bam" "$sorted_bam" chrX:1880000-2000000
        samtools index "$region_bam"
        
        # 清理中間文件
        rm -f "$temp_bam"
    done
    
    # 合併重複
    echo "Merging replicates for ${strain}..."
    samtools merge -f "${WORKDIR}/${strain}_chrX_region.bam" \
        "${TMPDIR}/${strain}"_*_chrX_region.bam
    samtools index "${WORKDIR}/${strain}_chrX_region.bam"
done

# 清理臨時文件
rm -rf "$TMPDIR"

echo "Pipeline completed successfully"
echo "Final outputs:"
ls -lh ${WORKDIR}/*_chrX_region.bam*
