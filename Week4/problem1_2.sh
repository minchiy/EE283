#!/bin/bash
#SBATCH --job-name=dna_seq_problem1
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

# 創建臨時目錄
mkdir -p "$TMPDIR"

# 清理舊的臨時文件
rm -f $WORKDIR/*.tmp.* $WORKDIR/*.temp.*

# 載入所需模組
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
    
    # 檢查是否已經有最終的 BAM 文件
    if [ -f "${WORKDIR}/${strain}_chrX_region.bam" ]; then
        echo "Final BAM for ${strain} already exists, skipping..."
        continue
    fi
    
    for rep in {1..3}; do
        # 檢查是否已經有區域 BAM 文件
        if [ -f "${WORKDIR}/${sample}_${rep}_region.bam" ]; then
            echo "Region BAM for ${sample} replicate ${rep} already exists, skipping..."
            continue
        fi
        
        echo "Processing replicate ${rep}..."
        R1="${DATADIR}/${sample}_${rep}_1.fq.gz"
        R2="${DATADIR}/${sample}_${rep}_2.fq.gz"
        temp_bam="${TMPDIR}/${sample}_${rep}.temp.bam"
        sorted_bam="${WORKDIR}/${sample}_${rep}.sorted.bam"
        
        # BWA 比對
        echo "Aligning replicate ${rep}..."
        bwa mem -t 16 "$REF" "$R1" "$R2" | \
            samtools view -b -h -o "$temp_bam" -
        
        # 排序 BAM
        echo "Sorting BAM file..."
        samtools sort -T "$TMPDIR/${sample}_${rep}" \
                     -o "$sorted_bam" "$temp_bam"
        samtools index "$sorted_bam"
        
        # 提取區域
        echo "Extracting region..."
        samtools view -b -h -q 30 -o "${WORKDIR}/${sample}_${rep}_region.bam" \
            "$sorted_bam" chrX:1880000-2000000
        
        # 清理文件
        rm -f "$temp_bam" "$sorted_bam" "$sorted_bam.bai"
    done
    
    # 合併重複
    echo "Merging replicates for ${strain}..."
    samtools merge -f "${WORKDIR}/${strain}.bam" \
        "${WORKDIR}/${sample}"_*_region.bam
    
    # 提取最終區域
    echo "Creating final BAM for ${strain}..."
    samtools view -b -h -q 30 -o "${WORKDIR}/${strain}_chrX_region.bam" \
        "${WORKDIR}/${strain}.bam" chrX:1880000-2000000
    samtools index "${WORKDIR}/${strain}_chrX_region.bam"
    
    # 清理中間文件
    rm -f "${WORKDIR}/${sample}"_*_region.bam "${WORKDIR}/${strain}.bam"
done

# 清理臨時目錄
rm -rf "$TMPDIR"

echo "Pipeline completed successfully"
echo "Final outputs:"
ls -lh ${WORKDIR}/*_chrX_region.bam*
