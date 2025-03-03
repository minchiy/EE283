#!/bin/bash
#SBATCH --job-name=demoRNAseq_multi   # 任務名稱
#SBATCH --output=demoRNAseq_multi.out # 標準輸出檔案
#SBATCH --error=demoRNAseq_multi.err  # 錯誤訊息檔案
#SBATCH --time=10:00:00               # 最長執行時間 4 小時 (依需要調整)
#SBATCH --mem=32G                     # 記憶體需求 (依需要調整)
#SBATCH --cpus-per-task=32            # 使用 CPU core 數量 (依需要調整)
#SBATCH --nodes=1

# ------------------------------------------------------------------------------
# (1) 載入所需模組
module load hisat2
module load samtools

# ------------------------------------------------------------------------------
# (2) 定義參考基因組索引
index="/pub/minchiy/EE283/Week4/dmel_r6.13_index"

# ------------------------------------------------------------------------------
# (3) 定義樣本列表
samples=("21148E0" "21286E0" "22162E0" "21297E0" "21029E0" "22052E0" "22031E0" "21293E0" "22378E0" "22390E0")

# ------------------------------------------------------------------------------
# (4) 設定 FASTQ 目錄
fastq_dir="/pub/minchiy/ee283_work/Raw_data/RNA_Seq"

# ------------------------------------------------------------------------------
# (5) 函式: 針對每個樣本執行 HISAT2 比對並處理 BAM
run_hisat2_and_sort() {
  local sname="$1"
  local R1="${fastq_dir}/*_${sname}_R1.fq.gz"
  local R2="${fastq_dir}/*_${sname}_R2.fq.gz"

  echo "=== [$(date)] Start HISAT2 mapping for sample: $sname ==="

  hisat2 \
    -x $index \
    -1 $R1 \
    -2 $R2 \
    --threads $SLURM_CPUS_PER_TASK \
    --rg-id $sname \
    --rg "SM:${sname}" \
  | samtools view -@ $SLURM_CPUS_PER_TASK -bS - \
  > ${sname}.bam

  # 排序
  samtools sort -@ $SLURM_CPUS_PER_TASK -o ${sname}.sorted.bam ${sname}.bam

  # 建立索引
  samtools index ${sname}.sorted.bam

  echo "=== [$(date)] Done processing: $sname ==="
}

# ------------------------------------------------------------------------------
# (6) 針對所有樣本執行處理
for sample in "${samples[@]}"; do
  run_hisat2_and_sort "$sample"
done

echo "=== All Samples Processed Successfully! ==="

