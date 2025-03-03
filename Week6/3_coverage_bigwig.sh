#!/bin/bash
#SBATCH -J coverage_bw
#SBATCH -p standard
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -t 6:00:00
#SBATCH -o coverage_bw.out
#SBATCH -e coverage_bw.err

# 3_coverage_bigwig.sh
# --------------------------------------------------
# 功能:
#   1) 計算 coverage (RPM normalized) for single .bam
#   2) 輸出 BigWig
# --------------------------------------------------

module load samtools
module load bedtools2
module load ucsc-tools

CHROMSIZE="/pub/minchiy/EE283/Week4/dm6.chrom.sizes"  # 請改成實際路徑
BAM="week6_prep_output/A7_WD_1.dedup.chrX.bam"         # 範例: 單一檔案

# (A) 計算總 reads 數 (Q30, dedup 後)
Nreads=$(samtools view -c "$BAM")
echo "Total reads in $BAM => $Nreads"

# (B) RPM scale factor = 1 / (Nreads/1e6)
Scale=$(echo "scale=6; 1000000 / $Nreads" | bc -l)
echo "Scale factor = $Scale"

# (C) 產生 coverage bedGraph
COV="A7_WD_1.coverage.bedGraph"
genomeCoverageBed -ibam "$BAM" \
  -g "$CHROMSIZE" \
  -bg \
  -scale "$Scale" \
  > "$COV"

# (D) 排序 + bigWig
sort -k1,1 -k2,2n "$COV" > "${COV%.bedGraph}.sorted.bedGraph"
bedGraphToBigWig "${COV%.bedGraph}.sorted.bedGraph" "$CHROMSIZE" "A7_WD_1.coverage.bw"

echo "Done. coverage BigWig => A7_WD_1.coverage.bw"


