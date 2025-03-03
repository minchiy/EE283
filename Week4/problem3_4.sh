#!/bin/bash
#SBATCH --job-name=problem3_4
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=6:00:00

set -euo pipefail

cd /pub/minchiy/EE283/Week4
exec > >(tee problem3_4.log) 2>&1

# 設置 conda 
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh

# 創建並設置 deeptools 環境
if ! conda env list | grep -q "deeptools_env"; then
    echo "Creating deeptools environment..."
    mamba create -y -n deeptools_env deeptools
fi

conda activate deeptools_env

echo "Starting deeptools analysis..."

# 使用第一題生成的 BAM 文件
A4="A4_chrX_region.bam"
A5="A5_chrX_region.bam"

# 第3題：計算每個鹼基的讀序覆蓋度（使用RPKM標準化）
echo "Computing per base coverage..."
bamCoverage -b "$A4" \
    -o A4_coverage.bw \
    --normalizeUsing RPKM \
    -p 16

bamCoverage -b "$A5" \
    -o A5_coverage.bw \
    --normalizeUsing RPKM \
    -p 16

# 繪製覆蓋度比較圖
plotCoverage -b "$A4" "$A5" \
    -o coverage_comparison.png \
    --plotTitle "Coverage Comparison A4 vs A5" \
    -p 16

# 第4題：使用片段覆蓋度（--extendReads）
echo "Computing fragment coverage..."
bamCoverage -b "$A4" \
    -o A4_fragment_coverage.bw \
    --normalizeUsing RPKM \
    --extendReads \
    -p 16

bamCoverage -b "$A5" \
    -o A5_fragment_coverage.bw \
    --normalizeUsing RPKM \
    --extendReads \
    -p 16

# 繪製片段覆蓋度比較圖
plotCoverage -b "$A4" "$A5" \
    -o fragment_coverage_comparison.png \
    --plotTitle "Fragment Coverage Comparison A4 vs A5" \
    --extendReads \
    -p 16

# 特別關注 1,904,042 位置附近的區域
echo "Creating detailed view around position 1,904,042..."
plotProfile -b A4_coverage.bw A5_coverage.bw \
    -o coverage_profile_1904042.png \
    --plotTitle "Coverage Profile around position 1,904,042" \
    -r chrX:1903542-1904542 \
    --yMin 0

plotProfile -b A4_fragment_coverage.bw A5_fragment_coverage.bw \
    -o fragment_coverage_profile_1904042.png \
    --plotTitle "Fragment Coverage Profile around position 1,904,042" \
    -r chrX:1903542-1904542 \
    --yMin 0

# 生成覆蓋度比較圖
computeMatrix scale-regions -S A4_coverage.bw A5_coverage.bw \
    -R chrX:1880000-2000000 \
    -o matrix.gz \
    --regionBodyLength 1000 \
    -p 16

plotHeatmap -m matrix.gz \
    -o coverage_heatmap.png \
    --colorMap RdYlBu \
    --whatToShow 'heatmap and colorbar' \
    --zMin 0

conda deactivate

echo "Analysis complete. Output files:"
ls -lh *.bw *.png

# 第5題的結果可以在 Santa Cruz Genome Browser 查看：
echo "For comparison with reference data, visit:"
echo "http://goo.gl/LLpoNH"
