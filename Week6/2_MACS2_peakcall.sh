#!/bin/bash
#SBATCH -J macs2_peak
#SBATCH -p standard      # 請改成您可用的 partition
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -o macs2_peak.out
#SBATCH -e macs2_peak.err

# 2_MACS2_peakcall.sh
# -------------------------------------------------------
# 功能:
#   1) 合併 "ED" 組織的 bam => ED_merged.bam
#   2) 轉 bed 並做 TN5 offset => ED_merged.tn5.bed
#   3) macs2 callpeak => 產生 broad peaks
#   4) 產生 bigWig (pileup)
# -------------------------------------------------------

module load samtools
module load bedtools2
module load macs2
module load ucsc-tools

# (A) 準備檔案
OUTDIR="macs2_ED"
mkdir -p "$OUTDIR"

# 假設您要合併 5 個 ED 樣本
# (根據檔案: A4_ED_2, A4_ED_3, A4_ED_4, A5_ED_1, A7_ED_2)
WORKDIR="week6_prep_output"
bams=(
  "A4_ED_2.dedup.chrX.bam"
  "A4_ED_3.dedup.chrX.bam"
  "A4_ED_4.dedup.chrX.bam"
  "A5_ED_1.dedup.chrX.bam"
  "A7_ED_2.dedup.chrX.bam"
)

# (B) 合併 BAM
mergedBAM="$OUTDIR/ED_merged.chrX.bam"
if [ ! -f "$mergedBAM" ]; then
    echo "[Merging BAM for ED samples]"
    samtools merge -o "$mergedBAM" $(printf "$WORKDIR/%s " "${bams[@]}")
    samtools index "$mergedBAM"
fi

# (C) 轉為 BED + Tn5 offset
tn5BED="$OUTDIR/ED_merged.tn5.bed"
if [ ! -f "$tn5BED" ]; then
    echo "[Converting to BED and do Tn5 offset]"
    bedtools bamtobed -i "$mergedBAM" \
    | awk 'BEGIN{OFS="\t"} {
        if($6=="+") $2=$2+4;
        else if($6=="-") $3=$3-5;
        print $0;
      }' > "$tn5BED"
fi

# (D) 呼叫peak => broad peak模式
#     注意 -g dm => fruitfly genome size? (或改 mm / hs 都只是計算P值時用)
#     --nomodel --shift -75 --extsize 150 => ATAC常用參數
cd "$OUTDIR"

macs2 callpeak \
  -t $(basename "$tn5BED") \
  -n ED \
  -f BED -g dm \
  -q 0.01 \
  --broad \
  --nomodel --shift -75 --extsize 150 --keep-dup all -B --call-summits

# (E) 轉為 bigWig
#     Macs2 會輸出 .broad_treat_pileup.bdg 檔 (或者 .bdg)
#     先 sort，再轉 bigWig
sort -k1,1 -k2,2n ED_broad_treat_pileup.bdg > ED_broad_treat_pileup.sorted.bdg

CHROMSIZE="/pub/minchiy/EE283/Week4/dm6.chrom.sizes"  # 請改成您實際的位置
bedGraphToBigWig ED_broad_treat_pileup.sorted.bdg "$CHROMSIZE" ED_broad_peaks.bw

echo "MACS2 + BigWig done. Result in $OUTDIR"

