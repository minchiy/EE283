#!/bin/bash
#SBATCH -J featureCounts
#SBATCH -n 8          # 要 8 cores
#SBATCH --mem=16G     # 依照資料大小調整
#SBATCH -t 12:00:00   # 時間可視需要調整

module load subread/2.0.3

# 註: dmel-all-r6.13.gtf 跟課堂範例相近
gtf="dmel-all-r6.13.gtf"

# 讀取 shortRNAseq.names.txt 內容，將換行替換成空白，以便用在 featureCounts
myfile=`cat shortRNAseq.names.txt | tr "\n" " "`

echo "Running featureCounts on the following BAMs:"
echo $myfile

featureCounts \
  -p \
  -T 8 \
  -t exon \
  -g gene_id \
  -Q 30 \
  -F GTF \
  -a $gtf \
  -o fly_counts.txt \
  $myfile

echo "featureCounts 完成，輸出檔案: fly_counts.txt"

