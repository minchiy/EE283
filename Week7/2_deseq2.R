#!/usr/bin/env Rscript

# ---------------------------------------------------------------------
# 2_deseq2.R
#
# 功能：
#   1. 讀入 fly_counts.txt (featureCounts 結果)
#   2. 讀入 shortRNAseq.txt (篩選後的樣本資訊)
#   3. 建立 DESeq2 物件並運行 DESeq 分析
#   4. 輸出結果表與部分檢視圖
# ---------------------------------------------------------------------

# 若還沒安裝 DESeq2，請先透過 BiocManager 安裝
# if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")

library(DESeq2)

# 讀入樣本設計表
sampleInfo <- read.table("shortRNAseq.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# 確保有 FullSampleName 這個欄位
sampleInfo$FullSampleName <- as.character(sampleInfo$FullSampleName)

# 讀入 featureCounts 產生的表格
countdata <- read.table("fly_counts.txt", header=TRUE, row.names=1, check.names=FALSE)
# 去除前 5 欄 (Chr, Start, End, Strand, Length)
countdata <- countdata[ , 6:ncol(countdata)]

# 整理 colnames 使其和 sampleInfo$FullSampleName 一一對應
temp <- colnames(countdata)
# 假設 featureCounts 出來的欄位是 "xxx.bam"
temp <- gsub(".bam", "", temp)
colnames(countdata) <- temp

# 確認順序是否與 sampleInfo 對應
# 如果不一樣，後面 DESeqDataSetFromMatrix() 會依名稱自行對應，但習慣上可以查看
cat("檢查對應是否一一符合:\n")
print(cbind(FeatureCountsCol=temp, SampleTableCol=sampleInfo$FullSampleName))

# 建立 DESeq2 物件，設計式這邊以 TissueCode 為示範
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData   = sampleInfo,
                              design    = ~ TissueCode)

# 執行 DESeq
dds <- DESeq(dds)

# 取得差異表達結果 (預設比較：若 TissueCode 有兩個以上水準，會以字母排序)
res <- results(dds)
# 看一下結果前 6 筆
cat("DESeq2 結果前 6 筆:\n")
print(head(res))

# 輸出整張結果表
write.csv(as.data.frame(res), file="deseq2_results.csv", row.names=TRUE)
cat("已輸出差異表達結果: deseq2_results.csv\n")

# 簡單畫幾張圖:
pdf("deseq2_QCplots.pdf")

# 1) MA plot
plotMA(res, main="DESeq2 MA-plot", ylim=c(-2,2))

# 2) p-value histogram
hist(res$pvalue, breaks=20, col="grey", main="p-value histogram")

# 3) dispersion plot
plotDispEsts(dds, main="Dispersion Estimates")

# 4) PCA plot
#   先做 rlog or vst 轉換 (此例用 rlog)
rld <- rlog(dds)
plotPCA(rld, intgroup="TissueCode")

dev.off()
cat("已繪製 QC 圖到檔案: deseq2_QCplots.pdf\n")

