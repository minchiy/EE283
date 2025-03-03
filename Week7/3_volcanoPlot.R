#!/usr/bin/env Rscript

# ---------------------------------------------------------------------
# 3_volcanoPlot.R
#
# 功能：
#   1. 讀入上一階段的 DESeq2 結果表 (deseq2_results.csv)
#   2. 繪製火山圖 (volcano plot)
#   3. 可視情況加上閾值線 (threshold)
# ---------------------------------------------------------------------

# 請依需要自行調整繪圖參數

# 讀入差異表達結果
res <- read.csv("deseq2_results.csv", header=TRUE, row.names=1)

# 若要以 padj 來畫 -log10(padj)，要先確保不為 NA 且不為 0
# 做個簡單的處理：把 NA 或 0 變成某個極小值
res$padj[is.na(res$padj)] <- 1
res$padj[res$padj == 0]   <- 1e-300

# 設定 x, y 軸
log2FC <- res$log2FoldChange
negLog10Padj <- -log10(res$padj)

pdf("VolcanoPlot.pdf")

plot(log2FC,
     negLog10Padj,
     pch=20,
     main="Volcano Plot",
     xlab="log2 Fold Change",
     ylab="-log10 (padj)")

# 顯著閾值(舉例) : padj < 0.05 與 |log2FC| > 1
alpha_cutoff <- 0.05
lfc_cutoff   <- 1

abline(h = -log10(alpha_cutoff), col="red", lty=2)   # padj=0.05
abline(v =  c(-lfc_cutoff, lfc_cutoff), col="blue", lty=2)

dev.off()

cat("已輸出火山圖: VolcanoPlot.pdf\n")

