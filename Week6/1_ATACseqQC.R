#!/usr/bin/env Rscript

# 1_ATACseqQC_allTSS.R
# -----------------------------------------------------------------
# 功能:
#   1) 對多個 ATAC-seq BAM (已去重 + 只保留 chrX) 做QC
#   2) 計算 fragment size distribution、library complexity
#   3) 對「所有樣本」各別進行 TSS enrichment (Tn5 shift, split nucleosomes, heatmap)
# -----------------------------------------------------------------
library(ATACseqQC)
library(Rsamtools)
library(GenomicFeatures)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene) # 若您使用dm6
library(ChIPpeakAnno)
library(BiocParallel)
library(BiocIO) # for export(gal1, ...)

# 平行運算 (可改成自己想用的核心數)
register(MulticoreParam(4))  

# (A) 設定檔案路徑: 假設您在 'week6_prep_output/' 有下列樣本
bamDir <- "week6_prep_output"
samples <- c("A4_ED_2", "A4_ED_3", "A4_ED_4",
             "A5_ED_1", "A7_ED_2", "A7_WD_1") 
bams <- file.path(bamDir, paste0(samples, ".dedup.chrX.bam"))
stopifnot(all(file.exists(bams)))  # 檔案存在檢查

# (B) 片段大小分布 & 圖形輸出
pdf("week6_fragSizeDist.pdf", width=7, height=6)
for(i in seq_along(bams)){
  cat("Plotting fragment size for:", bams[i], "\n")
  label <- samples[i]
  fragSizeDist(bams[i], bamfile.labels = label)
}
dev.off()
cat("已輸出片段大小分布: week6_fragSizeDist.pdf\n")

# (C) Library Complexity (dupFreq => estimateLibComplexity)
cat("\n==== Library Complexity ====\n")
for(i in seq_along(bams)){
  cat("\nSample:", samples[i], "\n")
  dfreq <- readsDupFreq(bams[i])
  est <- estimateLibComplexity(dfreq)
  print(est)
}
cat("\n===========================\n")

# (D) TSS enrichment for **each sample** 
# 建立 TxDb => 取得 TSS
txs <- transcripts(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)

# 建立一個主資料夾放所有 sample 的 TSS 結果
mainOutPath <- "week6_TSS_QC_all"
dir.create(mainOutPath, showWarnings=FALSE)

for(i in seq_along(bams)){
  sampleName <- samples[i]
  bamfile <- bams[i]
  cat("\n=== TSS Enrichment on:", sampleName, "===\n")

  # 1) 建立輸出資料夾
  outPath <- file.path(mainOutPath, sampleName)
  dir.create(outPath, showWarnings=FALSE, recursive=TRUE)

  # 2) 讀入 BAM, shift by Tn5 offsets
  gal <- readBamFile(bamfile, asMates=TRUE, bigFile=TRUE)
  gal1 <- shiftGAlignmentsList(gal)
  
  # 輸出 shift 後的 BAM
  shiftedBamfile <- file.path(outPath, paste0(sampleName, "_shifted.bam"))
  export(gal1, shiftedBamfile)
  
  # 3) split nucleosome-free, mono, di, ...
  objs <- splitGAlignmentsByCut(gal1, txs=txs, outPath=outPath)
  null <- writeListOfGAlignments(objs, outPath=outPath)

  # 4) enrichedFragments => TSS heatmap
  bamfiles <- file.path(outPath, c("NucleosomeFree.bam", 
                                   "mononucleosome.bam", 
                                   "dinucleosome.bam", 
                                   "trinucleosome.bam"))
  librarySize <- estLibSize(bamfiles)
  NTILE <- 101
  ups   <- 1010
  dws   <- 1010
  
  sigs <- enrichedFragments(bamfiles, TSS=TSS, librarySize=librarySize,
                            seqlev="chrX",  # 或 c("chrX","chr2L","chr2R","chr3L","chr3R")
                            TSS.filter=0.5,
                            n.tile=NTILE,
                            upstream=ups, 
                            downstream=dws)

  sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
  
  pdf(file.path(outPath, paste0(sampleName, "_TSSenrichment.pdf")), width=8, height=6)
  featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                        zeroAt=0.5, n.tile=NTILE)
  dev.off()

  cat("Done sample:", sampleName, " =>", outPath, "\n")
}

cat("\nAll samples TSS enrichment done! See:", mainOutPath, "\n")

