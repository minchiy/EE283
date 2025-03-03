#!/usr/bin/env Rscript

# ---------------------------------------------------------------------
# 1_subset.R
#
# 功能：
#   1. 讀入 "RNAseq.samcode.txt"
#   2. 使用 base R 進行資料篩選
#   3. 輸出 "shortRNAseq.txt" (包含子集樣本的實驗設計表)
#   4. 輸出 "shortRNAseq.names.txt" (包含對應 .bam 的清單)
#
# 注意：請確認您的檔名和路徑是否正確
# ---------------------------------------------------------------------

# 讀入 samcode 資料表
mytab <- read.table("RNAseq.samcode.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# 只選取這些 RILcode
ril_keep <- c(21148, 21286, 22162, 21297, 21029, 22052, 22031, 21293, 22378, 22390)

# 使用 base R 條件篩選
# 1) RILcode 在 ril_keep 裡
# 2) TissueCode 在 B 或 E
# 3) Replicate == 0
sub_mytab <- mytab[
  mytab$RILcode %in% ril_keep &
  mytab$TissueCode %in% c("B","E") &
  mytab$Replicate == 0,
  c("RILcode","TissueCode","Replicate","FullSampleName")  # 只取這四欄
]

# 檢查篩選結果
cat("篩選後的資料筆數：", nrow(sub_mytab), "\n")
cat("篩選後的前幾筆資料：\n")
print(head(sub_mytab))

# 輸出 shortRNAseq.names.txt
# 內容是一行一個 .bam 路徑 (假設 .bam 跟當前目錄下或相對路徑)
# 若您存放 bam 檔時，是在同一目錄下，直接用 "FullSampleName.bam"
# 如果像範例中是在同一層，就這樣寫：
file_out_names <- "shortRNAseq.names.txt"
if(file.exists(file_out_names)) {
  # 如果要保險，可以先移除舊檔或中斷
  file.remove(file_out_names)
}

for(i in seq_len(nrow(sub_mytab))) {
  # 每行輸出一個 .bam 檔名
  # 若您的 bam 檔實際路徑不是在同一層，要改成 "RNAseq/bam/xxx.bam" 之類
  # 課程示範常用 "RNAseq/bam/xxxx.bam"；您可依自己路徑調整
  cat(sub_mytab$FullSampleName[i], ".bam\n", sep="",
      file=file_out_names, append=TRUE)
}

cat("已寫出：", file_out_names, "\n")

# 輸出 shortRNAseq.txt
file_out_info <- "shortRNAseq.txt"
write.table(sub_mytab,
            file=file_out_info,
            sep="\t", quote=FALSE, row.names=FALSE)
cat("已寫出：", file_out_info, "\n")

