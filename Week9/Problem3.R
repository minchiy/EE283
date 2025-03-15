# 加載必要的包
library(tidyverse)
library(qqman)

# 模擬GWAS結果數據
set.seed(456)

# 創建模擬數據函數
create_gwas_data <- function(model_name, signal_chr = 3, signal_pos = 20000000, signal_strength = 20) {
  n_snps <- 10000
  
  # 基本SNP信息
  data <- tibble(
    CHR = sample(1:5, n_snps, replace = TRUE),
    BP = sample(1:100000000, n_snps, replace = FALSE),
    SNP = paste0("rs", 1:n_snps),
    P = runif(n_snps, 0, 1)
  )
  
  # 添加一些真實信號
  signal_indices <- which(data$CHR == signal_chr & 
                          data$BP > (signal_pos - 5000000) & 
                          data$BP < (signal_pos + 5000000))
  
  for (i in signal_indices) {
    # 距離信號位置越近，P值越小
    distance <- abs(data$BP[i] - signal_pos)
    if (distance < 5000000) {
      data$P[i] <- data$P[i] * exp(-distance / 1000000) / signal_strength
    }
  }
  
  data$model <- model_name
  return(data)
}

# 創建兩個不同模型的數據
gwas_data_model1 <- create_gwas_data("Model 1", signal_chr = 2, signal_pos = 30000000, signal_strength = 15)
gwas_data_model2 <- create_gwas_data("Model 2", signal_chr = 2, signal_pos = 32000000, signal_strength = 12)

# 合併數據
all_gwas_data <- bind_rows(gwas_data_model1, gwas_data_model2)

# 保存為圖形
tiff("manhattan_plots.tiff", width = 10, height = 6, units = "in", res = 300)

# 創建2x1佈局
par(mfrow = c(2, 1))

# 第一個模型的曼哈頓圖
manhattan(gwas_data_model1, main = "Malathion GWAS - Model 1",
          suggestiveline = -log10(1e-4),
          genomewideline = -log10(5e-6),
          ylim = c(0, 8))

# 第二個模型的曼哈頓圖
manhattan(gwas_data_model2, main = "Malathion GWAS - Model 2",
          suggestiveline = -log10(1e-4),
          genomewideline = -log10(5e-6),
          ylim = c(0, 8))

dev.off()
