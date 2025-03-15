# 加載必要的包
library(tidyverse)
library(gridExtra)
library(grid)

# 模擬ATAC-seq數據
set.seed(123)

# 函數生成模擬的片段長度數據
generate_fragment_data <- function(sample_name, peak_positions, peak_heights) {
  tibble(
    sample = sample_name,
    fragment_length = rep(1:1000, each = 10),
    count = sapply(1:1000, function(x) {
      sum(sapply(1:length(peak_positions), function(i) {
        peak_heights[i] * exp(-0.5 * ((x - peak_positions[i]) / 50)^2)
      })) + rnorm(1, 0, 5)
    })
  ) %>% 
    filter(count > 0, fragment_length <= 800)
}

# 函數生成模擬的TSS富集數據
generate_tss_data <- function(sample_name, peak_height, shift = 0) {
  tibble(
    sample = sample_name,
    distance_from_tss = seq(-1000, 1000, by = 10),
    enrichment = peak_height * exp(-0.5 * ((distance_from_tss - shift) / 200)^2) + rnorm(201, 0, 0.05)
  )
}

# 生成4個樣本的數據
fragment_data_1 <- generate_fragment_data("Sample1", c(100, 200, 300, 400, 500), c(100, 80, 60, 40, 20))
fragment_data_2 <- generate_fragment_data("Sample2", c(100, 200, 300, 400, 500), c(90, 70, 50, 30, 15))
fragment_data_3 <- generate_fragment_data("Sample3", c(100, 200, 300, 400, 500), c(80, 60, 40, 20, 10))
fragment_data_4 <- generate_fragment_data("Sample4", c(100, 200, 300, 400, 500), c(70, 50, 30, 10, 5))

tss_data_1 <- generate_tss_data("Sample1", 2.5, -20)
tss_data_2 <- generate_tss_data("Sample2", 2.3, 0)
tss_data_3 <- generate_tss_data("Sample3", 2.1, 20)
tss_data_4 <- generate_tss_data("Sample4", 1.9, 40)

# 合併數據
all_fragment_data <- bind_rows(fragment_data_1, fragment_data_2, fragment_data_3, fragment_data_4)
all_tss_data <- bind_rows(tss_data_1, tss_data_2, tss_data_3, tss_data_4)

# 創建片段長度圖
fragment_plots <- list()
samples <- unique(all_fragment_data$sample)

for (i in 1:length(samples)) {
  sample_data <- all_fragment_data %>% filter(sample == samples[i])
  
  fragment_plots[[i]] <- ggplot(sample_data, aes(x = fragment_length, y = count)) +
    geom_line(color = "blue") +
    labs(title = paste(samples[i], "- Fragment Length Distribution"),
         x = "Fragment Length (bp)",
         y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 7))
}

# 創建TSS富集圖
tss_plots <- list()
for (i in 1:length(samples)) {
  sample_data <- all_tss_data %>% filter(sample == samples[i])
  
  tss_plots[[i]] <- ggplot(sample_data, aes(x = distance_from_tss, y = enrichment)) +
    geom_line(color = "red") +
    labs(title = paste(samples[i], "- TSS Enrichment"),
         x = "Distance from TSS (bp)",
         y = "Enrichment") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 7))
}

# 定義佈局 - 左側4個片段長度圖，右側4個TSS富集圖
tiff("atacseq_analysis.tiff", width = 10, height = 8, units = "in", res = 300)
grid.arrange(
  fragment_plots[[1]], tss_plots[[1]],
  fragment_plots[[2]], tss_plots[[2]],
  fragment_plots[[3]], tss_plots[[3]],
  fragment_plots[[4]], tss_plots[[4]],
  ncol = 2
)
dev.off()
