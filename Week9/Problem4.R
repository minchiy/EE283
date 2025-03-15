# 加載必要的包
library(tidyverse)
library(qqman)
library(gridExtra)
library(grid)

# 使用前面創建的GWAS數據

# 將兩個模型的數據整合在一起用於比較
comparison_data <- inner_join(
  gwas_data_model1 %>% select(CHR, BP, SNP, P) %>% rename(P1 = P),
  gwas_data_model2 %>% select(CHR, BP, SNP, P) %>% rename(P2 = P),
  by = c("CHR", "BP", "SNP")
)

# 創建比較散點圖
p_comparison <- ggplot(comparison_data, aes(x = -log10(P1), y = -log10(P2))) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "比較兩個模型的-log10(p)值",
       x = "Model 1: -log10(p)",
       y = "Model 2: -log10(p)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

# 保存按照指定佈局的三面板圖
tiff("three_panel_comparison.tiff", width = 10, height = 8, units = "in", res = 300)

# 創建layout matrix按照要求：上面寬的兩個曼哈頓圖，下面是方形的散點圖
layout_matrix <- rbind(
  c(1, 1, 1),
  c(2, 2, 2),
  c(3, 3, 3)
)

# 創建臨時GG對象
manhattan_plot1 <- function() {
  manhattan(gwas_data_model1, main = "Malathion GWAS - Model 1",
            suggestiveline = -log10(1e-4),
            genomewideline = -log10(5e-6),
            ylim = c(0, 8))
}

manhattan_plot2 <- function() {
  manhattan(gwas_data_model2, main = "Malathion GWAS - Model 2",
            suggestiveline = -log10(1e-4),
            genomewideline = -log10(5e-6),
            ylim = c(0, 8))
}

# 使用grid.arrange和arrangeGrob佈局
grid.arrange(
  ggplotGrob(ggplot() + annotation_custom(grob = grid::textGrob("")) + theme_void()),
  ggplotGrob(ggplot() + annotation_custom(grob = grid::textGrob("")) + theme_void()),
  p_comparison,
  layout_matrix = layout_matrix
)

# 在第一個和第二個面板繪製曼哈頓圖
vp1 <- viewport(layout.pos.row = 1, layout.pos.col = 1:3)
pushViewport(vp1)
manhattan_plot1()
popViewport()

vp2 <- viewport(layout.pos.row = 2, layout.pos.col = 1:3)
pushViewport(vp2)
manhattan_plot2()
popViewport()

dev.off()
