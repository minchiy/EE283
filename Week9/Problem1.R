# 加載必要的包
library(tidyverse)
library(gridExtra)
library(grid)
library(nycflights13)  # 包含flights數據集

# 創建第一個圖：距離與到達延遲的散點圖
P1 <- ggplot(flights, aes(x = distance, y = arr_delay)) +
  geom_point(alpha = 0.1) +  # 降低點的透明度以處理重疊
  geom_smooth(color = "red") +
  labs(title = "距離與到達延遲的關係",
       x = "距離 (英里)",
       y = "到達延遲 (分鐘)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

# 計算每家航空公司的平均延遲時間
temp_flights <- flights %>%
  group_by(carrier) %>%
  summarize(m_arr_delay = mean(arr_delay, na.rm = TRUE)) %>%
  arrange(desc(m_arr_delay))

# 創建第二個圖：航空公司平均延遲的條形圖
P2 <- ggplot(temp_flights, aes(x = reorder(carrier, m_arr_delay), y = m_arr_delay)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "各航空公司平均延遲",
       x = "航空公司",
       y = "平均延遲 (分鐘)") +
  coord_flip() +  # 翻轉坐標軸使標籤更易讀
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

# 創建第三個圖：航空公司延遲的箱形圖
P3 <- ggplot(flights, aes(x = reorder(carrier, arr_delay, na.rm = TRUE, FUN = median), 
                          y = arr_delay)) +
  geom_boxplot(fill = "lightgreen") +
  labs(title = "各航空公司延遲分布",
       x = "航空公司",
       y = "到達延遲 (分鐘)") +
  coord_flip() +  # 翻轉坐標軸使標籤更易讀
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

# 創建第四個圖：延遲時間的直方圖
P4 <- ggplot(flights, aes(x = arr_delay)) +
  geom_histogram(bins = 30, fill = "coral", color = "white") +
  labs(title = "到達延遲分布",
       x = "到達延遲 (分鐘)",
       y = "頻率") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

# 定義佈局
lay <- rbind(c(1, 1, 1, 2),
             c(1, 1, 1, 3),
             c(1, 1, 1, 4))

# 保存為高品質TIFF文件
tiff("figure1_improved.tiff", width = 7, height = 6, units = "in", res = 600)
grid.arrange(P1, P2, P3, P4, layout_matrix = lay)
dev.off()  # 使用dev.off()代替graphics.off()更常見
