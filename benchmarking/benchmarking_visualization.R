metric <- read.csv("/Users/zhangxinyao/0515/metrics_result1.csv")
metric <- metric[,-ncol(metric)]
library(ggplot2)
library(tidyr)
df_long <- metric %>%
  pivot_longer(cols = -method, names_to = "Metric", values_to = "Score")

ggplot(df_long, aes(x = method, y = Score, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Benchmarking of scMVP, Schema, and Seurat",
       x = "Method", y = "Score") +
  scale_fill_brewer(palette = "Set2") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5))