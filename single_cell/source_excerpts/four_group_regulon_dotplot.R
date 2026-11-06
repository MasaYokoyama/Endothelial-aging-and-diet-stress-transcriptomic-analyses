# Source: SCENIC_Capillary/SCENIC_total_Adipose/Adipose regulon comparison.ipynb cell 5.
# Requires preloaded rss and zscore matrices from the paired original CSVs.
library(ggplot2)
library(scales)

# Fig.4g Adipose regulons
target_regulons <- c(
  "Irf7", "Irf9", "Stat1", "Stat2",
  "Pparg", "Jun", "Junb", "Jund", "Cebpd", "Klf2"
)

# 4群の順番
group_order <- c(
  "Normal_Young", "Normal_Aged",
  "HFD_Young", "HFD_Aged"
)

# Dot plot用データ
plot_df <- expand.grid(
  Regulon = target_regulons,
  Group = group_order,
  stringsAsFactors = FALSE
)

plot_df$RSS <- mapply(
  function(r, g) rss[r, g],
  plot_df$Regulon, plot_df$Group
)

plot_df$Z <- mapply(
  function(r, g) zscore[r, g],
  plot_df$Regulon, plot_df$Group
)

# 表示順
plot_df$Regulon <- factor(
  plot_df$Regulon,
  levels = rev(target_regulons)
)

plot_df$Group <- factor(
  plot_df$Group,
  levels = group_order,
  labels = c("3m N", "12m N", "3m HFD", "12m HFD")
)

# Dot plot
p_adipose_green <- ggplot(
  plot_df,
  aes(x = Group, y = Regulon)
) +
  geom_point(
    aes(size = RSS, color = Z)
  ) +
  scale_color_gradient(
    low = "#E5E5E5",
    high = "#1B7F1B",
    limits = c(0, 3),
    breaks = c(0, 1, 2, 3),
    oob = squish,
    name = "Z"
  ) +
  scale_size_continuous(
    breaks = c(0.1, 0.2, 0.3),
    range = c(1.5, 6),
    name = "RSS"
  ) +
  labs(
    title = "Adipose Capillary",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(
      color = "#D9D9D9",
      linewidth = 0.5
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "#4D4D4D",
      fill = NA,
      linewidth = 0.8
    ),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      color = "#4D4D4D"
    ),
    axis.text.y = element_text(color = "#4D4D4D"),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  )

print(p_adipose_green)
