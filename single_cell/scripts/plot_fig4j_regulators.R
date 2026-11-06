# Run from single_cell/. Usage: Rscript scripts/plot_fig4j_regulators.R output.pdf
# Source: Junb high regulators_bubble plot.ipynb, cells 5–6.
# Input is the author-provided Metascape/TRRUST summary, not a new enrichment run.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
output_pdf <- args[1]
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
# ---- ライブラリ ----
library(tidyverse)
library(scales)

# ---- データ読み込み ----
df <- read_csv("data/regulator_table.csv", show_col_types = FALSE)

# Trp53を一番上にして反転（Percent降順）
df <- df %>%
  mutate(Regulated_by = forcats::fct_reorder(Regulated_by, Percent, .desc = TRUE)) %>%
  mutate(Regulated_by = forcats::fct_rev(Regulated_by))

# ---- 図作成（幅を狭く・文字大きく）----
p <- ggplot(df, aes(x = Percent, y = Regulated_by)) +
  geom_segment(aes(x = 0, xend = Percent, y = Regulated_by, yend = Regulated_by),
               linewidth = 0.8, color = "grey80") +
  geom_point(aes(size = Count, color = Count), alpha = 0.9) +
  scale_color_gradient(low = "#FEE08B", high = "#D73027", name = "Count") +
  scale_size_continuous(name = "Count", range = c(4, 14),
                        breaks = sort(unique(df$Count))) +
  scale_x_continuous(labels = label_percent(accuracy = 0.1, scale = 1),
                     expand = expansion(mult = c(0.01, 0.08))) +
  labs(x = "Percent (%)", y = "Regulated by",
       title = "Regulators: Percent vs Count") +
  theme_minimal(base_size = 16) +  # ← フォントサイズを大きく
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 13),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    legend.position = "right"
  )

# 表示
p
# ---- 保存（横幅狭く・縦長で）----
width <- 5.2
height <- max(6, nrow(df) * 0.55)
ggsave(output_pdf, p, width = width, height = height)

