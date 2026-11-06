# Source excerpt: Endothelial cell analysis/scGSVA_EachOrgan/Adipose/scGSVA_Adipose_Reactome_correlation_with regulons-Spearman.ipynb
# Exact network visualization block. Requires All with AUC assay and Reactome metadata and Reactome_path vector.
# Cell numbers are one-based notebook positions, not execution counts.

# ---- Original cell 58 ----
# 必要なライブラリを読み込む
library(Seurat)
library(igraph)

# ターゲットのリスト
targets <- c("Irf7", "Irf9", "Stat1", "Stat2", "Junb", "Pparg", "Cebpd", "Jund", "Jun", "Klf2")

# 結果を格納するリスト
result_list <- list()
top20_list <- list()

# 各ターゲットについてループ
for (target in targets) {
  # regulonのデータをRNAアッセイから取得
  data <- GetAssayData(All, assay = "regulonAUC_extended", slot = "data")[target, ]
  
  # 相関係数を計算し、結果をリストに格納
  correlations <- list()
  
  # パスウェイ名のリストReactome_pathsを使用
  for (path in Reactome_path) {
    pathway_data <- All@meta.data[[path]]  # 仮のパスウェイデータを取得
    correlation <- cor(data, pathway_data, method = "spearman")  # 相関を計算
    
    # 相関が欠損値でない場合のみリストに追加
    if (!is.na(correlation)) {
      correlations[[path]] <- correlation  # 相関係数をリストに追加
    }
  }
  
  # 相関係数の値でソート
  correlation_sorted <- sort(unlist(correlations), decreasing = TRUE)
  
  # 選択されたパスウェイでデータフレームを作成
  df <- data.frame(Pathway = names(correlation_sorted), Coefficient = correlation_sorted)
  
  # トップ20の結果に絞り込む
  top20_df <- head(df, 20)
  
  # 結果をリストに格納
  result_list[[target]] <- df
  top20_list[[target]] <- top20_df
  
  # 動的に変数名を設定
  assign(paste0("df_", target), df)
  assign(paste0("top20_df_", target), top20_df)
}

# 各データフレームからPathway列のベクトルを取得
pathways_list <- lapply(top20_list, function(df) df$Pathway)
names(pathways_list) <- targets

# 全てのPathwayのリストを結合して、ユニークなPathwayを識別
all_pathways <- unique(unlist(pathways_list))

# エッジリストを初期化
edges <- data.frame(from = character(), to = character(), Coefficient = numeric(), stringsAsFactors = FALSE)

# 各データセットとPathwayの関連をエッジリストに追加
for (dataset_name in names(top20_list)) {
  for (i in 1:nrow(top20_list[[dataset_name]])) {
    edges <- rbind(edges, data.frame(from = dataset_name, to = top20_list[[dataset_name]]$Pathway[i], Coefficient = top20_list[[dataset_name]]$Coefficient[i]))
  }
}

# 重複を削除
edges <- unique(edges)

# エッジリストからグラフオブジェクトを作成
network_graph <- graph_from_data_frame(edges, directed = FALSE)

# 相関係数が0.3より大きいノードを特定
significant_nodes <- unique(edges$to[edges$Coefficient > 0.3])

# 全ノードのラベルを初期化（デフォルトでは非表示）
V(network_graph)$label <- NA

# 相関係数が0.3より大きいノードとターゲットノードだけラベルを表示
V(network_graph)$label[V(network_graph)$name %in% c(significant_nodes, targets)] <- V(network_graph)$name[V(network_graph)$name %in% c(significant_nodes, targets)]

# ノードの外周を初期化
V(network_graph)$frame.color <- "grey"
V(network_graph)$frame.width <- 1

# 相関係数が0.3より大きいノードの外周を太くする
V(network_graph)$frame.width[V(network_graph)$name %in% significant_nodes] <- 2

dataset_colors <- c("Irf7" = "royalblue1", "Irf9" = "royalblue1", "Stat1" = "royalblue1", "Stat2" = "royalblue1", "Junb" = "royalblue1",
                    "Pparg" = "coral", "Cebpd" = "coral", "Jund" = "coral", "Jun" = "coral", "Klf2" = "coral")

V(network_graph)$color <- "pink"
V(network_graph)$color[V(network_graph)$name %in% names(dataset_colors)] <- dataset_colors[V(network_graph)$name]

# 共通Pathwayのノードの色をskyblueに設定
common_pathways <- names(table(all_pathways)[table(all_pathways) > 1])
V(network_graph)$color[V(network_graph)$name %in% common_pathways] <- "skyblue"

# Pathwayノードのサイズを小さく設定
V(network_graph)$size <- 10  # 一般的なノードのサイズ
V(network_graph)$size[V(network_graph)$name %in% all_pathways] <- 2  # Pathwayノードのサイズを小さく設定

# top20_pathwaysのエッジの太さを調整
E(network_graph)$width <- 1
E(network_graph)$width[E(network_graph)$Coefficient > 0.3] <- 3

# エッジの色を調整
num_edges <- length(E(network_graph)$Coefficient)
E(network_graph)$color <- gray.colors(num_edges, start = 0.9, end = 0.1)[rank(E(network_graph)$Coefficient)]

# ネットワークをプロット
plot(network_graph, layout = layout_with_fr(network_graph),
     vertex.label = V(network_graph)$label,  # ラベル表示の制御
     vertex.color = V(network_graph)$color,  # ノードの色
     vertex.size = V(network_graph)$size,  # ノードのサイズ
     edge.width = E(network_graph)$width,  # エッジの太さ
     edge.color = E(network_graph)$color,  # エッジの色
     edge.arrow.size = .5,
     vertex.frame.color = V(network_graph)$frame.color,  # 外周の色
     vertex.frame.width = V(network_graph)$frame.width)  # 外周の太さ

