# Selected original calculation cells; publication reorganizing only.
# Source: SCENIC_Capillary/SCENIC_data_analysis_Aging_Diet/SCENIC_AgingComparison_reSCENIC_CorrealtionHeatmap(TF-top10x4tissues).ipynb
# Requires the source CSVs under SOURCE_ROOT. Cell numbers are one-based.
library(dplyr)
library(pheatmap)
source_root <- Sys.getenv("SC_SOURCE_ROOT")
stopifnot(nzchar(source_root), dir.exists(source_root))
source_file <- function(path) file.path(source_root, path)
my_palette <- colorRampPalette(c("gray95", "pink", "red", "darkred", "red4"))(101)

# ---- Original cell 2 ----
Adipose<-read.csv(source_file("SCENIC_Capillary/SCENIC_Aging_Adipose/rss_Aging_Adipose_extended.csv"))
Heart<-read.csv(source_file("SCENIC_Capillary/SCENIC_Aging_Heart/rss_Aging_Heart_extended.csv"))
Muscle<-read.csv(source_file("SCENIC_Capillary/SCENIC_Aging_Muscle/rss_Aging_Muscle_extended.csv"))
Liver<-read.csv(source_file("SCENIC_Capillary/SCENIC_Aging_Liver/rss_Aging_Liver_extended.csv"))

# ---- Original cell 9 ----
# 各データフレームのTop 10遺伝子を抽出する関数
get_top_genes <- function(df) {
  df %>% 
    arrange(desc(Aged)) %>% 
    slice(1:10) %>% 
    pull(X)
}

# 各データフレームからTop 10の遺伝子名を抽出
top_TFs_adipose <- get_top_genes(Adipose)
top_TFs_heart <- get_top_genes(Heart)
top_TFs_muscle <- get_top_genes(Muscle)
top_TFs_liver <- get_top_genes(Liver)

# 一意の遺伝子名のリストを作成
unique_TFs <- unique(c(top_TFs_adipose, top_TFs_heart, top_TFs_muscle, top_TFs_liver))

# 一意の遺伝子名を表示
unique_TFs


# ---- Original cell 15 ----
List<-read.csv(source_file("Endothelial cell analysis/Vascular cells_Capillary datail/BasicFood_AgingComparison/CorrelationMatrix_Top100_Aged_NewClustering_OrganSort.csv"))

# ---- Original cell 16 ----
# ListデータフレームからX列のみを抜き出し、as.vectorを使用してベクターに変換
X_list <- as.vector(List$gene)

# 結果を表示
print(X_list)

# ---- Original cell 17 ----
# 遺伝子リストを逆順にする
genes_reversed <- rev(X_list)

# ---- Original cell 19 ----
df<-read.csv(source_file("SCENIC_Capillary/SCENIC_Aging_Adipose/output/Step2_regulonTargetsInfo.csv"), header = TRUE)


# ---- Original cell 21 ----
# TF列がtarget_stringsのいずれかを含む行を抽出
result <- df %>% filter(TF %in% unique_TFs)

# ---- Original cell 23 ----
# さらにその中でgenesに含まれる分子がターゲットに入っているものに絞る
result_gene <- result %>% filter(gene %in% genes_reversed)
result_gene

# ---- Original cell 24 ----
# データフレームをマトリックスに変換
coex_matrix <- as.matrix(result_gene)

# 行名と列名を設定
row_names <- genes_reversed
col_names <- unique(result_gene$TF)

# 空のマトリックスを作成
result_matrix <- matrix(0, nrow = length(row_names), ncol = length(col_names))

# 行名と列名を設定
rownames(result_matrix) <- row_names
colnames(result_matrix) <- col_names

# マトリックスにCoexWeightの値を入れる
for (i in 1:nrow(result_gene)) {
  row_index <- match(result_gene[i, "gene"], row_names)
  col_index <- match(result_gene[i, "TF"], col_names)
  result_matrix[row_index, col_index] <- result_gene[i, "CoexWeight"]
}

# 結果を表示
result_matrix

# ---- Original cell 25 ----
Adipose_organized<-result_matrix

# ---- Original cell 27 ----
df<-read.csv(source_file("SCENIC_Capillary/SCENIC_Aging_Heart/output/Step2_regulonTargetsInfo.csv"), header = TRUE)


# ---- Original cell 29 ----
# TF列がtarget_stringsのいずれかを含む行を抽出
result <- df %>% filter(TF %in% unique_TFs)
result

# ---- Original cell 30 ----
# さらにその中でgenesに含まれる分子がターゲットに入っているものに絞る
result_gene <- result %>% filter(gene %in% genes_reversed)
result_gene

# ---- Original cell 31 ----
# データフレームをマトリックスに変換
coex_matrix <- as.matrix(result_gene)

# 行名と列名を設定
row_names <- genes_reversed
col_names <- unique(result_gene$TF)

# 空のマトリックスを作成
result_matrix <- matrix(0, nrow = length(row_names), ncol = length(col_names))

# 行名と列名を設定
rownames(result_matrix) <- row_names
colnames(result_matrix) <- col_names

# マトリックスにCoexWeightの値を入れる
for (i in 1:nrow(result_gene)) {
  row_index <- match(result_gene[i, "gene"], row_names)
  col_index <- match(result_gene[i, "TF"], col_names)
  result_matrix[row_index, col_index] <- result_gene[i, "CoexWeight"]
}

# 結果を表示
result_matrix

# ---- Original cell 32 ----
Heart_organized<-result_matrix

# ---- Original cell 34 ----
df<-read.csv(source_file("SCENIC_Capillary/SCENIC_Aging_Muscle/output/Step2_regulonTargetsInfo.csv"), header = TRUE)


# ---- Original cell 35 ----
# TF列がtarget_stringsのいずれかを含む行を抽出
result <- df %>% filter(TF %in% unique_TFs)
result

# ---- Original cell 36 ----
# さらにその中でgenesに含まれる分子がターゲットに入っているものに絞る
result_gene <- result %>% filter(gene %in% genes_reversed)
result_gene

# ---- Original cell 37 ----
# データフレームをマトリックスに変換
coex_matrix <- as.matrix(result_gene)

# 行名と列名を設定
row_names <- genes_reversed
col_names <- unique(result_gene$TF)

# 空のマトリックスを作成
result_matrix <- matrix(0, nrow = length(row_names), ncol = length(col_names))

# 行名と列名を設定
rownames(result_matrix) <- row_names
colnames(result_matrix) <- col_names

# マトリックスにCoexWeightの値を入れる
for (i in 1:nrow(result_gene)) {
  row_index <- match(result_gene[i, "gene"], row_names)
  col_index <- match(result_gene[i, "TF"], col_names)
  result_matrix[row_index, col_index] <- result_gene[i, "CoexWeight"]
}

# 結果を表示
result_matrix

# ---- Original cell 38 ----
Muscle_organized<-result_matrix

# ---- Original cell 40 ----
df<-read.csv(source_file("SCENIC_Capillary/SCENIC_Aging_Liver/output/Step2_regulonTargetsInfo.csv"), header = TRUE)


# ---- Original cell 41 ----
# TF列がtarget_stringsのいずれかを含む行を抽出
result <- df %>% filter(TF %in% unique_TFs)
result

# ---- Original cell 42 ----
# さらにその中でgenesに含まれる分子がターゲットに入っているものに絞る
result_gene <- result %>% filter(gene %in% genes_reversed)
result_gene

# ---- Original cell 43 ----
# データフレームをマトリックスに変換
coex_matrix <- as.matrix(result_gene)

# 行名と列名を設定
row_names <- genes_reversed
col_names <- unique(result_gene$TF)

# 空のマトリックスを作成
result_matrix <- matrix(0, nrow = length(row_names), ncol = length(col_names))

# 行名と列名を設定
rownames(result_matrix) <- row_names
colnames(result_matrix) <- col_names

# マトリックスにCoexWeightの値を入れる
for (i in 1:nrow(result_gene)) {
  row_index <- match(result_gene[i, "gene"], row_names)
  col_index <- match(result_gene[i, "TF"], col_names)
  result_matrix[row_index, col_index] <- result_gene[i, "CoexWeight"]
}

# 結果を表示
result_matrix

# ---- Original cell 44 ----
Liver_organized<-result_matrix

# ---- Original cell 45 ----
result_matrix_1<-Adipose_organized
result_matrix_2<-Heart_organized
result_matrix_3<-Muscle_organized
result_matrix_4<-Liver_organized

# ---- Original cell 47 ----
# 新しい列名リスト
new_col_names <- setdiff(unique_TFs, colnames(result_matrix_1))

# 新しい列を追加して欠損値を0で埋める
if (length(new_col_names) > 0) {
  result_matrix_1 <- cbind(result_matrix_1, matrix(0, nrow = nrow(result_matrix_1), ncol = length(new_col_names)))
  colnames(result_matrix_1)[(ncol(result_matrix_1) - length(new_col_names) + 1):ncol(result_matrix_1)] <- new_col_names
}
# 列名をABC順にソートしてマトリックスを再構築
result_matrix_1 <- result_matrix_1[, order(colnames(result_matrix_1))]

result_matrix_1

# ---- Original cell 48 ----
# 新しい列名リスト
new_col_names <- setdiff(unique_TFs, colnames(result_matrix_2))

# 新しい列を追加して欠損値を0で埋める
if (length(new_col_names) > 0) {
  result_matrix_2 <- cbind(result_matrix_2, matrix(0, nrow = nrow(result_matrix_2), ncol = length(new_col_names)))
  colnames(result_matrix_2)[(ncol(result_matrix_2) - length(new_col_names) + 1):ncol(result_matrix_2)] <- new_col_names
}

# 列名をABC順にソートしてマトリックスを再構築
result_matrix_2 <- result_matrix_2[, order(colnames(result_matrix_2))]

result_matrix_2

# ---- Original cell 49 ----
# 新しい列名リスト
new_col_names <- setdiff(unique_TFs, colnames(result_matrix_3))

# 新しい列を追加して欠損値を0で埋める
if (length(new_col_names) > 0) {
  result_matrix_3 <- cbind(result_matrix_3, matrix(0, nrow = nrow(result_matrix_3), ncol = length(new_col_names)))
  colnames(result_matrix_3)[(ncol(result_matrix_3) - length(new_col_names) + 1):ncol(result_matrix_3)] <- new_col_names
}

# 列名をABC順にソートしてマトリックスを再構築
result_matrix_3 <- result_matrix_3[, order(colnames(result_matrix_3))]

result_matrix_3

# ---- Original cell 50 ----
# 新しい列名リスト
new_col_names <- setdiff(unique_TFs, colnames(result_matrix_4))

# 新しい列を追加して欠損値を0で埋める
if (length(new_col_names) > 0) {
  result_matrix_4 <- cbind(result_matrix_4, matrix(0, nrow = nrow(result_matrix_4), ncol = length(new_col_names)))
  colnames(result_matrix_4)[(ncol(result_matrix_4) - length(new_col_names) + 1):ncol(result_matrix_4)] <- new_col_names
}

# 列名をABC順にソートしてマトリックスを再構築
result_matrix_4 <- result_matrix_4[, order(colnames(result_matrix_4))]

result_matrix_4

# ---- Original cell 51 ----
# マトリックスの数
num_matrices <- 4

# 新しいマトリックスのサイズ（行数と列数）
num_rows <- nrow(result_matrix_1)
num_cols <- ncol(result_matrix_1)

# 空のマトリックスを作成
average_matrix <- matrix(0, nrow = num_rows, ncol = num_cols)

# マトリックスごとに要素を合計して平均値を計算
for (i in 1:num_matrices) {
  matrix_name <- paste0("result_matrix_", i)
  matrix_data <- get(matrix_name)
  average_matrix <- average_matrix + matrix_data
}

average_matrix <- average_matrix / num_matrices

# 結果を表示
average_matrix

# ---- Original cell 59 ----
# アノテーションの論理データフレームを作成
create_annotations <- function(tfs, data) {
  # average_matrix の列名を取得
  all_colnames <- colnames(data)
  
  # 全ての列名に対して、リスト内のTFが含まれているかチェックし、1と0でマッピング
  annotation_vector <- sapply(all_colnames, function(column) as.numeric(column %in% tfs))
  
  # 名前付きベクトルを返す
  names(annotation_vector) <- all_colnames
  return(annotation_vector)
}

# 各組織のアノテーションを作成
ann_adipose <- create_annotations(top_TFs_adipose, average_matrix)
ann_heart <- create_annotations(top_TFs_heart, average_matrix)
ann_muscle <- create_annotations(top_TFs_muscle, average_matrix)
ann_liver <- create_annotations(top_TFs_liver, average_matrix)

# アノテーションデータフレームの統合
annotations <- data.frame(
  Adipose = ann_adipose,
  Heart = ann_heart,
  Muscle = ann_muscle,
  Liver = ann_liver
)

# アノテーションカラーの定義
ann_colors <- list(
  Adipose = c("gray92", '#F8766D'),
  Heart = c("gray92", '#7CAE00'),
  Muscle = c("gray92", '#00BFC4'),
  Liver = c("gray92", '#C77CFF')
)

# ヒートマップの描画
pheatmap::pheatmap(average_matrix, cluster_rows=TRUE, cluster_cols=TRUE,
                   col = my_palette, breaks = seq(0, 0.02, length.out = 101),
                   annotation_col = annotations, annotation_colors = ann_colors,
                  fontsize_row=4, fontsize_col=18,fontsize = 17,
         cellwidth = 15, angle_col =90, annotation_height = 3, main="Aging_CoexWeight")


# ---- Original cell 63 ----
library(pheatmap)

# ヒートマップを作成し、結果をオブジェクトに保存
heatmap_result <- pheatmap(average_matrix, 
                           cluster_rows = TRUE, 
                           cluster_cols = TRUE,
                           col = my_palette, 
                           breaks = seq(0, 0.03, length.out = 101),
                           annotation_col = annotations, 
                           annotation_colors = ann_colors,
                           fontsize_row = 4, 
                           fontsize_col = 18,
                           fontsize = 17,
                           cellwidth = 15, 
                           angle_col = 90, 
                           annotation_height = 3, 
                           main = "Aging_CoexWeight",
                           clustering_distance_rows = "euclidean",
                           clustering_distance_cols = "euclidean",
                           clustering_method = "ward.D2")

# 行のクラスタリング順序を抽出
row_order <- heatmap_result$tree_row$order

# 行の名前を順序に従って抽出
row_names_ordered <- rownames(average_matrix)[row_order]

# 順序付けられた行の名前を表示
print(row_names_ordered)


# ---- Original cell 64 ----
# Listデータフレームの行をrow_names_orderedの順に並べ替え
ordered_index <- match(row_names_ordered, List$gene)
List_ordered <- List[ordered_index, ]

# 並び替えたデータフレームを確認
print(List_ordered)

# ---- Original cell 72 ----
List_1 <- row_names_ordered[1:31]
List_1

# ---- Original cell 73 ----
List_2 <- row_names_ordered[32:45]
List_2

# ---- Original cell 74 ----
List_3 <- row_names_ordered[227:260]
List_3

# ---- Original cell 75 ----
List_4 <- row_names_ordered[261:266]
List_4
