library(parallel)  #* 加载 parallel 包，用于并行计算
library(adegenet)  #* 加载 adegenet 包，用于遗传数据分析
library(ggplot2)   #* 加载 ggplot2 包，用于数据可视化
library(RColorBrewer)  #* 加载 RColorBrewer 包，用于补充颜色
library(patchwork)  #* 加载 patchwork 包，用于组合多个图形
library(tidyverse) #* 加载 tidyverse 包，用于数据处理和可视化
#! 注意，本代码用于运行PCA，而非DAPC
#* 定义文件路径
fasta_path <- "/home/luolintao/0-tmp/2-DAPC/2-Global-20K.aln.snp-sites.fasta"
#* 注意大小写，请一定保证这个数据第一列是"ID"，第二列是"Group"
group_info_path <- "/home/luolintao/0-tmp/2-DAPC/2-group.csv"
#* 结果输出目录
output_dir <- "/home/luolintao/0-tmp/2-DAPC/output/2-Global" #! 请单独创建一个文件夹，避免文件覆盖
#* 创建输出目录（如果不存在）
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
#* 使用逗号分割，表内必须存在`ID`列
if (!file.exists(group_info_path)) {
  stop(paste("分组信息文件不存在:", group_info_path))
}

group_info <- read.table(group_info_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# 检查 group_info 的列名，确保存在 'ID' 列
if (!"ID" %in% colnames(group_info)) {
  stop("The 'group_info' data frame must contain a column named 'ID'.")
}

if (!"Group" %in% colnames(group_info)) {
  stop("The 'group_info' data frame must contain a column named 'Group'.")
}

cat("成功读取分组信息，包含", nrow(group_info), "个样本\n")
cat("包含", length(unique(group_info$Group)), "个分组:", paste(unique(group_info$Group), collapse = ", "), "\n")

# 读取 fasta 文件并转换为 genlight 对象
if (!file.exists(fasta_path)) {
  stop(paste("FASTA文件不存在:", fasta_path))
}

cat("正在读取FASTA文件并转换为genlight对象...\n")
#* 使用 fasta2genlight 函数读取 fasta 文件
flu <- fasta2genlight(fasta_path, 
                      n.cores = NULL, #* n.cores 代表使用的CPU核心数
                      chunkSize = 100000,  #* chunksize代表每个块的大小，默认为1000
                      parallel = TRUE) #* 打开并行处理
cat("成功读取", nInd(flu), "个个体，", nLoc(flu), "个位点\n")

# 绘制 SNP 位置图
p1 <- snpposi.plot(position(flu), genome.size = 16569, codon = FALSE) + theme_bw() + theme(panel.grid = element_blank())
pX <- snpposi.plot(position(flu), genome.size = 16569, codon = TRUE) + theme_bw() + theme(panel.grid = element_blank())

# PCA 分析
num_cores <- detectCores()

# 并行处理
df.pca <- glPca(flu,
                nf = 10,  #* nf 代表主成分的数量
                center = TRUE,  #* center 代表是否中心化数据
                parallel = TRUE, #* 并行处理
                n.cores = num_cores)

df.pca.scores <- as.data.frame(df.pca$scores)
df.pca.scores$ID <- rownames(df.pca.scores)

# 计算每个主成分的标准差和方差解释率
sdev <- apply(df.pca.scores[, -ncol(df.pca.scores)], 2, sd)
variance_explained <- sdev^2 / sum(sdev^2) * 100

# 确保特征值数量与方差解释率数量一致
eigenvalues <- df.pca$eig
n_components <- length(variance_explained)
if(length(eigenvalues) > n_components) {
  eigenvalues <- eigenvalues[1:n_components]
}

# 将分组信息与 PCA 结果合并
df.pca.scores <- merge(df.pca.scores, group_info, by = "ID")
df.pca.scores |> as_tibble() -> df_pca_score
print(df_pca_score,n = 10)

# 1. 数值统计摘要
cat("\n=== PCA 分析结果摘要 ===\n")

# 特征值和累积方差解释率
# eigenvalues 已在上面处理，确保数量一致
cumulative_var <- cumsum(variance_explained)

summary_table <- data.frame(
  PC = paste0("PC", 1:length(eigenvalues)),
  Eigenvalue = round(eigenvalues, 4),
  Variance_Explained = round(variance_explained, 2),
  Cumulative_Variance = round(cumulative_var, 2)
)

print(summary_table)
write.csv(summary_table, file.path(output_dir, "pca_summary_statistics.csv"), row.names = FALSE)

# Kaiser准则分析
kaiser_threshold <- mean(eigenvalues)
significant_pcs <- sum(eigenvalues > kaiser_threshold)
cat("\nKaiser准则分析:\n")
cat("平均特征值:", round(kaiser_threshold, 4), "\n")
cat("有意义的主成分数量:", significant_pcs, "\n")

# 2. 载荷矩阵和贡献度分析
if(!is.null(df.pca$loadings)) {
  # 前几个PC的主要贡献变量
  loadings_matrix <- df.pca$loadings[, 1:min(4, ncol(df.pca$loadings))]
  loadings_df <- as.data.frame(loadings_matrix)
  loadings_df$SNP <- rownames(loadings_df)
  
  # 确保载荷值是数值型
  for(i in 1:(ncol(loadings_df)-1)) {
    loadings_df[, i] <- as.numeric(loadings_df[, i])
  }

  # 保存载荷矩阵
  write.csv(loadings_df, file.path(output_dir, "pca_loadings.csv"), row.names = FALSE)

  # 找出对PC1和PC2贡献最大的SNP
  if("PC1" %in% colnames(loadings_df)) {
    top_contributors_pc1 <- head(loadings_df[order(abs(loadings_df$PC1), decreasing = TRUE), ], 10)
    cat("\n对PC1贡献最大的10个SNP:\n")
    print(top_contributors_pc1[, c("SNP", "PC1")])
  }
  
  if("PC2" %in% colnames(loadings_df)) {
    top_contributors_pc2 <- head(loadings_df[order(abs(loadings_df$PC2), decreasing = TRUE), ], 10)
    cat("\n对PC2贡献最大的10个SNP:\n")
    print(top_contributors_pc2[, c("SNP", "PC2")])
  }
} else {
  cat("\n注意: 该PCA方法未返回载荷矩阵信息\n")
}

# 3. 群体间距离分析
group_centroids <- df_pca_score %>%
  group_by(Group) %>%
  summarise(
    PC1_mean = mean(PC1, na.rm = TRUE),
    PC2_mean = mean(PC2, na.rm = TRUE),
    PC1_sd = sd(PC1, na.rm = TRUE),
    PC2_sd = sd(PC2, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  )

print("群体质心坐标:")
print(group_centroids)
write.csv(group_centroids, file.path(output_dir, "group_centroids.csv"), row.names = FALSE)

# 计算群体间欧氏距离
if(nrow(group_centroids) > 1) {
  centroid_coords <- group_centroids[, c("PC1_mean", "PC2_mean")]
  rownames(centroid_coords) <- group_centroids$Group

  group_distances <- as.matrix(dist(centroid_coords))
  cat("\n群体间距离矩阵 (基于PC1-PC2):\n")
  print(round(group_distances, 4))
  write.csv(group_distances, file.path(output_dir, "group_distances.csv"))
}

# 4. 数据质量评估
sample_quality <- df_pca_score %>%
  mutate(
    distance_from_origin = sqrt(PC1^2 + PC2^2),
    PC1_abs = abs(PC1),
    PC2_abs = abs(PC2)
  ) %>%
  arrange(desc(distance_from_origin))

cat("\n距离原点最远的10个样本（可能的离群值）:\n")
print(sample_quality[1:min(10, nrow(sample_quality)), c("ID", "Group", "distance_from_origin")])

# 识别潜在离群值（距离超过2个标准差）
outlier_threshold <- mean(sample_quality$distance_from_origin) + 2 * sd(sample_quality$distance_from_origin)
outliers <- sample_quality[sample_quality$distance_from_origin > outlier_threshold, ]

if(nrow(outliers) > 0) {
  cat("\n潜在离群样本:\n")
  print(outliers[, c("ID", "Group", "distance_from_origin")])
  write.csv(outliers, file.path(output_dir, "potential_outliers.csv"), row.names = FALSE)
}

#* 将 PCA 结果保存为 CSV 文件
write.csv(df_pca_score, file.path(output_dir, "df_pca_score.csv"), row.names = FALSE)
# 检查合并后的数据
cat("合并后的数据维度:", dim(df.pca.scores), "\n")
cat("各组样本数量:\n")
print(table(df.pca.scores$Group))

# 创建一个形状向量，为每个组别分配不同的形状
unique_groups <- unique(df.pca.scores$Group)
num_groups <- length(unique_groups)

# 确保有足够的形状可用（R中可用的形状有限）
max_shapes <- 25  # R中可用的形状数量
if (num_groups > max_shapes) {
  # 如果组别太多，重复使用形状
  shapes <- rep(1:max_shapes, length.out = num_groups)
  warning(paste("组别数量(", num_groups, ")超过可用形状数量(", max_shapes, ")，将重复使用形状"))
} else {
  shapes <- seq(1, num_groups)
}

# 指定初始颜色
custom_colors <- c('#EA1F1F', '#E88421', '#E5C923', '#a5ce4c', '#35ca32', '#41fcb1', 
                   '#449657', '#176555', '#369BA8', '#2B7EBC', '#3626D1', '#A128CE', 
                   '#999999')

# 如果指定颜色不够，用 RColorBrewer 中的色盘补充
if (num_groups > length(custom_colors)) {
  cat("需要额外的", num_groups - length(custom_colors), "种颜色\n")
  # 使用多个调色板确保颜色足够且区分度高
  additional_colors <- c()
  palettes <- c("Set1", "Set2", "Set3", "Dark2", "Paired")
  
  for (palette in palettes) {
    if (length(additional_colors) >= (num_groups - length(custom_colors))) break
    max_colors <- brewer.pal.info[palette, "maxcolors"]
    new_colors <- brewer.pal(min(max_colors, num_groups - length(custom_colors) - length(additional_colors)), palette)
    additional_colors <- c(additional_colors, new_colors)
  }
  
  # 如果还不够，使用渐变色
  if (length(additional_colors) < (num_groups - length(custom_colors))) {
    remaining_needed <- num_groups - length(custom_colors) - length(additional_colors)
    gradient_colors <- colorRampPalette(c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7"))(remaining_needed)
    additional_colors <- c(additional_colors, gradient_colors)
  }
  
  custom_colors <- c(custom_colors, additional_colors[1:(num_groups - length(custom_colors))])
}

# 确保颜色数量与组别数量匹配
custom_colors <- custom_colors[1:num_groups]
names(custom_colors) <- unique_groups

# 调试输出：检查颜色映射
cat("组别名称:", paste(unique_groups, collapse = ", "), "\n")
cat("颜色映射检查完成\n")

# 定义绘图函数，带有置信区间（稳健版本）
plot_pca <- function(data, show_id = FALSE, show_ellipse = TRUE) {
  # 检查每个组的样本数量
  group_counts <- table(data$Group)
  groups_with_enough_points <- names(group_counts[group_counts >= 3])
  
  # 创建坐标轴标签，包含方差解释率
  x_label <- paste0("PC1 (", round(variance_explained[1], 2), "%)")
  y_label <- paste0("PC2 (", round(variance_explained[2], 2), "%)")
  
  p <- ggplot(data, aes(x = PC1, y = PC2, color = Group, shape = Group)) +  
    geom_point(size = 3) +  
    theme_bw() +  
    labs(x = x_label, y = y_label) +
    scale_shape_manual(values = shapes) +  
    scale_color_manual(values = custom_colors) +  
    scale_fill_manual(values = custom_colors) +
    theme(panel.grid = element_blank())
  
  # 只为有足够样本的组添加椭圆
  if (show_ellipse && length(groups_with_enough_points) > 0) {
    data_for_ellipse <- data[data$Group %in% groups_with_enough_points, ]
    if (nrow(data_for_ellipse) > 0) {
      # 使用 tryCatch 捕获可能的错误
      tryCatch({
        p <- p + stat_ellipse(data = data_for_ellipse, 
                              aes(fill = Group, color = Group), 
                              geom = "polygon", 
                              alpha = 0.2, 
                              level = 0.95,
                              type = "norm")  # 使用 "norm" 方法更稳定
      }, warning = function(w) {
        message("Warning in ellipse calculation: ", w$message)
      }, error = function(e) {
        message("Error in ellipse calculation, skipping ellipses: ", e$message)
      })
    }
  }
  
  if (show_id) {
    p <- p + geom_text(aes(label = ID), vjust = -0.5, size = 3, show.legend = FALSE)
  }
  return(p)
}

# 绘制 PCA 结果，带 ID 和不带 ID（第一次调用，带置信椭圆）
# 输出各组样本数量信息
group_counts <- table(df.pca.scores$Group)
cat("各组样本数量:\n")
print(group_counts)
cat("\n样本数少于3的组将不显示置信椭圆\n")

p3 <- plot_pca(df.pca.scores, show_id = TRUE, show_ellipse = TRUE)
p4 <- plot_pca(df.pca.scores, show_id = FALSE, show_ellipse = TRUE)


# 定义备用绘图函数，不带置信区间（如果需要简化显示）
plot_pca_simple <- function(data, show_id = FALSE) {
  # 创建坐标轴标签，包含方差解释率
  x_label <- paste0("PC1 (", round(variance_explained[1], 2), "%)")
  y_label <- paste0("PC2 (", round(variance_explained[2], 2), "%)")
  
  p <- ggplot(data, aes(x = PC1, y = PC2, color = Group, shape = Group)) +  
    geom_point(size = 3) +  
    theme_bw() +  
    labs(x = x_label, y = y_label) +
    scale_shape_manual(values = shapes) +  
    scale_color_manual(values = custom_colors) +  
    scale_fill_manual(values = custom_colors) +
    theme(panel.grid = element_blank())
  if (show_id) {
    p <- p + geom_text(aes(label = ID), vjust = -0.5, size = 3, show.legend = FALSE)
  }
  return(p)
}

plot_pca_3_4_simple <- function(data, show_id = FALSE) {
  # 创建坐标轴标签，包含方差解释率
  x_label <- paste0("PC3 (", round(variance_explained[3], 2), "%)")
  y_label <- paste0("PC4 (", round(variance_explained[4], 2), "%)")

  p <- ggplot(data, aes(x = PC3, y = PC4, color = Group, shape = Group)) +  
    geom_point(size = 3) +  
    theme_bw() +  
    labs(x = x_label, y = y_label) +
    scale_shape_manual(values = shapes) +  
    scale_color_manual(values = custom_colors) +  
    scale_fill_manual(values = custom_colors) +
    theme(panel.grid = element_blank())
  if (show_id) {
    p <- p + geom_text(aes(label = ID), vjust = -0.5, size = 3, show.legend = FALSE)
  }
  return(p)
}


# 绘制 PCA 结果，不带置信区间（第二次调用，简化版）
cat("\n不带置信椭圆的简化版本:\n")
p5 <- plot_pca_simple(df.pca.scores, show_id = TRUE)
p6 <- plot_pca_simple(df.pca.scores, show_id = FALSE)
p7 <- plot_pca_3_4_simple(df.pca.scores, show_id = TRUE)
p8 <- plot_pca_3_4_simple(df.pca.scores, show_id = FALSE)

# 绘制variance_explained的柱状图
p2 <- variance_explained |> 
  as_tibble() |>
  mutate(PC = row_number()) |>
  select(PC, value = value) |>
  ggplot(aes(x = factor(PC), y = value)) +
  geom_bar(stat = "identity", fill = "#228c93") +
  labs(title = "Variance explained by each principal component",
       x = "Principal components",
       y = "Variance explained (%)") +
  geom_text(aes(label = round(value, 2)), vjust = -0.5, size = 3 ,color ="#228c93") +
  theme_bw() + 
  theme(panel.grid = element_blank())

# 5. 增强的可视化输出
# 添加碎石图
scree_data <- summary_table[1:min(10, nrow(summary_table)), ]
scree_data$PC_num <- as.numeric(gsub("PC", "", scree_data$PC))  # 提取PC编号用于排序
scree_data <- scree_data[order(scree_data$PC_num), ]  # 按PC编号排序

scree_plot <- ggplot(scree_data, 
                     aes(x = factor(PC, levels = PC), y = Variance_Explained)) +
  geom_point(size = 3, color = "#228c93") +
  geom_line(group = 1, color = "#228c93") +
  geom_hline(yintercept = kaiser_threshold/sum(eigenvalues)*100, 
             linetype = "dashed", color = "red", alpha = 0.7) +
  labs(title = "Scree Plot", 
       x = "Principal Components", 
       y = "Variance Explained (%)") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

# 双标图（如果需要显示变量载荷）
biplot <- NULL
if(!is.null(df.pca$loadings) && ncol(df.pca$loadings) >= 2) {
  biplot_data <- data.frame(
    PC1 = df.pca$loadings[, 1] * 5, # 缩放因子
    PC2 = df.pca$loadings[, 2] * 5,
    SNP = rownames(df.pca$loadings)
  )

  # 只显示贡献最大的变量
  top_vars <- head(biplot_data[order(sqrt(biplot_data$PC1^2 + biplot_data$PC2^2), decreasing = TRUE), ], 20)

  biplot <- p6 + 
    geom_segment(data = top_vars, 
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.1, "cm")), 
                 color = "red", alpha = 0.6, inherit.aes = FALSE) +
    geom_text(data = top_vars, 
              aes(x = PC1*1.1, y = PC2*1.1, label = SNP),
              size = 2, color = "red", inherit.aes = FALSE) +
    labs(title = "PCA Biplot with Variable Loadings")
}

# 使用 patchwork 组合图形
if(!is.null(biplot)) {
  # 包含双标图的完整布局
  p_all <-  (p1 + p2 + scree_plot + p5 + p6 + biplot + p7 + p8) + 
    plot_layout(ncol = 2) + 
    plot_annotation(title = "PCA analysis", 
                    subtitle = "Comprehensive Analysis with Enhanced Visualizations",
                    tag_levels = "A") 
} else {
  # 不包含双标图的布局
  p_all <-  (p1 + p2 + scree_plot + p5 + p6 + p7 + p8) + 
    plot_layout(ncol = 2) + 
    plot_annotation(title = "PCA analysis", 
                    subtitle = "Comprehensive Analysis with Enhanced Visualizations",
                    tag_levels = "A") 
}

p_all

# 6. 完整的报告生成
report <- list(
  analysis_date = Sys.time(),
  input_file = fasta_path,
  n_individuals = nInd(flu),
  n_loci = nLoc(flu),
  n_groups = length(unique(df_pca_score$Group)),
  group_names = paste(unique(df_pca_score$Group), collapse = ", "),
  variance_pc1_pc2 = sum(variance_explained[1:2]),
  significant_pcs = significant_pcs,
  potential_outliers = if(exists("outliers")) nrow(outliers) else 0
)

# 保存报告
report_text <- paste(names(report), ":", report, collapse = "\n")
writeLines(report_text, file.path(output_dir, "pca_analysis_report.txt"))

# 保存图形
ggsave(p_all,
       filename = file.path(output_dir, "pca_comprehensive_plot.pdf"),
       height = 16,
       width = 12)

# 单独保存碎石图
ggsave(scree_plot,
       filename = file.path(output_dir, "scree_plot.pdf"),
       height = 6,
       width = 8)

# 如果有双标图，单独保存
if(!is.null(biplot)) {
  ggsave(biplot,
         filename = file.path(output_dir, "pca_biplot.pdf"),
         height = 8,
         width = 10)
}

cat("\n=== 分析完成 ===\n")
cat("所有结果文件已保存到目录:", output_dir, "\n")
cat("生成的文件:\n")
cat("- PCA scores: df_pca_score.csv\n")
cat("- 统计摘要: pca_summary_statistics.csv\n")
cat("- 群体质心: group_centroids.csv\n")
cat("- 群体距离: group_distances.csv\n")
cat("- 分析报告: pca_analysis_report.txt\n")
cat("- 综合图形: pca_comprehensive_plot.pdf\n")
cat("- 碎石图: scree_plot.pdf\n")
if(!is.null(biplot)) {
  cat("- 双标图: pca_biplot.pdf\n")
}
if(exists("outliers") && nrow(outliers) > 0) {
  cat("- 离群样本: potential_outliers.csv\n")
}
if(!is.null(df.pca$loadings)) {
  cat("- 载荷矩阵: pca_loadings.csv\n")
}
