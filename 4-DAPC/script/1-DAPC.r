library(parallel)  #* 加载 parallel 包，用于并行计算
library(adegenet)  #* 加载 adegenet 包，用于遗传数据分析
library(ggplot2)   #* 加载 ggplot2 包，用于数据可视化
library(RColorBrewer)  #* 加载 RColorBrewer 包，用于补充颜色
library(patchwork)  #* 加载 patchwork 包，用于组合多个图形
library(tidyverse) #* 加载 tidyverse 包，用于数据处理和可视化

#* 定义文件路径
fasta_path <- "/mnt/f/OneDrive/文档（科研）/脚本/Download/2-Dimensionality-analysis/4-DAPC/example/example.aln.fasta"
#* 注意大小写，请一定保证这个数据第一列是"ID"，第二列是"Group"
group_info_path <- "/mnt/f/OneDrive/文档（科研）/脚本/Download/2-Dimensionality-analysis/4-DAPC/conf/group.csv"
#* 结果csv文件
result_csv_path <- "/mnt/c/Users/Administrator/Desktop/df_pca_score.csv"
#* 绘图结果保存路径
pdf_file_path <- "/mnt/c/Users/Administrator/Desktop/pca_plot.pdf"
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
                      chunkSize = 1000,  #* chunksize代表每个块的大小，默认为1000
                      parallel = TRUE) #* 打开并行处理
cat("成功读取", nInd(flu), "个个体，", nLoc(flu), "个位点\n")

# 绘制 SNP 位置图
p1 <- snpposi.plot(position(flu), genome.size = 4176, codon = FALSE) + theme_bw() + theme(panel.grid = element_blank())
pX <- snpposi.plot(position(flu), genome.size = 4176, codon = TRUE) + theme_bw() + theme(panel.grid = element_blank())

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

# 将分组信息与 PCA 结果合并
df.pca.scores <- merge(df.pca.scores, group_info, by = "ID")
df.pca.scores |> as_tibble() -> df_pca_score
print(df_pca_score,n = 10)
#* 将 PCA 结果保存为 CSV 文件
write.csv(df_pca_score, file = result_csv_path, row.names = FALSE)
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

# 使用 patchwork 组合图形
p_all <-  (p1 + p2 + p5 + p6 + p7 + p8) + 
  plot_layout(ncol = 2) + 
  plot_annotation(title = "DAPC analysis", 
                  subtitle = "With and Without Confidence Ellipses",
                  tag_levels = "A") 
p_all                  
ggsave(p_all ,
       filename = pdf_file_path ,
       height = 12,
       width = 10)
