# ============================================================================
# PCA结果散点图绘制
# 使用传统R基础绘图函数，仿照经典的PCA可视化风格
# ============================================================================

# 加载必要的包
library(tidyverse)
library(RColorBrewer)
# ============================================================================
# 1. 数据读取和预处理
# ============================================================================

# 读取PCA结果数据
cat("正在读取PCA数据...\n")
# !注意：请根据实际数据路径修改以下路径
data <- read.csv("/mnt/f/OneDrive/文档（科研）/脚本/Download/2-Dimensionality-analysis/3-PCA可视化/data/example_pca.csv", header = TRUE)

# 数据基本信息
cat("数据概览:\n")
cat("样本数量:", nrow(data), "\n")
cat("变量数量:", ncol(data), "\n")
print(str(data))

# 提取分类信息
Class_big_levels <- unique(data$Class_big) #todo 注意csv文件中Class_big和Class_small的列名
class_small_levels <- unique(data$Class_small) #todo 注意csv文件中Class_big和Class_small的列名

cat("\nClass_big类别 (", length(Class_big_levels), "个):", paste(Class_big_levels, collapse=", "), "\n")
cat("Class_small类别 (", length(class_small_levels), "个):", paste(class_small_levels, collapse=", "), "\n")

# ============================================================================
# 2. 数据准备
# ============================================================================

# 提取主成分数据
PC1 <- data$PC1
PC2 <- data$PC2  
PC3 <- data$PC3

# 创建完整的数据框
frame <- data.frame(
  ID = data$ID,
  PC1 = PC1,
  PC2 = PC2,
  PC3 = PC3, #! 如果没有PC3数据，可以注释掉这一行
  Class_big = data$Class_big,
  Class_small = data$Class_small
)

cat("\n数据框创建完成，包含", nrow(frame), "个样本\n")

# ============================================================================
# 3. 颜色和形状映射设置
# ============================================================================


# ========== 新增：从 color.csv 读取颜色映射 ===========
color_map_file <- "/mnt/f/OneDrive/文档（科研）/脚本/Download/2-Dimensionality-analysis/3-PCA可视化/conf/color.csv"
color_df <- read.csv(color_map_file, header = TRUE, stringsAsFactors = FALSE)
color_map_from_file <- setNames(color_df$color, color_df$Class_big)
# 自动分配颜色的调色板
auto_palette <- brewer.pal(max(8, length(Class_big_levels)), "Set2")
if (length(Class_big_levels) > length(auto_palette)) {
    auto_palette <- colorRampPalette(brewer.pal(8, "Set2"))(length(Class_big_levels))
}
auto_palette_idx <- 1

color_mapping <- c()
for (big_class in Class_big_levels) {
  if (!is.na(color_map_from_file[big_class])) {
    color_mapping[big_class] <- color_map_from_file[big_class]
  } else {
    # 自动分配颜色
    color_mapping[big_class] <- auto_palette[auto_palette_idx]
    auto_palette_idx <- auto_palette_idx + 1
    cat(sprintf("%s类别没有分配颜色，请注意\n", big_class))
  }
}

# 定义形状值（R基础绘图支持的点形状）
shape_palette <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
                  16, 17, 18, 19, 20, 21, 22, 23, 24, 25)

# 为每个Class_big内的Class_small分配形状
create_shape_mapping <- function(data, big_levels, shape_values) {
  shape_map <- c()
  
  for (big_class in big_levels) {
    small_classes <- unique(data$Class_small[data$Class_big == big_class])
    small_classes <- sort(small_classes)  # 排序确保一致性
    
    # 每个大类别内部重新从第1个形状开始分配
    shape_idx <- 1
    
    for (small_class in small_classes) {
      shape_map[small_class] <- shape_values[shape_idx]
      shape_idx <- shape_idx + 1
      if (shape_idx > length(shape_values)) shape_idx <- 1
    }
  }
  return(shape_map)
}

shape_mapping <- create_shape_mapping(data, Class_big_levels, shape_palette)

# 打印映射信息
cat("\n=== 颜色映射 ===\n")
for (i in seq_along(color_mapping)) {
  cat(sprintf("%-12s: %s\n", names(color_mapping)[i], color_mapping[i]))
}

cat("\n=== 形状映射 ===\n")
for (i in seq_along(shape_mapping)) {
  cat(sprintf("%-20s: %2d\n", names(shape_mapping)[i], shape_mapping[i]))
}

# ============================================================================
# 4. 图形参数设置
# ============================================================================


# 定义绘图参数
plot_params <- list(
  # 图形尺寸
  pdf_width = 35, #* PDF宽度
  pdf_height = 35, #* PDF高度
  
  # 坐标轴参数
  axis_text_size = 2,      # 坐标轴刻度文字大小
  axis_label_size = 3,     # 坐标轴标签大小
  title_size = 3,          # 标题大小
  
  # 散点参数
  point_size = 3,          # 散点大小
  point_alpha = 0.8,         # 散点透明度（通过颜色实现）
  point_lwd = 3,             # 散点描边粗细
  
  # 图例参数
  legend_text_size = 2,    # 图例文字大小
  legend_point_size = 2,   # 图例中点的大小
  legend_ncol = 5            # 图例列数
)

cat("\n=== 绘图参数 ===\n")
cat("散点大小:", plot_params$point_size, "\n")
cat("坐标轴文字大小:", plot_params$axis_text_size, "\n")
cat("坐标轴标签大小:", plot_params$axis_label_size, "\n")

# ============================================================================
# 5. 主图绘制
# ============================================================================

cat("\n开始绘制PCA散点图...\n")

# 创建PDF文件
#! 注意：请根据实际数据路径修改以下路径
pdf("/mnt/f/OneDrive/文档（科研）/脚本/Download/2-Dimensionality-analysis/3-PCA可视化/output/PCA_R_plot.pdf", width = plot_params$pdf_width, height = plot_params$pdf_height)

# 设置图形布局：上方主图，下方图例
layout(matrix(c(1, 2)), widths = c(1, 1), heights = c(4, 1))

# === 绘制主散点图 ===
par(mar = c(5, 5, 4, 2))  # 增加边距以容纳更大的文字

# 创建空白图形框架
plot(PC1, PC2, 
     type = "n",                           # 不绘制点，只创建框架
     xlab = "PC1",                         # X轴标签
     ylab = "PC2",                         # Y轴标签
     main = "PCA",           # 标题
     cex.axis = plot_params$axis_text_size,    # 坐标轴刻度文字大小
     cex.lab = plot_params$axis_label_size,    # 坐标轴标签大小
     cex.main = plot_params$title_size,        # 标题大小
     mgp = c(3, 1, 0))                     # 坐标轴位置调整



# 按分类绘制散点

# 按分类绘制散点，支持描边粗细参数
draw_points_by_class <- function(frame, big_levels, color_map, shape_map, point_size, point_lwd = 1) {
  points_drawn <- 0
  
  for (big_class in big_levels) {
    small_classes <- unique(frame$Class_small[frame$Class_big == big_class])
    small_classes <- sort(small_classes)
    
    cat("绘制", big_class, "类别，包含", length(small_classes), "个子类别\n")
    
    for (small_class in small_classes) {
      # 筛选当前类别的数据
      current_data <- subset(frame, Class_big == big_class & Class_small == small_class)
      
      if (nrow(current_data) > 0) {
        # 绘制散点
        points(current_data$PC1, current_data$PC2,
               pch = shape_map[small_class],      # 点形状
               bg = color_map[big_class],         # 填充颜色
               col = color_map[big_class],        # 边框颜色
               cex = point_size,                  # 点大小
               lwd = point_lwd)                   # 描边粗细
        
        points_drawn <- points_drawn + nrow(current_data)
        cat("  -", small_class, ":", nrow(current_data), "个点\n")
      }
    }
  }
  
  cat("总共绘制了", points_drawn, "个散点\n")
}

# 执行散点绘制，增加描边粗细参数
draw_points_by_class(frame, Class_big_levels, color_mapping, shape_mapping, plot_params$point_size, plot_params$point_lwd)

# ============================================================================
# 6. 图例绘制
# ============================================================================

# === 绘制图例 ===
par(mar = c(1, 1, 1, 1))  # 减少图例区域的边距
plot.new()

# 创建图例数据
create_legend_data <- function(data, big_levels, color_map, shape_map) {
  legend_items <- list(
    labels = c(),
    colors = c(),
    shapes = c(),
    fonts = c()
  )
  
  for (big_class in big_levels) {
    # 添加大类别标题
    legend_items$labels <- c(legend_items$labels, big_class)
    legend_items$colors <- c(legend_items$colors, NA)
    legend_items$shapes <- c(legend_items$shapes, NA)
    legend_items$fonts <- c(legend_items$fonts, 2)  # 粗体
    
    # 添加小类别
    small_classes <- unique(data$Class_small[data$Class_big == big_class])
    small_classes <- sort(small_classes)
    
    for (small_class in small_classes) {
      legend_items$labels <- c(legend_items$labels, paste("  ", small_class))
      legend_items$colors <- c(legend_items$colors, color_map[big_class])
      legend_items$shapes <- c(legend_items$shapes, shape_map[small_class])
      legend_items$fonts <- c(legend_items$fonts, 1)  # 正常字体
    }
  }
  
  return(legend_items)
}

legend_data <- create_legend_data(data, Class_big_levels, color_mapping, shape_mapping)

# 绘制图例
legend("center", 
       legend = legend_data$labels,
       pch = legend_data$shapes,
       pt.bg = legend_data$colors,
       col = legend_data$colors,
       ncol = plot_params$legend_ncol,
       cex = plot_params$legend_text_size,
       pt.cex = plot_params$legend_point_size,
       bty = "n",                          # 无边框
       text.font = legend_data$fonts,
       xpd = TRUE)

# 关闭PDF设备
dev.off()

cat("\nPDF图形已保存为: PCA_plot_菌株分类.pdf\n")

# ============================================================================
# 7. 数据汇总和导出
# ============================================================================

# 创建详细的分类信息表
create_classification_summary <- function(data, color_map, shape_map) {
  summary_table <- data %>%
    select(ID, Class_big, Class_small, PC1, PC2, PC3) %>%
    arrange(Class_big, Class_small, ID) %>%
    mutate(
      Color = color_map[Class_big],
      Shape = shape_map[Class_small],
      Color_Name = names(color_map)[match(Class_big, names(color_map))]
    )
  
  return(summary_table)
}

classification_summary <- create_classification_summary(data, color_mapping, shape_mapping)

# 打印汇总信息
cat("\n=== 分类汇总信息 ===\n")
summary_stats <- classification_summary %>%
  group_by(Class_big, Class_small) %>%
  summarise(
    Count = n(),
    Color = first(Color),
    Shape = first(Shape),
    .groups = "drop"
  ) %>%
  arrange(Class_big, Class_small)

print(summary_stats)

