# 安装和加载所需的包，如果包未安装则进行安装
required_packages <- c("adegenet", "ggplot2", "RColorBrewer")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# install.packages("parallel") 
library(parallel) 
library(adegenet)  # 加载 adegenet 包，用于遗传数据分析
library(ggplot2)   # 加载 ggplot2 包，用于数据可视化
library(RColorBrewer)  # 加载 RColorBrewer 包，用于补充颜色

# 定义文件路径
dfpath <- "/mnt/c/Users/Administrator/Desktop/1.fa"
group_info_path <- "/mnt/c/Users/Administrator/Desktop/1.txt" # 注意大小写，请一定保证这个数据第第一列是"ID"，第二列是"Group"

# 从 TXT 文件中读取分组信息
#* 使用逗号分割，表内必须存在`ID`列
group_info <- read.table(group_info_path, header = TRUE, sep = ",")

# 检查 group_info 的列名，确保存在 'ID' 列
if (!"ID" %in% colnames(group_info)) {
  stop("The 'group_info' data frame must contain a column named 'ID'.")
}

# 读取 fasta 文件并转换为 genlight 对象
flu <- fasta2genlight(dfpath, chunkSize = 1000, parallel = FALSE)

# 绘制 SNP 位置图
snpposi.plot(position(flu), genome.size = 4176, codon = FALSE) + theme_bw()
snpposi.plot(position(flu), genome.size = 4176, codon = TRUE) + theme_bw()

# PCA 分析
num_cores <- detectCores()
# # 不并行处理
# df.pca <- glPca(flu, nf = 3) #这一步耗费很多时间
# 并行处理
df.pca <- glPca(flu, nf = 3, parallel = TRUE, n.cores = num_cores)

df.pca.scores <- as.data.frame(df.pca$scores)
df.pca.scores$ID <- rownames(df.pca.scores)

# 计算每个主成分的标准差和方差解释率
sdev <- apply(df.pca.scores[, -ncol(df.pca.scores)], 2, sd)
variance_explained <- sdev^2 / sum(sdev^2) * 100

# 将分组信息与 PCA 结果合并
df.pca.scores <- merge(df.pca.scores, group_info, by = "ID")
df.pca.scores
# 创建一个形状向量，为每个组别分配不同的形状
unique_groups <- unique(df.pca.scores$Group)
shapes <- seq(1, length(unique_groups))

# 指定初始颜色
custom_colors <- c('#EA1F1F', '#E88421', '#E5C923', '#a5ce4c', '#35ca32', '#41fcb1', 
                   '#449657', '#176555', '#369BA8', '#2B7EBC', '#3626D1', '#A128CE', 
                   '#999999')

# 如果指定颜色不够，用 RColorBrewer 中的色盘补充
if (length(unique_groups) > length(custom_colors)) {
  additional_colors <- colorRampPalette(brewer.pal(9, "Set3"))(length(unique_groups) - length(custom_colors))
  custom_colors <- c(custom_colors, additional_colors)
}

# 定义绘图函数，带有置信区间
plot_pca <- function(data, show_id = FALSE) {
  p <- ggplot(data, aes(x = PC1, y = PC2, color = Group, shape = Group)) +  
    geom_point(size = 2) +  
    theme_bw() +  
    stat_ellipse(aes(fill = Group, color = Group), geom = "polygon", alpha = 0.2, level = 0.95) +  
    scale_shape_manual(values = shapes) +  
    scale_color_manual(values = custom_colors) +  
    scale_fill_manual(values = custom_colors)
  if (show_id) {
    p <- p + geom_text(aes(label = ID), vjust = -0.5, size = 3, show.legend = FALSE)
  }
  return(p)
}

# 绘制 PCA 结果，带 ID 和不带 ID
plot_pca(df.pca.scores, show_id = TRUE)
plot_pca(df.pca.scores, show_id = FALSE)

write.table(df.pca$scores, file = "C:/Users/victo/Desktop/result.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

# 定义绘图函数，不带置信区间
plot_pca <- function(data, show_id = FALSE) {
  p <- ggplot(data, aes(x = PC1, y = PC2, color = Group, shape = Group)) +  
    geom_point(size = 2) +  
    theme_bw() +  
    scale_shape_manual(values = shapes) +  
    scale_color_manual(values = custom_colors) +  
    scale_fill_manual(values = custom_colors)
  if (show_id) {
    p <- p + geom_text(aes(label = ID), vjust = -0.5, size = 3, show.legend = FALSE)
  }
  return(p)
}

# 绘制 PCA 结果，带 ID 和不带 ID
plot_pca(df.pca.scores, show_id = TRUE)
plot_pca(df.pca.scores, show_id = FALSE)
