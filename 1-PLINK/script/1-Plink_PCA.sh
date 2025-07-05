#!/bin/bash

# 脚本名称: simple_pca.sh
# 功能: 使用 PLINK (v1.9) 对 VCF 文件进行主成分分析（PCA）
# 输入文件: /mnt/d/幽门螺旋杆菌/Script/分析结果/merged_vcf/filtered.vcf.gz

# 设置错误处理
set -e

# 定义输入和输出路径
INPUT_VCF="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/merged_biallelic_7544_re.NoN.vcf.gz"
OUTPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/3-PCA/7544_整体"
mkdir -p "$OUTPUT_DIR"
PLINK_PREFIX="${OUTPUT_DIR}/filtered"
PCA_OUTPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/3-PCA/7544_整体"
mkdir -p "$PCA_OUTPUT_DIR"

# 将 VCF 转换为 PLINK 二进制格式
echo "将 VCF 转换为 PLINK 二进制格式"
plink --vcf "$INPUT_VCF" --make-bed \
      --double-id  \
      --allow-extra-chr \
      --out "$PLINK_PREFIX"

# 进行连锁不平衡（LD）修剪
echo "进行连锁不平衡（LD）修剪"
plink --bfile "$PLINK_PREFIX" \
      --indep-pairwise 50 10 0.1 \
	--allow-extra-chr \
      --out "${PLINK_PREFIX}_ld_pruned"

# 运行主成分分析（PCA）
echo "运行主成分分析（PCA）"
plink --bfile "$PLINK_PREFIX" \
      --extract "${PLINK_PREFIX}_ld_pruned.prune.in" \
      --pca 20 \
	  --allow-extra-chr \
      --out "${PCA_OUTPUT_DIR}/hpglobal_LD_PCA"

echo "PCA 分析完成！所有输出文件保存在 ${PCA_OUTPUT_DIR}/ 目录下。"
