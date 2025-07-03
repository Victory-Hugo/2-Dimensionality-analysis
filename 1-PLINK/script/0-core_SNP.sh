: '
脚本功能说明：
本脚本用于对VCF文件中的SNP位点进行过滤和统计，主要包括以下步骤：

1. 定义输入（INPUT_VCF）和输出（OUTPUT_PREFIX）文件路径。
2. 使用VCFtools对VCF文件进行过滤：
    - 去除次等位基因频率（MAF）小于0.01的位点（即去除超低频SNP）。
    - 去除缺失率大于1%的位点（即去除大范围缺失的位点）。
    - 输出过滤后的VCF文件。
3. 对过滤后的VCF文件进行压缩（bgzip）和索引（tabix），便于后续分析。
4. 使用bcftools统计原始和过滤后VCF文件的变异信息，输出统计结果。

适用场景：
适用于群体遗传学分析前的SNP数据预处理，确保下游分析数据质量。
'

# 定义输入输出路径 
INPUT_VCF="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/\output/merge/merged_biallelic_N_correct.vcf.gz"
OUTPUT_PREFIX="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/merged_biallelic_N_correct.filtered.maf0.01.miss0.99"

# 运行 VCFtools 进行过滤
#! 去掉 MAF < 0.01 的位点(去掉超低频 SNP)
#! 去掉缺失率 > 1% 的位点(去掉大范围缺失的位点)
vcftools --gzvcf $INPUT_VCF \
         --maf 0.01 \
         --max-missing 0.99 \
         --recode \
         --recode-INFO-all \
         --out $OUTPUT_PREFIX

# 压缩和索引过滤后的 VCF
bgzip ${OUTPUT_PREFIX}.recode.vcf
tabix -p vcf ${OUTPUT_PREFIX}.recode.vcf.gz

#todo 统计变异数量
bcftools stats ${INPUT_VCF} > ${INPUT_VCF}.stats.txt
bcftools stats ${OUTPUT_PREFIX}.recode.vcf.gz > ${OUTPUT_PREFIX}.stats.txt