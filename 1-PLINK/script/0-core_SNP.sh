# 定义输入输出路径 
INPUT_VCF="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/merged_biallelic_N_correct.vcf.gz"
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