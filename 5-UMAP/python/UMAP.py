#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UMAP 降维主成分分析脚本
---------------------------------------
此脚本用于：
1. 从 PLINK 输出的 .eigenvec 文件中读取样本的前 n 个主成分（PCs）。
2. 对选定的主成分运行 UMAP 降维至二维空间。
3. 计算并打印降维后两维（V1、V2）的方差及其解释度（variance explained）。
4. 将结果保存为 CSV（可选）。

使用方法：
    python hpglobal_umap.py \
        --input /path/to/hpglobal_LD_PCA.eigenvec \
        --n_pcs 10 \
        --n_neighbors 10 \
        --min_dist 0.1 \
        --output result.csv
"""
# 输入文件类似如下，原始PLINK输出的 .eigenvec 文件格式，未经过任何修改：
# SAMPLE1 SAMPLE1 0.123456 0.234567 0.345678 ...
# SAMPLE2 SAMPLE2 0.234567 0.345678 0.456789 ...
# SAMPLE3 SAMPLE3 0.345678 0.456789 0.567890 ...

import argparse
import logging
import warnings

import numpy as np
import pandas as pd
from umap import UMAP


def parse_args():
    """
    解析命令行参数。

    Returns
    -------
    argparse.Namespace
        包含输入文件路径、主成分数量、UMAP 参数和输出路径等配置项。
    """
    parser = argparse.ArgumentParser(
        description="从 .eigenvec 文件读取主成分并执行 UMAP 降维"
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="输入 .eigenvec 文件路径（PLINK 输出）"
    )
    parser.add_argument(
        "--n_pcs", "-k", type=int, default=10,
        help="用于 UMAP 的主成分数量（默认为前 10 个）"
    )
    parser.add_argument(
        "--n_neighbors", type=int, default=10,
        help="UMAP 的邻居数（n_neighbors，默认为 10）"
    )
    parser.add_argument(
        "--min_dist", type=float, default=0.1,
        help="UMAP 的最小距离（min_dist，默认为 0.1）"
    )
    parser.add_argument(
        "--output", "-o", default=None,
        help="结果保存路径（CSV 格式），若不指定则不保存"
    )
    return parser.parse_args()


def read_eigenvec(path: str) -> pd.DataFrame:
    """
    读取 PLINK 输出的 .eigenvec 文件，并添加列名。

    Parameters
    ----------
    path : str
        .eigenvec 文件的完整路径。

    Returns
    -------
    pandas.DataFrame
        列名为 ['FID', 'IID', 'PC1', 'PC2', …] 的数据框。
    """
    # 默认两个前置列为 FID、IID，后续列为各主成分
    df = pd.read_csv(path, delim_whitespace=True, header=None)
    n_pcs = df.shape[1] - 2
    cols = ['FID', 'IID'] + [f'PC{i}' for i in range(1, n_pcs + 1)]
    df.columns = cols
    logging.info(f"读取 {path}，共 {df.shape[0]} 个样本，{n_pcs} 个主成分")
    return df


def run_umap(
    X: np.ndarray,
    n_neighbors: int = 10,
    min_dist: float = 0.1,
    random_state: int = 42
) -> np.ndarray:
    """
    对输入矩阵 X 运行 UMAP 降维至二维。

    Parameters
    ----------
    X : numpy.ndarray
        输入数据矩阵，形状为 (样本数, 特征数)。
    n_neighbors : int
        UMAP 参数，邻居数量。
    min_dist : float
        UMAP 参数，最小距离。
    random_state : int
        随机种子，用于结果可重复。

    Returns
    -------
    numpy.ndarray
        降维后坐标，形状为 (样本数, 2)。
    """
    reducer = UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=2,
        metric='euclidean',
        init='spectral',
        random_state=random_state
    )
    embedding = reducer.fit_transform(X)
    logging.info("UMAP 降维完成")
    return embedding


def compute_variance_explained(embedding: np.ndarray) -> tuple:
    """
    计算降维后两维坐标的方差及解释度。

    Parameters
    ----------
    embedding : numpy.ndarray
        UMAP 输出坐标，形状为 (样本数, 2)。

    Returns
    -------
    tuple
        (v1_var, v2_var, ratio1, ratio2):
        - v1_var: V1 的方差
        - v2_var: V2 的方差
        - ratio1: V1 的解释度（百分比）
        - ratio2: V2 的解释度（百分比）
    """
    v1_var = np.var(embedding[:, 0])
    v2_var = np.var(embedding[:, 1])
    tot = v1_var + v2_var
    ratio1 = v1_var / tot * 100
    ratio2 = v2_var / tot * 100
    logging.info(
        f"V1 方差={v1_var:.6f}, V2 方差={v2_var:.6f}, "
        f"解释度 V1={ratio1:.2f}%, V2={ratio2:.2f}%"
    )
    return v1_var, v2_var, ratio1, ratio2


def main():
    """
    主函数：串联读取、降维、计算并输出结果。
    """
    # 配置日志
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s"
    )
    # 忽略 FutureWarning 与 UMAP UserWarning
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)

    args = parse_args()

    # 读取数据
    df = read_eigenvec(args.input)

    # 提取前 k 个主成分
    pcs = [f'PC{i}' for i in range(1, args.n_pcs + 1)]
    X = df.loc[:, pcs].values

    # UMAP 降维
    embedding = run_umap(
        X,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist
    )

    # 构造输出结果 DataFrame
    out_df = pd.DataFrame({
        'FID': df['FID'],
        'IID': df['IID'],
        'V1': embedding[:, 0],
        'V2': embedding[:, 1]
    })

    # 计算解释度
    v1_var, v2_var, ratio1, ratio2 = compute_variance_explained(embedding)

    # 打印并（可选）保存
    print(out_df.to_string(index=False))
    print(f"\nV1 方差: {v1_var:.6f}, V2 方差: {v2_var:.6f}")
    print(f"解释度: V1: {ratio1:.2f}%, V2: {ratio2:.2f}%")

    if args.output:
        out_df.to_csv(args.output, index=False)
        logging.info(f"结果已保存至 {args.output}")


if __name__ == "__main__":
    main()
