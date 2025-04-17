
import pandas as pd
import bioframe as bf
import re


def parse_attribute(attr_str):
    """
    将单个属性字符串解析为 pd.Series，格式：
    - 索引（name）: 引号外的字段（如 gene_id）
    - 值: 引号内的内容（如 ENSG00000243485）
    
    参数:
        attr_str: 输入字符串，例如 'gene_id "ENSG00000243485"; gene_version "5"; ...'
    
    返回:
        pd.Series
    """
    # 使用正则表达式提取键值对
    matches = re.findall(r'(\w+)\s+"([^"]+)"', attr_str)
    
    # 转换为字典 → 再转为 Series
    return pd.Series(
        data=[value for _, value in matches],  # 值（引号内）
        index=[key for key, _ in matches],     # 索引（引号外）
        name="attributes"  # Series 的名称（可选）
    )

def gene_region(GTF, gid="gene_id", gname="gene_name"):
    GTF = GTF[GTF[2]=="gene"]
    GTF[9] = None
    GTF[10] = None
    GTF_cp = GTF.copy()
    
    for i, row in GTF_cp.iterrows():
        attrs = parse_attribute(row[8])
        GTF.loc[i, 9] = attrs[gid]
        if gname in attrs:
            GTF.loc[i, 10] = attrs[gname]
        else:
            GTF.loc[i, 10] = attrs[gid]

    gene = GTF[[0,3,4,9,10,6]]
    gene.columns = ["chrom", "start", "end", "gene_id", "gene_name", "strand"]
    gene['chrom'] = gene['chrom'].astype(str)
    return gene

def transcript_region(GTF, gid="gene_id", tid="transcript_id"):
    GTF = GTF[GTF[2]=="transcript"]
    GTF[9] = None
    GTF[10] = None
    GTF_cp = GTF.copy()
    
    for i, row in GTF_cp.iterrows():
        attrs = parse_attribute(row[8])
        GTF.loc[i, 9] = attrs[gid]
        GTF.loc[i, 10] = attrs[tid]

    transcript = GTF[[0,3,4,10,9,6]]
    transcript.columns = ["chrom", "start", "end", "transcript_id", "gene_id", "strand"]
    transcript['chrom'] = transcript['chrom'].astype(str)
    return transcript


def exon_region(GTF, gid="gene_id", tid="transcript_id"):
    GTF = GTF[GTF[2]=="exon"]
    GTF[9] = None
    GTF[10] = "."
    GTF_cp = GTF.copy()
    for i, row in GTF_cp.iterrows():
        attrs = parse_attribute(row[8])
        GTF.loc[i, 9] = attrs[gid]
        GTF.loc[i, 10] = attrs[tid]

    exon = GTF[[0,3,4,10,9,6]]
    exon.columns = ["chrom", "start", "end", "transcript_id", "gene_id", "strand"]
    exon['chrom'] = exon['chrom'].astype(str)
    return exon


