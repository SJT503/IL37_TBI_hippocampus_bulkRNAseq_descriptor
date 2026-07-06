"""
Novogene 官方 QC / mapping 真实统计 (ground-truth, 已从报告核实)
来源:
  data_table.xls   -> reads / Q20 / Q30 / GC
  align_pct.xls    -> total_map / unique_map / multi_map
  align_region.xls -> exon / intron / intergenic
  insert_mean.xls  -> insert size
禁止修改数值; 任何绘图/表格只读本文件.
"""

SAMPLES = ["sham_1","sham_2","sham_3",
           "RFP_CCI_1","RFP_CCI_2","RFP_CCI_3",
           "IL37_RFP_CCI_1","IL37_RFP_CCI_2","IL37_RFP_CCI_3"]

GROUP = {
    "sham_1":"Sham","sham_2":"Sham","sham_3":"Sham",
    "RFP_CCI_1":"RFP_CCI","RFP_CCI_2":"RFP_CCI","RFP_CCI_3":"RFP_CCI",
    "IL37_RFP_CCI_1":"IL37_RFP_CCI","IL37_RFP_CCI_2":"IL37_RFP_CCI","IL37_RFP_CCI_3":"IL37_RFP_CCI",
}

# raw_reads, clean_reads, error_rate(%), Q20(%), Q30(%), GC(%)
QC = {
    "sham_1":          dict(raw=42412330, clean=41155714, err=0.01, q20=99.24, q30=97.24, gc=48.20),
    "sham_2":          dict(raw=41634586, clean=40313282, err=0.01, q20=99.25, q30=97.23, gc=49.21),
    "sham_3":          dict(raw=41249834, clean=40154324, err=0.01, q20=99.27, q30=97.30, gc=48.89),
    "RFP_CCI_1":       dict(raw=41695500, clean=40661774, err=0.01, q20=99.26, q30=97.21, gc=48.68),
    "RFP_CCI_2":       dict(raw=41507996, clean=40569548, err=0.01, q20=99.29, q30=97.26, gc=49.20),
    "RFP_CCI_3":       dict(raw=41997624, clean=41133008, err=0.01, q20=99.31, q30=97.33, gc=49.69),
    "IL37_RFP_CCI_1":  dict(raw=42213154, clean=41397458, err=0.01, q20=99.26, q30=97.26, gc=48.26),
    "IL37_RFP_CCI_2":  dict(raw=42996972, clean=42238806, err=0.01, q20=99.26, q30=97.29, gc=48.24),
    "IL37_RFP_CCI_3":  dict(raw=46363964, clean=45388410, err=0.01, q20=99.19, q30=96.47, gc=49.11),
}

# total_map%, unique_map%, multi_map%
MAP = {
    "sham_1":          dict(total=88.18, unique=86.03, multi=2.15),
    "sham_2":          dict(total=90.94, unique=88.75, multi=2.19),
    "sham_3":          dict(total=89.82, unique=87.64, multi=2.18),
    "RFP_CCI_1":       dict(total=88.91, unique=86.63, multi=2.28),
    "RFP_CCI_2":       dict(total=91.10, unique=88.60, multi=2.50),
    "RFP_CCI_3":       dict(total=91.73, unique=89.38, multi=2.35),
    "IL37_RFP_CCI_1":  dict(total=89.33, unique=86.67, multi=2.67),
    "IL37_RFP_CCI_2":  dict(total=89.28, unique=87.11, multi=2.17),
    "IL37_RFP_CCI_3":  dict(total=91.03, unique=89.01, multi=2.02),
}

# exon%, intron%, intergenic%
REGION = {
    "sham_1":          dict(exon=87.8651, intron=5.3069, intergenic=6.8279),
    "sham_2":          dict(exon=87.1542, intron=6.6285, intergenic=6.2172),
    "sham_3":          dict(exon=87.0851, intron=6.2946, intergenic=6.6204),
    "RFP_CCI_1":       dict(exon=86.1567, intron=7.0859, intergenic=6.7575),
    "RFP_CCI_2":       dict(exon=88.8317, intron=5.0683, intergenic=6.1000),
    "RFP_CCI_3":       dict(exon=88.9850, intron=5.2727, intergenic=5.7423),
    "IL37_RFP_CCI_1":  dict(exon=87.9906, intron=4.1419, intergenic=7.8674),
    "IL37_RFP_CCI_2":  dict(exon=87.6218, intron=5.6213, intergenic=6.7570),
    "IL37_RFP_CCI_3":  dict(exon=82.6201, intron=10.9609, intergenic=6.4189),
}

# insert size mean, sd
INSERT = {
    "sham_1":          dict(mean=297.38, sd=69.81),
    "sham_2":          dict(mean=297.91, sd=70.48),
    "sham_3":          dict(mean=302.77, sd=69.92),
    "RFP_CCI_1":       dict(mean=313.84, sd=72.76),
    "RFP_CCI_2":       dict(mean=314.72, sd=72.91),
    "RFP_CCI_3":       dict(mean=307.14, sd=70.08),
    "IL37_RFP_CCI_1":  dict(mean=302.58, sd=68.30),
    "IL37_RFP_CCI_2":  dict(mean=308.09, sd=69.54),
    "IL37_RFP_CCI_3":  dict(mean=287.80, sd=63.56),
}

# Novogene 官方 DEG 概况 (DESeq2 pvalue<=0.05 & |log2FC|>=1.0)
DEG_OFFICIAL = {
    "RFP_CCI_vs_Sham":       dict(total=2094, up=1751, down=343),
    "IL37_RFP_CCI_vs_RFP_CCI":dict(total=2175, up=449,  down=1726),
    "IL37_RFP_CCI_vs_Sham":  dict(total=293,  up=147,  down=146),
}

# 流程与参考 (以 Novogene 报告为准, 覆盖 JNI 手稿的 mm10 占位笔误)
PIPELINE = dict(
    platform="Illumina NovaSeq 6000",
    read_layout="paired-end 150 bp",
    aligner="HISAT2 v2.0.5",
    quantifier="featureCounts (Subread)",
    reference="Ensembl release 107, Mus musculus C57BL/6NJ (mus_musculus_c57bl6nj)",
    vendor="Novogene (Beijing, China)",
    library_batch="X101SC25026107-Z01-J003",
)

# PCA (Novogene PCA2D.xls) — variance explained
PCA_VAR = dict(PC1=41.05, PC2=20.49, PC3=12.29)
