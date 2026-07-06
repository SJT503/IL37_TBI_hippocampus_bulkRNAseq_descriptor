"""
扩展分析真实统计 (ground-truth, 来自诺禾结果包, 已核实).
仅客观描述数据规模, 不做机制解读 (数据论文红线).

来源:
  7.AS/1.daslist/*        rMATS 5类可变剪接事件 (全量 + 显著)
  8.SNP/1.snpsite/all_*   SNP/Indel 合集
  8.SNP/3.stat/*          逐样本 SNP/Indel 分类统计
  3.Assemble/*novel*      新转录本/新基因
"""

SAMPLES = ["sham_1","sham_2","sham_3","RFP_CCI_1","RFP_CCI_2","RFP_CCI_3",
           "IL37_RFP_CCI_1","IL37_RFP_CCI_2","IL37_RFP_CCI_3"]

# ---- 7.AS rMATS 可变剪接事件 (全量事件数; sig = FDR<0.05 显著数) ----
# 事件类型: SE跳跃外显子 MXE互斥外显子 A5SS/A3SS 5'/3'可变剪接位点 RI内含子保留
AS_EVENTS = {
    "IL37_vs_CCI":  dict(SE=(23454,357), MXE=(2201,15), A5SS=(1622,15), A3SS=(2900,10), RI=(2322,14)),
    "IL37_vs_Sham": dict(SE=(22902,36),  MXE=(2164,5),  A5SS=(1603,7),  A3SS=(2872,8),  RI=(2293,11)),
    "CCI_vs_Sham":  dict(SE=(23763,380), MXE=(2279,10), A5SS=(1632,20), A3SS=(2917,9),  RI=(2314,13)),
}
AS_TOOL = "rMATS"

# ---- 8.SNP/Indel 合集规模 ----
VARIANTS = dict(SNP_total=14116, INDEL_total=8584)
VARIANT_CALLER = "GATK (per Novogene RNA-seq variant pipeline), annotated with SnpEff"

# ---- 8.SNP/3.stat SNP 逐样本按基因组区域分类 (真实) ----
SNP_REGION = {
    "sham_1":         dict(EXON=2309, INTRON=3523, UTR3=2142, UTR5=129, INTERGENIC=2344, UP=2000, DOWN=4343),
    "sham_2":         dict(EXON=2258, INTRON=3775, UTR3=2086, UTR5=127, INTERGENIC=2316, UP=1997, DOWN=4264),
    "sham_3":         dict(EXON=2299, INTRON=3907, UTR3=2077, UTR5=133, INTERGENIC=2350, UP=2041, DOWN=4369),
    "RFP_CCI_1":      dict(EXON=2608, INTRON=4253, UTR3=2495, UTR5=129, INTERGENIC=2506, UP=2290, DOWN=5129),
    "RFP_CCI_2":      dict(EXON=2693, INTRON=3994, UTR3=2854, UTR5=125, INTERGENIC=2459, UP=2155, DOWN=5240),
    "RFP_CCI_3":      dict(EXON=2430, INTRON=3749, UTR3=2383, UTR5=139, INTERGENIC=2386, UP=2037, DOWN=4563),
    "IL37_RFP_CCI_1": dict(EXON=2006, INTRON=3170, UTR3=2033, UTR5=104, INTERGENIC=2278, UP=1953, DOWN=4543),
    "IL37_RFP_CCI_2": dict(EXON=2272, INTRON=3713, UTR3=2197, UTR5=125, INTERGENIC=2292, UP=1922, DOWN=4401),
    "IL37_RFP_CCI_3": dict(EXON=2640, INTRON=7308, UTR3=2472, UTR5=159, INTERGENIC=2979, UP=3219, DOWN=6071),
}

# ---- SNP 功能后果 (真实) ----
SNP_FUNCTION = {
    "sham_1":         dict(MISSENSE=384, NONSENSE=2, SILENT=239),
    "sham_2":         dict(MISSENSE=357, NONSENSE=5, SILENT=228),
    "sham_3":         dict(MISSENSE=392, NONSENSE=3, SILENT=229),
    "RFP_CCI_1":      dict(MISSENSE=475, NONSENSE=5, SILENT=282),
    "RFP_CCI_2":      dict(MISSENSE=476, NONSENSE=6, SILENT=303),
    "RFP_CCI_3":      dict(MISSENSE=375, NONSENSE=3, SILENT=275),
    "IL37_RFP_CCI_1": dict(MISSENSE=255, NONSENSE=1, SILENT=138),
    "IL37_RFP_CCI_2": dict(MISSENSE=378, NONSENSE=5, SILENT=237),
    "IL37_RFP_CCI_3": dict(MISSENSE=439, NONSENSE=8, SILENT=257),
}

# ---- 3.Assemble 新转录本 (StringTie) ----
NOVEL = dict(novel_genes=1026, novel_transcripts=1396, novel_exons=3856)
ASSEMBLE_TOOL = "StringTie"
