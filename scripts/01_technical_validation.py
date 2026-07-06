"""
IL-37 CCI hippocampus bulk RNA-seq — Technical Validation 分析
数据论文 (Scientific Data). 只做数据质控, 不做机制解读.

产出:
  results/tables/Table1_sequencing_qc.csv     测序/比对质控表
  results/tables/Table2_sample_info.csv       样本信息表
  results/tables/deg_official_summary.csv      Novogene DEG 概况
  results/tables/marker_response.csv           损伤/受体 marker 表达
  results/figures/Fig1_seq_qc.png              测序质控综合图
  results/figures/Fig2_pca.png                 PCA 分组分离
  results/figures/Fig3_correlation.png         样本间相关性热图
  results/figures/Fig4_marker_response.png     marker 响应

所有输入真实, 零模拟.
"""
import os, sys
os.environ["PYTHONUTF8"] = "1"
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import novogene_qc_stats as NQ

ROOT = os.path.abspath(os.path.join(HERE, ".."))
COUNTS = r"F:\IL-37\analysis\deconvolution\results\bulk_counts.tsv"
TBL = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "figures")
os.makedirs(TBL, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

SAMPLES = NQ.SAMPLES
GROUP = NQ.GROUP
GROUP_ORDER = ["Sham", "RFP_CCI", "IL37_RFP_CCI"]
GROUP_LABEL = {"Sham": "Sham", "RFP_CCI": "CCI", "IL37_RFP_CCI": "CCI + IL-37"}

# ---------------------------------------------------------------------------
# 1. 载入 count 矩阵 + 纯度/完整性验证 (输入数据验证铁律)
# ---------------------------------------------------------------------------
def load_counts():
    df = pd.read_csv(COUNTS, sep="\t")
    assert df.columns[0] == "gene_name", f"unexpected header: {df.columns[0]}"
    # 去掉首行 '-'(featureCounts 的 unassigned 汇总行), 保留真实基因
    df = df[df["gene_name"] != "-"].copy()
    # 列顺序对齐
    missing = [s for s in SAMPLES if s not in df.columns]
    assert not missing, f"missing sample columns: {missing}"
    df = df.set_index("gene_name")[SAMPLES]
    # 去重基因名 (featureCounts 偶有重复 symbol) — 求和合并
    df = df.groupby(level=0).sum()
    # 纯度验证: 脑组织必须神经元/胶质 marker 高表达
    for mk in ["Snap25", "Gfap"]:
        assert mk in df.index, f"marker {mk} absent — wrong matrix?"
    assert df.loc["Snap25"].mean() > 1000, "Snap25 too low — not brain tissue?"
    print(f"[load] genes={df.shape[0]} samples={df.shape[1]} "
          f"total_counts_range={df.sum().min():,}-{df.sum().max():,}")
    return df

# ---------------------------------------------------------------------------
# 2. CPM / logCPM 归一化
# ---------------------------------------------------------------------------
def logcpm(counts):
    lib = counts.sum(axis=0)
    cpm = counts.divide(lib, axis=1) * 1e6
    return np.log2(cpm + 1)

# ---------------------------------------------------------------------------
# 3. 质控表 Table 1 (测序 + 比对) — 全部来自 Novogene 真实统计
# ---------------------------------------------------------------------------
def build_qc_table(counts):
    rows = []
    for s in SAMPLES:
        q, m, r, ins = NQ.QC[s], NQ.MAP[s], NQ.REGION[s], NQ.INSERT[s]
        rows.append(dict(
            Sample=s, Group=GROUP_LABEL[GROUP[s]],
            Raw_reads=q["raw"], Clean_reads=q["clean"],
            Clean_read_pairs=q["clean"],  # paired-end: reads reported per Novogene
            Q20_pct=q["q20"], Q30_pct=q["q30"], GC_pct=q["gc"], Error_pct=q["err"],
            Total_mapped_pct=m["total"], Unique_mapped_pct=m["unique"], Multi_mapped_pct=m["multi"],
            Exon_pct=round(r["exon"], 2), Intron_pct=round(r["intron"], 2),
            Intergenic_pct=round(r["intergenic"], 2),
            Insert_mean_bp=ins["mean"],
            Assigned_gene_counts=int(counts[s].sum()),
        ))
    return pd.DataFrame(rows)

# ---------------------------------------------------------------------------
# 4. 样本信息表 Table 2
# ---------------------------------------------------------------------------
def build_sample_info():
    vec = {"Sham": "None (sham surgery)",
           "RFP_CCI": "pAAV-T2A-mCherry (control)",
           "IL37_RFP_CCI": "pAAV-hIL-37-T2A-mCherry"}
    rows = []
    for i, s in enumerate(SAMPLES, 1):
        g = GROUP[s]
        rows.append(dict(
            Sample=s, Group=GROUP_LABEL[g], Biological_replicate=(i - 1) % 3 + 1,
            Species="Mus musculus (C57BL/6NJ)", Sex="Male",
            Tissue="Ipsilateral hippocampus", Injury_model="Controlled cortical impact (CCI)",
            Injection_site="Hippocampal CA3 (stereotaxic)", AAV_vector=vec[g],
            Collection_time="7 days post-CCI",
        ))
    return pd.DataFrame(rows)

# ---------------------------------------------------------------------------
# 5. marker 响应表 (损伤 marker + IL-37 受体机器) — 客观描述, 不做机制解读
# ---------------------------------------------------------------------------
def build_marker_table(counts):
    cpm = counts.divide(counts.sum(axis=0), axis=1) * 1e6
    injury = ["Gfap", "Aif1", "C1qa", "C1qb", "C3", "Cd68", "Trem2", "Tyrobp"]
    neuron = ["Snap25", "Syt1", "Rbfox3"]
    receptor = ["Il1rl2", "Sigirr", "Il18r1", "Il18rap", "Il1rap"]  # IL-37 受体机器(鼠保守)
    panels = [("Injury/glial", injury), ("Neuronal", neuron), ("IL-37 receptor", receptor)]
    rows = []
    for panel, genes in panels:
        for g in genes:
            if g not in cpm.index:
                continue
            rec = dict(Panel=panel, Gene=g)
            for grp in GROUP_ORDER:
                cols = [s for s in SAMPLES if GROUP[s] == grp]
                rec[f"{GROUP_LABEL[grp]}_meanCPM"] = round(cpm.loc[g, cols].mean(), 2)
            rows.append(rec)
    return pd.DataFrame(rows)

if __name__ == "__main__":
    counts = load_counts()
    lcpm = logcpm(counts)

    # 过滤低表达 (CPM>1 in >=3 samples) 用于 PCA/相关
    cpm = counts.divide(counts.sum(axis=0), axis=1) * 1e6
    keep = (cpm > 1).sum(axis=1) >= 3
    lcpm_f = lcpm.loc[keep]
    print(f"[filter] expressed genes kept for PCA/corr: {keep.sum()} / {len(keep)}")
    lcpm_f.to_csv(os.path.join(TBL, "logCPM_filtered.csv"))

    # ---- 表 ----
    qc = build_qc_table(counts)
    qc.to_csv(os.path.join(TBL, "Table1_sequencing_qc.csv"), index=False)
    print("[write] Table1_sequencing_qc.csv")

    si = build_sample_info()
    si.to_csv(os.path.join(TBL, "Table2_sample_info.csv"), index=False)
    print("[write] Table2_sample_info.csv")

    # Table 3 = marker 响应 (正文首次引用在 DEG 之前, 三锁编号)
    mk = build_marker_table(counts)
    mk.to_csv(os.path.join(TBL, "Table3_marker_response.csv"), index=False)
    print("[write] Table3_marker_response.csv")

    # Table 4 = DEG 概况
    deg = pd.DataFrame([dict(Comparison=k, **v) for k, v in NQ.DEG_OFFICIAL.items()])
    deg.to_csv(os.path.join(TBL, "Table4_deg_summary.csv"), index=False)
    print("[write] Table4_deg_summary.csv")

    # ---- 补充表 ----
    import supplementary as SUP
    SUP.build_supp_tables(counts, lcpm, TBL)

    # ---- 图 (主图 + 补充图) ----
    import make_figures as MF
    MF.all_figures(counts, lcpm_f, qc, mk, FIG)
    SUP.supp_figures(counts, lcpm_f, FIG)
    print("[done] tables + figures + supplementary")
