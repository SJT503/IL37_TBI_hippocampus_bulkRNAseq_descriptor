"""
补充材料 (Supplementary) 生成 — Scientific Data data descriptor.
补充表 S1-S3 + 补充图 S1-S2. 全部真实数据, 与手稿 Table/正文同源.

S1: 官方 DESeq2 全量 DEG 列表 (三组比较, 与 Table 4 计数同源 = 诺禾权威输出)
S2: 全基因 logCPM 表达矩阵
S3: marker 逐样本 CPM 完整表 (Table 3 的展开)
Fig S1: 每样本基因检测数 + 文库复杂度分布
Fig S2: DEG 火山图 (三组比较, 客观呈现不做机制解读)
"""
import os, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 诺禾官方 DEG 目录 (与 Table 4 同源)
NOVO_DEG = (r"F:\IL-37\转录组测序(2)\03.Result_X101SC25026107-Z01-J003_mus_musculus_c57bl6nj"
            r"\Result_X101SC25026107-Z01-J003_mus_musculus_c57bl6nj\5.Differential\1.deglist")
COMPARISONS = {
    "RFP_CCIvsSham": "CCI_vs_Sham",
    "IL37_RFP_CCIvsRFP_CCI": "IL37_vs_CCI",
    "IL37_RFP_CCIvsSham": "IL37_vs_Sham",
}
SAMPLES = ["sham_1","sham_2","sham_3","RFP_CCI_1","RFP_CCI_2","RFP_CCI_3",
           "IL37_RFP_CCI_1","IL37_RFP_CCI_2","IL37_RFP_CCI_3"]
GRP = {s:("Sham" if s.startswith("sham") else
          "CCI + IL-37" if s.startswith("IL37") else "CCI") for s in SAMPLES}
PALETTE = {"Sham": "#4C72B0", "CCI": "#C44E52", "CCI + IL-37": "#55A868"}


def _read_novo_deg_all(comp_dir):
    """读诺禾官方全量 DEG (deg_all.xls), 返回精简列 DataFrame."""
    f = os.path.join(NOVO_DEG, comp_dir, f"{comp_dir}_deg_all.xls")
    if not os.path.exists(f):
        return None
    df = pd.read_csv(f, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    keep = {}
    for want in ["gene_name", "gene_id", "log2foldchange", "pvalue", "padj", "gene_biotype"]:
        if want in cols:
            keep[cols[want]] = want
    out = df[list(keep.keys())].rename(columns=keep)
    return out


def build_supp_tables(counts, lcpm, tbl_dir):
    supp = os.path.join(tbl_dir, "supplementary")
    os.makedirs(supp, exist_ok=True)

    # ---- S1: 官方全量 DEG 三组比较 ----
    frames = []
    for comp_dir, label in COMPARISONS.items():
        d = _read_novo_deg_all(comp_dir)
        if d is None:
            print(f"[S1][WARN] missing {comp_dir}")
            continue
        d.insert(0, "Comparison", label)
        # 显著性标记 (与官方 Table 4 同阈值: pvalue<=0.05 & |log2FC|>=1)
        d["Significant"] = ((d["pvalue"] <= 0.05) & (d["log2foldchange"].abs() >= 1.0))
        frames.append(d)
    if frames:
        s1 = pd.concat(frames, ignore_index=True)
        s1.to_csv(os.path.join(supp, "TableS1_DEG_full_lists.csv"), index=False)
        # 一致性自检: 显著 DEG 数应吻合 Table 4
        chk = s1[s1.Significant].groupby("Comparison").size().to_dict()
        print(f"[S1] TableS1_DEG_full_lists.csv  sig_counts={chk}")

    # ---- S2: 全基因 logCPM 矩阵 ----
    lcpm.round(4).to_csv(os.path.join(supp, "TableS2_logCPM_all_genes.csv"))
    print(f"[S2] TableS2_logCPM_all_genes.csv  {lcpm.shape[0]} genes")

    # ---- S3: marker 逐样本 CPM 完整表 ----
    cpm = counts.divide(counts.sum(axis=0), axis=1) * 1e6
    injury = ["Gfap","Aif1","C1qa","C1qb","C3","Cd68","Trem2","Tyrobp","Il1b","Tnf"]
    neuron = ["Snap25","Syt1","Rbfox3","Tubb3"]
    receptor = ["Il1rl2","Sigirr","Il18r1","Il18rap","Il1rap"]
    panel_map = ([("Injury/glial", g) for g in injury] +
                 [("Neuronal", g) for g in neuron] +
                 [("IL-37 receptor", g) for g in receptor])
    rows = []
    for panel, g in panel_map:
        if g not in cpm.index:
            continue
        rec = {"Panel": panel, "Gene": g}
        for s in SAMPLES:
            rec[s] = round(cpm.loc[g, s], 3)
        rows.append(rec)
    pd.DataFrame(rows).to_csv(os.path.join(supp, "TableS3_marker_CPM_per_sample.csv"), index=False)
    print(f"[S3] TableS3_marker_CPM_per_sample.csv  {len(rows)} markers")

    # ---- S4: 可变剪接事件数 (rMATS 5类, 全量+显著) ----
    import extended_stats as EX
    rows = []
    for comp, ev in EX.AS_EVENTS.items():
        rec = {"Comparison": comp}
        tot = 0
        for typ in ["SE", "MXE", "A5SS", "A3SS", "RI"]:
            n, sig = ev[typ]
            rec[f"{typ}_total"] = n
            rec[f"{typ}_significant"] = sig
            tot += n
        rec["All_events_total"] = tot
        rows.append(rec)
    pd.DataFrame(rows).to_csv(os.path.join(supp, "TableS4_AS_events.csv"), index=False)
    print(f"[S4] TableS4_AS_events.csv  ({EX.AS_TOOL})")

    # ---- S5: SNP/Indel 逐样本按区域 + 功能 分类 ----
    rows = []
    for s in SAMPLES:
        r = EX.SNP_REGION[s]; fn = EX.SNP_FUNCTION[s]
        rows.append(dict(Sample=s, **{f"region_{k}": v for k, v in r.items()},
                         **{f"function_{k}": v for k, v in fn.items()}))
    pd.DataFrame(rows).to_csv(os.path.join(supp, "TableS5_SNP_annotation.csv"), index=False)
    print(f"[S5] TableS5_SNP_annotation.csv  (total {EX.VARIANTS['SNP_total']} SNP / "
          f"{EX.VARIANTS['INDEL_total']} Indel)")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9,
    "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})


def supp_figures(counts, lcpm_f, figdir):
    supp = os.path.join(figdir, "supplementary")
    os.makedirs(supp, exist_ok=True)
    _figS1_detection(counts, os.path.join(supp, "FigS1_gene_detection.png"))
    _figS2_volcano(os.path.join(supp, "FigS2_volcano.png"))
    _figS3_as_events(os.path.join(supp, "FigS3_AS_events.png"))
    _figS4_variants(os.path.join(supp, "FigS4_variant_annotation.png"))


# ---- Fig S3: 可变剪接事件类型分布 (rMATS 全量, 客观展示) ----
def _figS3_as_events(out):
    import extended_stats as EX
    types = ["SE", "MXE", "A5SS", "A3SS", "RI"]
    comps = list(EX.AS_EVENTS.keys())
    comp_lbl = {"IL37_vs_CCI": "CCI + IL-37\nvs CCI", "IL37_vs_Sham": "CCI + IL-37\nvs Sham",
                "CCI_vs_Sham": "CCI vs Sham"}
    tcolors = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3"]
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4), gridspec_kw={"width_ratios": [3, 3]})

    # (a) 全量事件数 (堆叠, 各比较接近 -> 展示数据集规模)
    ax = axes[0]
    x = np.arange(len(comps)); bottom = np.zeros(len(comps))
    for ti, typ in enumerate(types):
        vals = [EX.AS_EVENTS[c][typ][0] for c in comps]
        ax.bar(x, vals, bottom=bottom, color=tcolors[ti], edgecolor="black",
               linewidth=0.3, label=typ)
        bottom += np.array(vals)
    ax.set_xticks(x); ax.set_xticklabels([comp_lbl[c] for c in comps], fontsize=6.5)
    ax.set_ylabel("Alternative splicing events (total)")
    ax.set_title("a  AS events detected (rMATS)", loc="left", fontweight="bold")
    # 堆叠柱最高~33k; 抬 ylim 顶部留白, 图例放右上
    ax.set_ylim(0, bottom.max() * 1.32)
    ax.legend(loc="upper right", frameon=False, fontsize=6.5, ncol=2)

    # (b) 显著差异事件数 (log 尺度, 数量小)
    ax = axes[1]
    w = 0.15
    maxv = max(EX.AS_EVENTS[c][t][1] for c in comps for t in types)
    for ti, typ in enumerate(types):
        vals = [EX.AS_EVENTS[c][typ][1] for c in comps]
        ax.bar(x + (ti - 2) * w, vals, w, color=tcolors[ti], edgecolor="black",
               linewidth=0.3, label=typ)
    ax.set_xticks(x); ax.set_xticklabels([comp_lbl[c] for c in comps], fontsize=6.5)
    ax.set_ylabel("Significant differential AS events (FDR < 0.05)")
    ax.set_title("b  Significant differential AS events", loc="left", fontweight="bold")
    ax.set_ylim(0, maxv * 1.30)  # 顶部留白给图例
    ax.legend(loc="upper right", frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[figS3] {os.path.basename(out)}")


# ---- Fig S4: SNP 基因组区域分布 + 功能后果 (客观展示变异数据特征) ----
def _figS4_variants(out):
    import extended_stats as EX
    regions = ["EXON", "INTRON", "UTR3", "UTR5", "INTERGENIC", "UP", "DOWN"]
    rlabel = {"EXON": "Exon", "INTRON": "Intron", "UTR3": "3'UTR", "UTR5": "5'UTR",
              "INTERGENIC": "Intergenic", "UP": "Upstream", "DOWN": "Downstream"}
    colors = [PALETTE[GRP[s]] for s in SAMPLES]
    x = np.arange(len(SAMPLES))
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4), gridspec_kw={"width_ratios": [4, 3]})

    # (a) SNP 区域分布 (堆叠)
    ax = axes[0]
    rcolors = ["#DD8452", "#E8C4A0", "#55A868", "#A6D3A0", "#C44E52", "#B0B0B0", "#D9D9D9"]
    bottom = np.zeros(len(SAMPLES))
    for ri, reg in enumerate(regions):
        vals = [EX.SNP_REGION[s][reg] for s in SAMPLES]
        ax.bar(x, vals, bottom=bottom, color=rcolors[ri], edgecolor="black",
               linewidth=0.25, label=rlabel[reg])
        bottom += np.array(vals)
    ax.set_xticks(x); ax.set_xticklabels(SAMPLES, rotation=90)
    ax.set_ylabel("SNPs by genomic region")
    ax.set_title("a  SNP distribution by region", loc="left", fontweight="bold")
    # 堆叠柱; 抬 ylim 顶部留白, 图例放右上 (7类用2列紧凑)
    ax.set_ylim(0, bottom.max() * 1.42)
    ax.legend(loc="upper right", frameon=False, fontsize=6, ncol=2, columnspacing=1.0)

    # (b) SNP 功能后果 (missense/nonsense/silent)
    ax = axes[1]
    fcats = ["MISSENSE", "NONSENSE", "SILENT"]
    fcolors = ["#C44E52", "#8172B3", "#4C72B0"]
    bottom = np.zeros(len(SAMPLES))
    for fi, fc in enumerate(fcats):
        vals = [EX.SNP_FUNCTION[s][fc] for s in SAMPLES]
        ax.bar(x, vals, bottom=bottom, color=fcolors[fi], edgecolor="black",
               linewidth=0.25, label=fc.capitalize())
        bottom += np.array(vals)
    ax.set_xticks(x); ax.set_xticklabels(SAMPLES, rotation=90)
    ax.set_ylabel("Coding SNPs by functional consequence")
    ax.set_title("b  Coding SNP consequences", loc="left", fontweight="bold")
    # 堆叠柱; 抬 ylim 顶部留白, 图例放右上
    ax.set_ylim(0, bottom.max() * 1.32)
    ax.legend(loc="upper right", frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[figS4] {os.path.basename(out)}")


# ---- Fig S1: 每样本基因检测数 + 表达分布 (测序充分性) ----
def _figS1_detection(counts, out):
    colors = [PALETTE[GRP[s]] for s in SAMPLES]
    x = np.arange(len(SAMPLES))
    # 检测基因数 (count>0)
    n_genes = [(counts[s] > 0).sum() for s in SAMPLES]
    # 每样本文库大小 (总 assigned counts, 百万)
    libsize = [counts[s].sum() / 1e6 for s in SAMPLES]

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.2))
    ax = axes[0]
    ax.bar(x, n_genes, color=colors, edgecolor="black", linewidth=0.4)
    ax.set_ylabel("Genes detected (count > 0)")
    ax.set_title("a  Gene detection per sample", loc="left", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(SAMPLES, rotation=90)
    ax.set_ylim(0, max(n_genes) * 1.12)

    ax = axes[1]
    ax.bar(x, libsize, color=colors, edgecolor="black", linewidth=0.4)
    ax.set_ylabel("Assigned gene counts (million)")
    ax.set_title("b  Library size per sample", loc="left", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(SAMPLES, rotation=90)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[figS1] {os.path.basename(out)}  genes={min(n_genes)}-{max(n_genes)}")


# ---- Fig S2: 官方 DEG 火山图 三联 (客观呈现, 不做机制解读) ----
def _figS2_volcano(out):
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.4))
    titles = {"RFP_CCIvsSham": "CCI vs Sham",
              "IL37_RFP_CCIvsRFP_CCI": "CCI + IL-37 vs CCI",
              "IL37_RFP_CCIvsSham": "CCI + IL-37 vs Sham"}
    letters = "abc"
    for k, (comp_dir, title) in enumerate(titles.items()):
        ax = axes[k]
        f = os.path.join(NOVO_DEG, comp_dir, f"{comp_dir}_deg_all.xls")
        if not os.path.exists(f):
            ax.text(0.5, 0.5, "n/a", ha="center"); continue
        df = pd.read_csv(f, sep="\t")
        cols = {c.lower(): c for c in df.columns}
        lfc = pd.to_numeric(df[cols["log2foldchange"]], errors="coerce")
        pv = pd.to_numeric(df[cols["pvalue"]], errors="coerce")
        ok = lfc.notna() & pv.notna() & (pv > 0)
        lfc, pv = lfc[ok], pv[ok]
        neglog = -np.log10(pv)
        sig = (pv <= 0.05) & (lfc.abs() >= 1.0)
        up = sig & (lfc > 0); down = sig & (lfc < 0)
        ax.scatter(lfc[~sig], neglog[~sig], s=3, c="#B0B0B0", alpha=0.5, linewidths=0)
        ax.scatter(lfc[down], neglog[down], s=4, c="#4C72B0", alpha=0.7, linewidths=0, label=f"down ({down.sum()})")
        ax.scatter(lfc[up], neglog[up], s=4, c="#C44E52", alpha=0.7, linewidths=0, label=f"up ({up.sum()})")
        ax.axvline(-1, ls="--", lw=0.5, color="grey"); ax.axvline(1, ls="--", lw=0.5, color="grey")
        ax.axhline(-np.log10(0.05), ls="--", lw=0.5, color="grey")
        ax.set_xlabel("log2 fold change"); ax.set_ylabel("-log10 P" if k == 0 else "")
        ax.set_title(f"{letters[k]}  {title}", loc="left", fontweight="bold", fontsize=8)
        ax.legend(loc="upper right", frameon=False, fontsize=6, markerscale=2)
        lim = np.nanpercentile(lfc.abs(), 99.5)
        ax.set_xlim(-lim, lim)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[figS2] {os.path.basename(out)}")
