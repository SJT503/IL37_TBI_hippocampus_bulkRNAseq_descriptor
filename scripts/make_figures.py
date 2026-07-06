"""
Technical Validation 发表级质控图 (Scientific Data).
只画数据质控图 (测序质控/PCA/相关性/marker响应), 不搬 JNI 机制图.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})

GROUP_ORDER = ["Sham", "CCI", "CCI + IL-37"]
PALETTE = {"Sham": "#4C72B0", "CCI": "#C44E52", "CCI + IL-37": "#55A868"}
SAMPLES = ["sham_1","sham_2","sham_3","RFP_CCI_1","RFP_CCI_2","RFP_CCI_3",
           "IL37_RFP_CCI_1","IL37_RFP_CCI_2","IL37_RFP_CCI_3"]
GRP = {s:("Sham" if s.startswith("sham") else
          "CCI + IL-37" if s.startswith("IL37") else "CCI") for s in SAMPLES}
SHORT = {s: g[:-2].replace("_"," ")+"" for s,g in GRP.items()}


def _grp_colors():
    return [PALETTE[GRP[s]] for s in SAMPLES]


# ---------------------------------------------------------------------------
# Figure 1: 测序 & 比对质控
# ---------------------------------------------------------------------------
def fig_seq_qc(qc, out):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4))
    x = np.arange(len(qc))
    colors = [PALETTE[g] for g in qc["Group"]]

    # (a) Clean reads (M)
    ax = axes[0, 0]
    ax.bar(x, qc["Clean_reads"]/1e6, color=colors, edgecolor="black", linewidth=0.4)
    ax.set_ylabel("Clean reads (million pairs)")
    ax.set_title("a  Sequencing yield", loc="left", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(qc["Sample"], rotation=90)
    ax.set_ylim(0, qc["Clean_reads"].max()/1e6 * 1.15)

    # (b) Q30 / Q20
    ax = axes[0, 1]
    ax.bar(x-0.2, qc["Q30_pct"], width=0.4, color="#4C72B0", edgecolor="black",
           linewidth=0.4, label="Q30")
    ax.bar(x+0.2, qc["Q20_pct"], width=0.4, color="#9DB8D6", edgecolor="black",
           linewidth=0.4, label="Q20")
    ax.axhline(80, ls="--", lw=0.8, color="grey")
    ax.text(len(qc)-0.5, 81.5, "Illumina Q30 ≥ 80%", ha="right", va="bottom",
            fontsize=6.5, color="grey")
    # 柱实心到 96-99%; 抬高 ylim 到 108 在柱顶之上腾出留白, 图例放右上
    ax.set_ylabel("Base quality (%)"); ax.set_ylim(75, 108)
    ax.set_yticks([75, 80, 85, 90, 95, 100])
    ax.set_title("b  Base-calling quality", loc="left", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(qc["Sample"], rotation=90)
    ax.legend(loc="upper right", frameon=False, ncol=2)

    # (c) mapping rate
    ax = axes[1, 0]
    ax.bar(x-0.2, qc["Total_mapped_pct"], width=0.4, color="#55A868",
           edgecolor="black", linewidth=0.4, label="Total")
    ax.bar(x+0.2, qc["Unique_mapped_pct"], width=0.4, color="#A6D3A0",
           edgecolor="black", linewidth=0.4, label="Unique")
    # 柱最高 91.7; 抬 ylim 到 97 在柱顶之上留白, 图例放右上
    ax.set_ylabel("Reads mapped (%)"); ax.set_ylim(80, 97)
    ax.set_yticks([80, 82.5, 85, 87.5, 90, 92.5, 95])
    ax.set_title("c  Genome alignment rate", loc="left", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(qc["Sample"], rotation=90)
    ax.legend(loc="upper right", frameon=False, ncol=2)

    # (d) read distribution: exon/intron/intergenic (stacked)
    ax = axes[1, 1]
    ex, intr, inter = qc["Exon_pct"], qc["Intron_pct"], qc["Intergenic_pct"]
    ax.bar(x, ex, color="#DD8452", edgecolor="black", linewidth=0.4, label="Exon")
    ax.bar(x, intr, bottom=ex, color="#E8C4A0", edgecolor="black", linewidth=0.4,
           label="Intron")
    ax.bar(x, inter, bottom=ex+intr, color="#D9D9D9", edgecolor="black",
           linewidth=0.4, label="Intergenic")
    # 堆叠柱填满 100%; 抬 ylim 到 130 在柱顶之上留白, yticks 仍标到 100, 图例放右上
    ax.set_ylabel("Mapped reads (%)"); ax.set_ylim(0, 130)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_title("d  Genomic origin of reads", loc="left", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(qc["Sample"], rotation=90)
    ax.legend(loc="upper right", frameon=False, ncol=3, fontsize=6.5,
              handlelength=1.2, columnspacing=1.0)

    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[fig] {os.path.basename(out)}")


# ---------------------------------------------------------------------------
# Figure 2: PCA (从 logCPM 自算, 组间分离)
# ---------------------------------------------------------------------------
def fig_pca(lcpm_f, out):
    from sklearn.decomposition import PCA
    X = lcpm_f[SAMPLES].T.values
    X = X - X.mean(axis=0, keepdims=True)
    pca = PCA(n_components=3)
    sc = pca.fit_transform(X)
    var = pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    for g in GROUP_ORDER:
        idx = [i for i, s in enumerate(SAMPLES) if GRP[s] == g]
        ax.scatter(sc[idx, 0], sc[idx, 1], s=70, color=PALETTE[g],
                   edgecolor="black", linewidth=0.6, label=g, zorder=3)
    for i, s in enumerate(SAMPLES):
        ax.annotate(str((i % 3) + 1), (sc[i, 0], sc[i, 1]), fontsize=6,
                    ha="center", va="center", color="white", zorder=4, fontweight="bold")
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.set_title("Principal component analysis", loc="left", fontweight="bold")
    ax.axhline(0, ls=":", lw=0.5, color="grey", zorder=0)
    ax.axvline(0, ls=":", lw=0.5, color="grey", zorder=0)
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[fig] {os.path.basename(out)}  var={var[:3].round(1)}")
    return var


# ---------------------------------------------------------------------------
# Figure 3: 样本间相关性热图 (Spearman on logCPM)
# ---------------------------------------------------------------------------
def fig_corr(lcpm_f, out):
    corr = lcpm_f[SAMPLES].corr(method="spearman")
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(corr.values, cmap="RdYlBu_r", vmin=corr.values.min(), vmax=1.0)
    ax.set_xticks(range(len(SAMPLES))); ax.set_yticks(range(len(SAMPLES)))
    ax.set_xticklabels(SAMPLES, rotation=90); ax.set_yticklabels(SAMPLES)
    # 组色条 (左侧 + 顶部)
    for i, s in enumerate(SAMPLES):
        ax.add_patch(plt.Rectangle((-1.4, i-0.5), 0.8, 1, color=PALETTE[GRP[s]],
                     clip_on=False, lw=0))
        ax.add_patch(plt.Rectangle((i-0.5, -1.4), 1, 0.8, color=PALETTE[GRP[s]],
                     clip_on=False, lw=0))
    # 数值标注
    for i in range(len(SAMPLES)):
        for j in range(len(SAMPLES)):
            v = corr.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.2,
                    color="black" if v > 0.6 else "white")
    ax.set_xlim(-1.6, len(SAMPLES)-0.5); ax.set_ylim(len(SAMPLES)-0.5, -1.6)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Spearman correlation", fontsize=7)
    ax.set_title("Sample-to-sample correlation", loc="left", fontweight="bold")
    handles = [Patch(facecolor=PALETTE[g], label=g) for g in GROUP_ORDER]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.28, 1.0),
              frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[fig] {os.path.basename(out)}  intra_min={_intra_min(corr):.3f}")


def _intra_min(corr):
    vals = []
    for g in GROUP_ORDER:
        cols = [s for s in SAMPLES if GRP[s] == g]
        sub = corr.loc[cols, cols].values
        vals.append(sub[np.triu_indices(len(cols), 1)].min())
    return min(vals)


# ---------------------------------------------------------------------------
# Figure 4: marker 响应 (损伤 + 神经元 + IL-37 受体), 客观展示 CPM
# ---------------------------------------------------------------------------
def fig_marker(counts, mk, out):
    cpm = counts.divide(counts.sum(axis=0), axis=1) * 1e6
    panels = [("Injury / glial markers", ["Gfap", "C1qa", "C3", "Cd68", "Trem2", "Aif1"]),
              ("Neuronal markers", ["Snap25", "Syt1", "Rbfox3"]),
              ("IL-37 receptor machinery", ["Sigirr", "Il1rl2", "Il18r1", "Il18rap"])]
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.2),
                             gridspec_kw={"width_ratios": [6, 3, 4]})
    letters = "abc"
    for k, (title, genes) in enumerate(panels):
        ax = axes[k]
        genes = [g for g in genes if g in cpm.index]
        xb = np.arange(len(genes)); w = 0.26
        for gi, grp in enumerate(GROUP_ORDER):
            cols = [s for s in SAMPLES if GRP[s] == grp]
            means = np.array([cpm.loc[g, cols].mean() for g in genes])
            sds = np.array([cpm.loc[g, cols].std(ddof=1) for g in genes])
            # log 尺度下误差棒下界截断到 >0
            lower = np.minimum(sds, means * 0.999)
            ax.bar(xb + (gi-1)*w, means, w, yerr=[lower, sds], capsize=1.5,
                   color=PALETTE[grp], edgecolor="black", linewidth=0.4,
                   error_kw=dict(lw=0.5), label=grp)
        ax.set_yscale("log")
        # 顶部乘性留白 (log 尺度) 给右上图例腾空间
        ymax = max((cpm.loc[g, [s for s in SAMPLES if GRP[s]==grp]].mean()
                    for g in genes for grp in GROUP_ORDER))
        ax.set_ylim(0.1, ymax * (30 if k == 0 else 5))
        ax.set_xticks(xb); ax.set_xticklabels(genes, rotation=45, ha="right", style="italic")
        ax.set_ylabel("Mean CPM (log scale)")
        ax.set_title(f"{letters[k]}  {title}", loc="left", fontweight="bold", fontsize=8)
        if k == 0:
            ax.legend(loc="upper right", frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[fig] {os.path.basename(out)}")


def all_figures(counts, lcpm_f, qc, mk, figdir):
    fig_seq_qc(qc, os.path.join(figdir, "Fig1_seq_qc.png"))
    fig_pca(lcpm_f, os.path.join(figdir, "Fig2_pca.png"))
    fig_corr(lcpm_f, os.path.join(figdir, "Fig3_correlation.png"))
    fig_marker(counts, mk, os.path.join(figdir, "Fig4_marker_response.png"))
