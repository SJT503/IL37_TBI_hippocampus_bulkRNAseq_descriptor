"""
Fig 1 v2: Study design schematic — 重构布局消除箭头交叉
策略:
  - 时间轴顶部横贯,三组竖排为行
  - 每行水平流动: 组标签 → [AAV注射] → CI → 取材
  - 三组汇聚到一个漏斗"↓ × 9 libraries"
  - 然后 RNA-seq → 4个产出(垂直列表,无扇形交叉)
  - panel b: 简洁样本信息表
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})

OUT = r"F:\IL-37\SciData_descriptor\results\figures\Fig1_study_design.png"

# ── 颜色 ────────────────────────────────────────────────────────────────────
C = dict(
    sham   = "#4C72B0",
    cci    = "#C44E52",
    il37   = "#55A868",
    arrow  = "#555555",
    bg     = "#F5F5F5",
    bd     = "#AAAAAA",
    output = "#7B68B0",
    seq    = "#3A3A3A",
    tl     = "#888888",
)

fig = plt.figure(figsize=(11.0, 5.0))
# 上下分区: 上方 panel a (流程图), 下方 panel b (表格)
ax = fig.add_axes([0, 0.32, 1, 0.68])   # panel a
axb = fig.add_axes([0.02, 0.01, 0.96, 0.30])  # panel b
ax.set_xlim(0, 11); ax.set_ylim(0, 4.2); ax.axis("off")
axb.axis("off")

# ── 辅助 ─────────────────────────────────────────────────────────────────────
def rbox(ax, cx, cy, w, h, fc, text, fs=7.0, tc="white", bold=False, r=0.15):
    ax.add_patch(FancyBboxPatch((cx-w/2, cy-h/2), w, h,
        boxstyle=f"round,pad=0,rounding_size={r}",
        facecolor=fc, edgecolor=C["bd"], linewidth=0.7, zorder=3))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", color=tc, zorder=4,
            multialignment="center")

def arr(ax, x1, y1, x2, y2, col=None, lw=1.0, ms=8):
    col = col or C["arrow"]
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(arrowstyle="-|>", color=col, lw=lw, mutation_scale=ms), zorder=2)

# ── 时间轴 (panel a 顶部) ────────────────────────────────────────────────────
tl_y = 3.88
ax.annotate("", xy=(8.8, tl_y), xytext=(0.65, tl_y),
    arrowprops=dict(arrowstyle="-|>", color=C["tl"], lw=1.1, mutation_scale=9), zorder=2)
ax.text(0.55, tl_y, "Time", ha="right", va="center", fontsize=6.5,
        color=C["tl"], style="italic")
for tx, tlabel in [(1.55, "−2 weeks"), (3.45, "Day 0"), (5.6, "Day 7")]:
    ax.plot([tx, tx], [tl_y, tl_y-0.14], lw=0.9, color=C["tl"])
    ax.text(tx, tl_y-0.25, tlabel, ha="center", va="top", fontsize=6.2, color=C["tl"])

# ── 三组行 ───────────────────────────────────────────────────────────────────
groups = [
    dict(y=3.0,  col=C["sham"], lbl="Sham\n(n = 3)",
         aav=None,
         surg="Sham\nsurgery",
         tissue="Ipsilateral\nhippocampus"),
    dict(y=2.08, col=C["cci"],  lbl="RFP_CCI\n(n = 3)",
         aav="pAAV-T2A-\nmCherry",
         surg="CCI",
         tissue="Ipsilateral\nhippocampus"),
    dict(y=1.16, col=C["il37"], lbl="IL37_RFP_CCI\n(n = 3)",
         aav="pAAV-hIL-37-\nT2A-mCherry",
         surg="CCI",
         tissue="Ipsilateral\nhippocampus"),
]

for g in groups:
    gy = g["y"]
    # 组标签
    rbox(ax, 1.0, gy, 1.4, 0.60, g["col"], g["lbl"], fs=7.0, bold=True)
    # AAV 注射
    if g["aav"]:
        arr(ax, 1.72, gy, 2.08, gy)
        rbox(ax, 2.68, gy, 1.12, 0.52, C["bg"], g["aav"], fs=6.2, tc="#333333", r=0.12)
        arr(ax, 3.24, gy, 3.55, gy)
    else:
        # Sham: 实线穿过 AAV 区域，文字标注"No injection"
        ax.plot([1.72, 3.55], [gy, gy], lw=0.9, color="#CCCCCC",
                linestyle="--", zorder=1)
        ax.text(2.63, gy + 0.24, "No injection", ha="center", fontsize=6.0,
                color="#AAAAAA", style="italic")
    # 手术/CI
    rbox(ax, 3.97, gy, 0.80, 0.50, C["bg"], g["surg"], fs=6.5, tc="#333333", r=0.12)
    arr(ax, 4.37, gy, 4.68, gy)
    # 取材
    rbox(ax, 5.22, gy, 1.0, 0.50, C["bg"], g["tissue"], fs=6.2, tc="#333333", r=0.12)
    # 取材框右边缘 5.72 → 竖线 5.92，实线无缺口
    ax.plot([5.71, 5.92], [gy, gy], lw=1.0, color=C["arrow"], zorder=2)

# 垂直汇聚线 x=5.92
ax.plot([5.92, 5.92], [groups[-1]["y"], groups[0]["y"]], lw=1.1, color=C["arrow"], zorder=2)
# 垂直线中点 → 汇聚框（水平箭头）
arr(ax, 5.92, groups[1]["y"], 6.30, groups[1]["y"], lw=1.0)

# ── 汇聚 → RNA-seq ────────────────────────────────────────────────────────────
mid_y = groups[1]["y"]   # 2.08
# 汇聚标注框
rbox(ax, 6.62, mid_y, 0.90, 0.48, C["bg"],
     "9 libraries\n(3 / group)", fs=6.0, tc="#555555", r=0.12)
arr(ax, 7.07, mid_y, 7.30, mid_y, lw=1.0)

# RNA-seq 框
rbox(ax, 7.88, mid_y, 1.06, 0.58, C["seq"], "Bulk RNA-seq\n(Novogene)", fs=7.0, bold=True)
arr(ax, 8.41, mid_y, 8.72, mid_y, lw=1.0)

# ── 产出列表（4条，垂直，无扇形交叉）────────────────────────────────────────
out_ys   = [3.05, 2.55, 2.08, 1.58]
out_lbls = ["Gene-level\ncount matrix",
            "AS event\ntables (rMATS)",
            "SNV / Indel\ncalls (VCF)",
            "Novel transcripts\n(GTF / FASTA)"]
out_cols = ["#5B5EA6", C["output"], C["output"], C["output"]]

# 垂直连线 x=8.72，从最上到最下产出框
ax.plot([8.72, 8.72], [out_ys[-1]-0.24, out_ys[0]+0.24], lw=1.0, color=C["output"], zorder=1)
for oy, olbl, oc in zip(out_ys, out_lbls, out_cols):
    ax.plot([8.72, 9.00], [oy, oy], lw=1.0, color=C["output"], zorder=2)
    rbox(ax, 9.76, oy, 1.52, 0.42, oc, olbl, fs=6.0, bold=False, r=0.10)

# ── panel a 标签 ─────────────────────────────────────────────────────────────
ax.text(0.05, 4.1, "a", ha="left", va="top", fontsize=11, fontweight="bold")

# ── panel b: 样本信息表 ───────────────────────────────────────────────────────
axb.set_xlim(0, 10); axb.set_ylim(0, 1.4)
axb.text(0.1, 1.32, "b  Sample information", fontsize=8.5, fontweight="bold", va="top")

hdrs = ["Group (label)", "n", "AAV vector", "Tissue", "Collection timepoint"]
hx   = [0.8, 2.0, 3.55, 6.15, 8.0]
axb.axhline(1.10, xmin=0.03, xmax=0.99, lw=0.7, color="#AAAAAA")
for hxv, hd in zip(hx, hdrs):
    axb.text(hxv, 1.16, hd, ha="left", va="center", fontsize=6.5,
             fontweight="bold", color="#333333")

rows = [
    ("Sham",         "3", "None",                        "Ipsilateral hippocampus",  "7 d post-CCI",  C["sham"]),
    ("RFP_CCI",      "3", "pAAV-T2A-mCherry",            "Ipsilateral hippocampus",  "7 d post-CCI",  C["cci"]),
    ("IL37_RFP_CCI", "3", "pAAV-hIL-37-T2A-mCherry",    "Ipsilateral hippocampus",  "7 d post-CCI",  C["il37"]),
]
ry_vals = [0.76, 0.44, 0.12]
axb.axhline(1.02, xmin=0.03, xmax=0.99, lw=0.5, color="#DDDDDD")
for ry, row in zip(ry_vals, rows):
    grp, n, aav, tis, tp, gc = row
    axb.add_patch(plt.Rectangle((0.25, ry-0.15), 9.5, 0.30,
                                 facecolor=gc, alpha=0.09, zorder=0, lw=0))
    for xv, cell in zip(hx, [grp, n, aav, tis, tp]):
        axb.text(xv, ry, cell, ha="left", va="center", fontsize=6.2, color="#333333")
axb.axhline(-0.05, xmin=0.03, xmax=0.99, lw=0.7, color="#AAAAAA")
axb.text(0.05, 1.32, "b", ha="left", va="top", fontsize=11, fontweight="bold")

fig.savefig(OUT, dpi=300)
plt.close(fig)
print(f"[done] {OUT}")
