"""
生成投稿版打包 — Scientific Data data descriptor.
1. 清理正文内部标签(核实标签/PI note)-> 投稿版 md -> pandoc docx
2. cover letter -> docx
3. 4 主表 + 3 补充表 CSV -> 独立 docx (三线表)
4. 主图/补充图 复制到投稿目录
所有内容真实, 仅做格式转换与标签清理.
"""
import os, re, shutil, subprocess
import pandas as pd

ROOT = r"F:\IL-37\SciData_descriptor"
MAN = os.path.join(ROOT, "manuscript")
TBL = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "figures")
SUB = os.path.join(ROOT, "submission")

os.makedirs(os.path.join(SUB, "01_manuscript"), exist_ok=True)
os.makedirs(os.path.join(SUB, "02_figures"), exist_ok=True)
os.makedirs(os.path.join(SUB, "03_supplementary"), exist_ok=True)
os.makedirs(os.path.join(SUB, "04_tables"), exist_ok=True)


def clean_manuscript_md():
    """剥离内部核实标签/PI note, 生成投稿版 md."""
    src = os.path.join(MAN, "data_descriptor_IL37.md")
    with open(src, encoding="utf-8") as f:
        txt = f.read()
    # 移除 *[...]*  形式的内部备注 (核实标签/PI note/待定), 跨行
    txt = re.sub(r"\*\[[^\]]*?\]\*", "", txt, flags=re.S)
    # 移除 References 段开头的核实说明块 (斜体整段)
    txt = re.sub(r"\*All references below carry.*?citation-integrity rule\.\*\s*",
                 "", txt, flags=re.S)
    # 移除行内残留的 [reused-verified]/[PubMed-verified ...] 标记
    txt = re.sub(r"\s*\*\[(reused-verified|PubMed-verified[^\]]*)\]\*", "", txt)
    # 移除作者占位行 (投稿前 PI 填真实作者)
    txt = re.sub(r"^\*\*Author list[^\n]*\(draft placeholder[^\n]*\n?", "", txt, flags=re.M)
    # 移除 cover letter 中文草稿标注注释行
    txt = re.sub(r"^>?\s*草稿[^\n]*\n?", "", txt, flags=re.M)
    # 清理因删除产生的多余空行
    txt = re.sub(r"\n{3,}", "\n\n", txt)
    out = os.path.join(SUB, "01_manuscript", "data_descriptor_IL37_submission.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write(txt)
    print(f"[clean] {os.path.basename(out)}  ({len(txt.splitlines())} lines)")
    return out


def md_to_docx(md_path, docx_path):
    subprocess.run(["pandoc", md_path, "-o", docx_path,
                    "--from", "gfm", "--standalone"], check=True)
    print(f"[docx] {os.path.basename(docx_path)}")


def csv_to_docx_table(csv_path, docx_path, title, caption=""):
    """CSV -> 三线表 docx (Sci Data 惯例)."""
    from docx import Document
    from docx.shared import Pt
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    df = pd.read_csv(csv_path)
    doc = Document()
    h = doc.add_paragraph(); r = h.add_run(title); r.bold = True; r.font.size = Pt(10)
    if caption:
        c = doc.add_paragraph(); cr = c.add_run(caption); cr.font.size = Pt(8); cr.italic = True
    t = doc.add_table(rows=1, cols=len(df.columns)); t.style = "Light List Accent 1"
    for j, col in enumerate(df.columns):
        cell = t.rows[0].cells[j]; cell.text = str(col)
        for p in cell.paragraphs:
            for rr in p.runs: rr.bold = True; rr.font.size = Pt(8)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for j, col in enumerate(df.columns):
            cells[j].text = str(row[col])
            for p in cells[j].paragraphs:
                for rr in p.runs: rr.font.size = Pt(8)
    doc.save(docx_path)
    print(f"[table] {os.path.basename(docx_path)}  ({len(df)}x{len(df.columns)})")


MAIN_TABLES = [
    ("Table1_sequencing_qc.csv", "Table 1. Sequencing and alignment quality per library.", TBL),
    ("Table2_sample_info.csv", "Table 2. Sample information.", TBL),
    ("Table3_marker_response.csv", "Table 3. Canonical marker gene expression by condition (mean CPM).", TBL),
    ("Table4_deg_summary.csv", "Table 4. Differential-expression summary (DESeq2, |log2FC|>=1 & P<=0.05).", TBL),
]
SUPP_TABLES = [
    ("supplementary/TableS1_DEG_full_lists.csv", "Supplementary Table S1. Complete per-gene differential-expression results for all three comparisons.", TBL),
    ("supplementary/TableS2_logCPM_all_genes.csv", "Supplementary Table S2. log2-CPM expression matrix for all 41,045 genes.", TBL),
    ("supplementary/TableS3_marker_CPM_per_sample.csv", "Supplementary Table S3. Per-sample CPM values for marker genes.", TBL),
    ("supplementary/TableS4_AS_events.csv", "Supplementary Table S4. Alternative-splicing event counts (rMATS) per comparison.", TBL),
    ("supplementary/TableS5_SNP_annotation.csv", "Supplementary Table S5. Per-sample SNV counts by genomic region and coding consequence.", TBL),
]


if __name__ == "__main__":
    # 1. 手稿
    md = clean_manuscript_md()
    md_to_docx(md, os.path.join(SUB, "01_manuscript", "data_descriptor_IL37_submission.docx"))
    # 2. cover letter
    cl = os.path.join(MAN, "cover_letter_IL37.md")
    if os.path.exists(cl):
        md_to_docx(cl, os.path.join(SUB, "01_manuscript", "cover_letter_IL37.docx"))
    print("[stage1] manuscript + cover letter done")

    # 3. 主表 docx
    for fn, title, base in MAIN_TABLES:
        csv_to_docx_table(os.path.join(base, fn),
                          os.path.join(SUB, "04_tables", fn.split(".")[0] + ".docx"), title)
    # 4. 补充表 docx (S1/S2 太大, 只出前若干行示例 + 完整 CSV 随附)
    for fn, title, base in SUPP_TABLES:
        src = os.path.join(base, fn)
        base_name = os.path.basename(fn).split(".")[0]
        # 完整 CSV 复制到补充目录
        shutil.copy(src, os.path.join(SUB, "03_supplementary", os.path.basename(fn)))
        # S1(4562行)/S2(41045行) 只做示例 docx(前50行), S3 全量
        df = pd.read_csv(src)
        if len(df) > 200:
            df.head(50).to_csv(src + ".preview", index=False)
            csv_to_docx_table(src + ".preview",
                              os.path.join(SUB, "03_supplementary", base_name + "_preview.docx"),
                              title + " (first 50 rows; full data in accompanying CSV)")
            os.remove(src + ".preview")
        else:
            csv_to_docx_table(src, os.path.join(SUB, "03_supplementary", base_name + ".docx"), title)

    # 5. 复制图
    for f in os.listdir(FIG):
        if f.endswith(".png"):
            shutil.copy(os.path.join(FIG, f), os.path.join(SUB, "02_figures", f))
    supp_fig = os.path.join(FIG, "supplementary")
    if os.path.isdir(supp_fig):
        for f in os.listdir(supp_fig):
            if f.endswith(".png"):
                shutil.copy(os.path.join(supp_fig, f), os.path.join(SUB, "03_supplementary", f))
    print("[stage2] tables + figures packaged")
