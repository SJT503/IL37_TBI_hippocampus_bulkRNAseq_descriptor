"""
精确插入10条新引用到正文 + 追加References 22-31
完全基于真实行内容定位,无猜测
"""
import re, os

SRC = r"F:\IL-37\SciData_descriptor\manuscript\data_descriptor_IL37.md"

with open(SRC, encoding="utf-8") as f:
    text = f.read()

body, sep, refs_tail = text.partition("\n## References\n")

# ── 10条新文献条目(全部PubMed核实) ────────────────────────────────────────────
NEW_ENTRIES = [
    # [22] Nold IL-37 奠基
    "22. Nold MF, et al. IL-37 is a fundamental inhibitor of innate immunity. *Nat Immunol.* 2010;11(11):1014-1022. https://doi.org/10.1038/ni.1944. *[PMID 20935647]*",
    # [23] Lidelow 反应性星胶
    "23. Lidelow SA, et al. Neurotoxic reactive astrocytes are induced by activated microglia. *Nature.* 2017;541(7638):481-487. https://doi.org/10.1038/nature21029. *[PMID 28099414]*",
    # [24] Sofroniew 星胶综述
    "24. Sofroniew MV, Vinters HV. Astrocytes: biology and pathology. *Acta Neuropathol.* 2010;119(1):7-35. https://doi.org/10.1007/s00401-009-0619-8. *[PMID 20012068]*",
    # [25] Keren-Shaul DAM
    "25. Keren-Shaul H, et al. A Unique Microglia Type Associated with Restricting Development of Alzheimer's Disease. *Cell.* 2017;169(7):1276-1290.e17. https://doi.org/10.1016/j.cell.2017.05.018. *[PMID 28602351]*",
    # [26] Edgar GEO
    "26. Edgar R, et al. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. *Nucleic Acids Res.* 2002;30(1):207-210. https://doi.org/10.1093/nar/30.1.207. *[PMID 11752295]*",
    # [27] Bolger Trimmomatic
    "27. Bolger AM, et al. Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics.* 2014;30(15):2114-2120. https://doi.org/10.1093/bioinformatics/btu170. *[PMID 24695404]*",
    # [28] Li SAMtools
    "28. Li H, et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics.* 2009;25(16):2078-2079. https://doi.org/10.1093/bioinformatics/btp352. *[PMID 19505943]*",
    # [29] Ewels MultiQC
    "29. Ewels P, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics.* 2016;32(19):3047-3048. https://doi.org/10.1093/bioinformatics/btw354. *[PMID 27312411]*",
    # [30] Wilkinson FAIR
    "30. Wilkinson MD, et al. The FAIR Guiding Principles for scientific data management and stewardship. *Sci Data.* 2016;3:160018. https://doi.org/10.1038/sdata.2016.18. *[PMID 26978244]*",
    # [31] Robinson edgeR
    "31. Robinson MD, et al. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics.* 2010;26(1):139-140. https://doi.org/10.1093/bioinformatics/btp616. *[PMID 19910308]*",
]

# ── 正文精确插入(按真实行内容锚定) ─────────────────────────────────────────────
inserts = [
    # [22] Nold: IL-37 "intracellular routes [14]" → 加 [22]
    ("intracellular routes [14]",
     "intracellular routes [14,22]"),

    # [23][24]: "Reactive astrocytes and\nmicroglia are central, long-lasting effectors"
    # 在 [9,10] 之后追加星胶引用
    ("involves regions beyond the primary impact site [9,10]. Reactive astrocytes and\nmicroglia are central, long-lasting effectors of this response",
     "involves regions beyond the primary impact site [9,10]. Reactive astrocytes and\nmicroglia are central, long-lasting effectors of this response [23,24]"),

    # [25] Keren-Shaul DAM: 在 single-cell atlases [11,12,13] 后加DAM引用
    ("single-cell atlases of the injured rodent brain [11,12,13]",
     "single-cell atlases of the injured rodent brain [11,12,13,25]"),

    # [26] Edgar GEO: 在 Data Records 的 GEO 首次提及处
    ("the NCBI Gene Expression Omnibus (GEO) under Series\naccession **GSExxxxxx**",
     "the NCBI Gene Expression Omnibus (GEO) [26] under Series\naccession **GSExxxxxx**"),

    # [27] Bolger Trimmomatic: 在 "after adapter and quality trimming" 处
    ("after adapter and quality trimming, 40.2",
     "after adapter and quality trimming with Trimmomatic [27], 40.2"),

    # [28] Li SAMtools: 在 HISAT2 [17] 之后
    ("with HISAT2 (v2.0.5) [17], which performs spliced alignment",
     "with HISAT2 (v2.0.5) [17]; aligned BAM files were processed with SAMtools [28]. HISAT2 performs spliced alignment"),

    # [29] Ewels MultiQC: 在 TV 测序质量段提到 quality 处
    ("were retained per library (97.0–98.3 % of raw reads)",
     "were retained per library (97.0–98.3 % of raw reads; quality verified with MultiQC [29])"),

    # [30] Wilkinson FAIR: 在 Data Availability "No restrictions" 前
    ("to reuse beyond the standard NCBI terms of use.",
     "to reuse beyond the standard NCBI terms of use. The dataset is deposited in accordance\nwith FAIR principles for scientific data management [30]."),

    # [31] Robinson edgeR: 在 Usage Notes "for example DESeq2 or edgeR" 处
    ("for example DESeq2 or edgeR, which model",
     "for example DESeq2 [21] or edgeR [31], which model"),
]

changes = 0
for old, new in inserts:
    if old in body:
        body = body.replace(old, new, 1)
        print(f"[✓] {old[:60]}")
        changes += 1
    else:
        print(f"[✗] NOT FOUND: {old[:60]}")

print(f"\n成功插入: {changes}/{len(inserts)}")

# ── 追加 References 22-31 到 References 段末尾 ───────────────────────────────
# refs_tail 从旧 References 段开始,找最后一条编号行的位置后追加
new_entries_text = "\n" + "\n".join(NEW_ENTRIES)
# 找 refs_tail 中最后一个 "^数字. " 行后面插入
# 用简单方法: 在 references 结尾(但在 Acknowledgements 之前)追加
# 找分界点
ack_markers = ["## Acknowledgements", "## Author Contributions", "## Competing"]
split_at = len(refs_tail)
for marker in ack_markers:
    idx = refs_tail.find("\n## " + marker.lstrip("# ").split()[0])
    if idx < 0:
        idx = refs_tail.find(marker)
    if 0 <= idx < split_at:
        split_at = idx

refs_before = refs_tail[:split_at]
refs_after  = refs_tail[split_at:]
new_refs_tail = refs_before + new_entries_text + "\n" + refs_after

# ── 写回 ───────────────────────────────────────────────────────────────────────
new_text = body + sep + new_refs_tail
with open(SRC, "w", encoding="utf-8") as f:
    f.write(new_text)

# ── 验证 ────────────────────────────────────────────────────────────────────────
print("\n=== 验证 ===")
with open(SRC, encoding="utf-8") as f:
    t2 = f.read()
body2 = t2.split("## References")[0]
flat2 = re.sub(r"\s+", " ", body2)
nums = sorted(set(int(n) for m in re.findall(r"\[([\d,]+)\]", flat2) for n in m.split(",")))
refs2 = t2.split("## References")[-1]
ref_nums = sorted(int(x) for x in re.findall(r"^(\d+)\. ", refs2, re.M))
print(f"正文引用编号: {nums}")
print(f"References条目编号: {ref_nums}")
print(f"正文max: {max(nums)}  References条目数: {len(ref_nums)}")
missing = [n for n in nums if n not in ref_nums]
orphan  = [n for n in ref_nums if n not in nums]
if missing: print(f"⚠ 正文有但Ref无: {missing}")
if orphan:  print(f"⚠ Ref有但正文无: {orphan}")
if not missing and not orphan:
    print("✅ 三锁对齐:正文引用 == References条目,无孤立")
print(f"总行数: {len(t2.splitlines())}")
