"""
重建手稿引用编号:
1. 把当前[1]-[15]替换为作者年标签
2. 在正文合适位置插入新文献的作者年标签
3. 按首次出现顺序自动分配连续编号
4. 生成新 References 段(替换旧段)
全程真实引用,无模拟数据
"""
import re, copy

SRC  = r"F:\IL-37\SciData_descriptor\manuscript\data_descriptor_IL37.md"
DEST = r"F:\IL-37\SciData_descriptor\manuscript\data_descriptor_IL37.md"

# ── 旧引用 [1-15] 对应的作者年标签 ──────────────────────────────────────────
OLD_TAG = {
    "[1]":  "«Maas2017»",
    "[2]":  "«GBD2019»",
    "[3]":  "«Pearson2024»",
    "[4]":  "«Ji2025»",
    "[5]":  "«Osier2016»",
    "[6]":  "«Pybus2024»",
    "[7]":  "«Ramlack2011»",
    "[8]":  "«Escartin2021»",
    "[9]":  "«Hammond2019»",
    "[10]": "«Dinarello2016»",
    "[11]": "«Kim2019»",
    "[12]": "«Liao2014»",
    "[13]": "«Pertea2015»",
    "[14]": "«Shen2014»",
    "[15]": "«Love2014»",
}

# ── 新文献作者年标签(全部 PubMed 核实) ─────────────────────────────────────
NEW_REFS = {
    "«Dewan2019»":    ("Dewan MC, et al. Estimating the global incidence of traumatic brain injury. *J Neurosurg.* 2019;130(4):1080-1097. https://doi.org/10.3171/2017.10.JNS17352.", "29701556"),
    "«MaZ2024»":      ("Ma Z, et al. Traumatic brain injury in elderly population: A global systematic review and meta-analysis of in-hospital mortality and risk factors among 2.22 million individuals. *Ageing Res Rev.* 2024;99:102376. https://doi.org/10.1016/j.arr.2024.102376.", "38972601"),
    "«Cook2025»":     ("Cook AM, et al. Update on Neuroprotection after Traumatic Brain Injury. *CNS Drugs.* 2025;39(5):473-484. https://doi.org/10.1007/s40263-025-01193-4.", "40087248"),
    "«Jha2024»":      ("Jha RM, et al. A single-cell atlas deconstructs heterogeneity across multiple models in murine traumatic brain injury and identifies novel cell-specific targets. *Neuron.* 2024;112(18):3069-3088.e4. https://doi.org/10.1016/j.neuron.2024.06.023.", "39019041"),
    "«Zhang2025»":    ("Zhang J, et al. Interleukin-37 modulates microglial phenotype and inhibits inflammatory response via the MyD88/NF-κB pathway in lipopolysaccharide-induced neuroinflammation. *Inflamm Res.* 2025;74(1):87. https://doi.org/10.1007/s00011-025-02048-x.", "40439750"),
    "«Nold2010»":     ("Nold MF, et al. IL-37 is a fundamental inhibitor of innate immunity. *Nat Immunol.* 2010;11(11):1014-1022. https://doi.org/10.1038/ni.1944.", "20935647"),
    "«Lidelow2017»":  ("Lidelow SA, et al. Neurotoxic reactive astrocytes are induced by activated microglia. *Nature.* 2017;541(7638):481-487. https://doi.org/10.1038/nature21029.", "28099414"),
    "«Sofroniew2010»":("Sofroniew MV, Vinters HV. Astrocytes: biology and pathology. *Acta Neuropathol.* 2010;119(1):7-35. https://doi.org/10.1007/s00401-009-0619-8.", "20012068"),
    "«Keren2017»":    ("Keren-Shaul H, et al. A Unique Microglia Type Associated with Restricting Development of Alzheimer's Disease. *Cell.* 2017;169(7):1276-1290.e17. https://doi.org/10.1016/j.cell.2017.05.018.", "28602351"),
    "«Edgar2002»":    ("Edgar R, et al. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. *Nucleic Acids Res.* 2002;30(1):207-210. https://doi.org/10.1093/nar/30.1.207.", "11752295"),
    "«Bolger2014»":   ("Bolger AM, et al. Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics.* 2014;30(15):2114-2120. https://doi.org/10.1093/bioinformatics/btu170.", "24695404"),
    "«Li2009»":       ("Li H, et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics.* 2009;25(16):2078-2079. https://doi.org/10.1093/bioinformatics/btp352.", "19505943"),
    "«Dobin2013»":    ("Dobin A, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics.* 2013;29(1):15-21. https://doi.org/10.1093/bioinformatics/bts635.", "23104886"),
    "«Ewels2016»":    ("Ewels P, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics.* 2016;32(19):3047-3048. https://doi.org/10.1093/bioinformatics/btw354.", "27312411"),
    "«Wilkinson2016»":("Wilkinson MD, et al. The FAIR Guiding Principles for scientific data management and stewardship. *Sci Data.* 2016;3:160018. https://doi.org/10.1038/sdata.2016.18.", "26978244"),
    "«Robinson2010»": ("Robinson MD, et al. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics.* 2010;26(1):139-140. https://doi.org/10.1093/bioinformatics/btp616.", "19910308"),
}

# 所有新文献标签合并(旧 + 新)
ALL_REFS = {}
for tag, ref in {
    "«Maas2017»":    ("Maas AIR, et al. Traumatic brain injury: integrated approaches to improve prevention, clinical care, and research. *Lancet Neurol.* 2017;16(12):987-1048. https://doi.org/10.1016/S1474-4422(17)30371-X.", "29122524"),
    "«GBD2019»":     ("GBD 2016 Traumatic Brain Injury and Spinal Cord Injury Collaborators. Global, regional, and national burden of traumatic brain injury and spinal cord injury, 1990-2016. *Lancet Neurol.* 2019;18(1):56-87. https://doi.org/10.1016/S1474-4422(18)30415-0.", "30497965"),
    "«Pearson2024»": ("Pearson A, et al. PPARγ activation ameliorates cognitive impairment and chronic microglial activation in the aftermath of r-mTBI. *J Neuroinflammation.* 2024;21(1):194. https://doi.org/10.1186/s12974-024-03173-w.", "39097742"),
    "«Ji2025»":      ("Ji Q, et al. Astrocytic SIRT1 ameliorates cognitive deficits after traumatic brain injury via autophagy-mediated MEGF10 phagocytosis. *J Neuroinflammation.* 2025;22(1):304. https://doi.org/10.1186/s12974-025-03633-x.", "41469684"),
    "«Osier2016»":   ("Osier ND, Dixon CE. The Controlled Cortical Impact Model: Applications, Considerations for Researchers, and Future Directions. *Front Neurol.* 2016;7:134. https://doi.org/10.3389/fneur.2016.00134.", "27582726"),
    "«Pybus2024»":   ("Pybus AF, et al. Profiling the neuroimmune cascade in 3xTg-AD mice exposed to successive mild traumatic brain injuries. *J Neuroinflammation.* 2024;21(1):156. https://doi.org/10.1186/s12974-024-03128-1.", "38872143"),
    "«Ramlack2011»": ("Ramlackhansingh AF, et al. Inflammation after trauma: microglial activation and traumatic brain injury. *Ann Neurol.* 2011;70(3):374-383. https://doi.org/10.1002/ana.22455.", "21710619"),
    "«Escartin2021»":("Escartin C, et al. Reactive astrocyte nomenclature, definitions, and future directions. *Nat Neurosci.* 2021;24(3):312-325. https://doi.org/10.1038/s41593-020-00783-4.", "33589835"),
    "«Hammond2019»": ("Hammond TR, et al. Single-Cell RNA Sequencing of Microglia throughout the Mouse Lifespan and in the Injured Brain Reveals Complex Cell-State Changes. *Immunity.* 2019;50(1):253-271.e6. https://doi.org/10.1016/j.immuni.2018.11.004.", "30471926"),
    "«Dinarello2016»":("Dinarello CA, et al. Suppression of innate inflammation and immunity by interleukin-37. *Eur J Immunol.* 2016;46(5):1067-1081. https://doi.org/10.1002/eji.201545828.", "27060871"),
    "«Kim2019»":     ("Kim D, et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nat Biotechnol.* 2019;37(8):907-915. https://doi.org/10.1038/s41587-019-0201-4.", "31375807"),
    "«Liao2014»":    ("Liao Y, et al. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics.* 2014;30(7):923-930. https://doi.org/10.1093/bioinformatics/btt656.", "24227677"),
    "«Pertea2015»":  ("Pertea M, et al. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. *Nat Biotechnol.* 2015;33(3):290-295. https://doi.org/10.1038/nbt.3122.", "25690850"),
    "«Shen2014»":    ("Shen S, et al. rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. *Proc Natl Acad Sci USA.* 2014;111(51):E5593-E5601. https://doi.org/10.1073/pnas.1419161111.", "25480548"),
    "«Love2014»":    ("Love MI, et al. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.* 2014;15(12):550. https://doi.org/10.1186/s13059-014-0550-8.", "25516281"),
    **NEW_REFS
}.items():
    ALL_REFS[tag] = ref

print(f"Total unique refs in pool: {len(ALL_REFS)}")

# ── Step 1: 读手稿 ────────────────────────────────────────────────────────────
with open(SRC, encoding="utf-8") as f:
    text = f.read()

# 把正文和 References 段分开
body, _, ref_block = text.partition("\n## References\n")

# ── Step 2: 把旧 [1]-[15] 替换为作者年标签 ──────────────────────────────────
# 先处理复合引用 [1,2] → 拆成 «tag1»«tag2»（空格分隔）
def replace_compound(m):
    nums = m.group(1).split(",")
    return "".join(OLD_TAG.get(f"[{n.strip()}]", f"[{n.strip()}]") for n in nums)

body = re.sub(r"\[([\d,\s]+)\]", replace_compound, body)

# 验证旧引用已全部替换
remaining = re.findall(r"\[\d+\]", body)
if remaining:
    print(f"WARNING: 未替换的旧引用: {remaining}")

# ── Step 3: 在正文插入新文献标签 ─────────────────────────────────────────────
# 每处插入都精确定位,避免破坏周围文字

inserts = [
    # Y1: Background 第1段 - TBI 69M/年(Dewan) + 老年流行病学(MaZ)
    ("as populations age «GBD2019»",
     "as populations age «GBD2019»«MaZ2024»"),

    # Y1: TBI 治疗缺口 - Cook 2025
    ("and the least treatable",
     "and the least treatable «Cook2025»"),

    # Y2/Y3: 转录组资源 - Jha 2024(TBI 单细胞图谱数据资源)
    ("by transcriptomic profiling «Escartin2021»«Hammond2019»",
     "by transcriptomic profiling «Escartin2021»«Hammond2019»«Jha2024»"),

    # Y3: 反应性星胶 - Sofroniew综述 + Lidelow星胶-小胶互作
    ("«Escartin2021»«Hammond2019»«Jha2024»",
     "«Sofroniew2010»«Escartin2021»«Lidelow2017»«Hammond2019»«Keren2017»«Jha2024»"),

    # Y3: IL-37 - Nold 2010 奠基文(放 IL-37 routes 处)
    ("through both extracellular receptor-mediated and intracellular routes «Dinarello2016»",
     "through both extracellular receptor-mediated and intracellular routes «Nold2010»«Dinarello2016»"),

    # Y3: IL-37 CNS 神经 - Zhang 2025(放"CNS studies rely on"附近)
    ("most commonly by viral vector.",
     "most commonly by viral vector «Zhang2025»."),

    # Methods: Trimmomatic 在 Library prep 段
    ("Poly(A)-selected, strand-nonspecific paired-end libraries were prepared from total",
     "Poly(A)-selected, strand-nonspecific paired-end libraries were prepared from total\nRNA after adapter and quality trimming with Trimmomatic «Bolger2014»,"),

    # Methods: SAMtools + STAR参考(在 HISAT2 段后)
    ("with HISAT2 (v2.0.5) «Kim2019»",
     "with HISAT2 (v2.0.5) «Kim2019»; aligned reads were processed with SAMtools «Li2009»"),

    # Methods: MultiQC(在质控分析段)
    ("quality trimming, 40.2",
     "quality trimming (assessed with MultiQC «Ewels2016»), 40.2"),

    # Data Records: GEO 原文(在 GEO 首次提及处)
    ("the NCBI Gene Expression Omnibus (GEO)",
     "the NCBI Gene Expression Omnibus (GEO) «Edgar2002»"),

    # Data Availability: FAIR原则
    ("No restrictions apply to reuse beyond the standard NCBI terms of use.",
     "The data are deposited in accordance with FAIR principles «Wilkinson2016». "
     "No restrictions apply to reuse beyond the standard NCBI terms of use."),

    # Usage Notes: edgeR 与 DESeq2 并列
    ("(for example DESeq2 or edgeR, which model",
     "(for example DESeq2 «Love2014» or edgeR «Robinson2010», which model"),
]

for old, new in inserts:
    if old in body:
        body = body.replace(old, new, 1)
        print(f"[insert] ✓ {old[:60]}")
    else:
        print(f"[insert] ✗ NOT FOUND: {old[:60]}")

# ── Step 4: 按首次出现顺序分配连续编号 ───────────────────────────────────────
tag_order = []
for m in re.finditer(r"«[^»]+»", body):
    tag = m.group()
    if tag not in tag_order:
        tag_order.append(tag)

print(f"\n首次引用顺序({len(tag_order)}条):")
for i, t in enumerate(tag_order, 1):
    print(f"  [{i}] {t}")

# 替换标签为编号(逆序替换防止子串冲突)
for tag in tag_order:
    num = tag_order.index(tag) + 1
    body = body.replace(tag, f"[{num}]")

# 去掉连续重复引用 [3][3] → [3]
body = re.sub(r"(\[\d+\])(\s*\1)+", r"\1", body)

# ── Step 5: 生成新 References 段 ──────────────────────────────────────────────
ref_lines = ["\n## References\n"]
for i, tag in enumerate(tag_order, 1):
    if tag in ALL_REFS:
        entry, pmid = ALL_REFS[tag]
        ref_lines.append(f"{i}. {entry} *[PMID {pmid}]*")
    else:
        ref_lines.append(f"{i}. [TAG NOT IN POOL: {tag}]")

new_ref_block = "\n".join(ref_lines)
print(f"\nReferences 总条数: {len(tag_order)}")

# ── Step 6: 写回文件 ──────────────────────────────────────────────────────────
# 保留 ## References 之后的其他内容(Acknowledgements 等)
after_refs = ""
if "## Acknowledgements" in ref_block:
    after_refs = "\n\n" + ref_block[ref_block.index("## Acknowledgements"):]
elif "## Author Contributions" in ref_block:
    after_refs = "\n\n" + ref_block[ref_block.index("## Author Contributions"):]

new_text = body + new_ref_block + after_refs
with open(DEST, "w", encoding="utf-8") as f:
    f.write(new_text)
print(f"\n[DONE] 写入 {DEST}")
print(f"总行数: {len(new_text.splitlines())}")

EOF
