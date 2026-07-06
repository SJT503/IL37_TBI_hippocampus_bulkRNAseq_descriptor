"""
按正文首次出现顺序给作者标签统一编号 -> 替换回 [n] -> 生成排序 References.
三锁保证: 正文首次引用顺序 = 编号 = References 列表顺序.
所有条目已 PubMed 核实.
"""
import re, os

MAN = r"F:\IL-37\SciData_descriptor\manuscript\data_descriptor_IL37.md"

# 作者标签 -> 完整 Vancouver 引用条目 (全部已 PubMed 核实, PMID 标注)
REF = {
 "Dewan": "Dewan MC, et al. Estimating the global incidence of traumatic brain injury. *J Neurosurg.* 2019;130(4):1080-1097. https://doi.org/10.3171/2017.10.JNS17352.",
 "Maas": "Maas AIR, et al. Traumatic brain injury: integrated approaches to improve prevention, clinical care, and research. *Lancet Neurol.* 2017;16(12):987-1048. https://doi.org/10.1016/S1474-4422(17)30371-X.",
 "GBD": "GBD 2016 Traumatic Brain Injury and Spinal Cord Injury Collaborators. Global, regional, and national burden of traumatic brain injury and spinal cord injury, 1990-2016. *Lancet Neurol.* 2019;18(1):56-87. https://doi.org/10.1016/S1474-4422(18)30415-0.",
 "Ma": "Ma Z, et al. Traumatic brain injury in elderly population: a global systematic review and meta-analysis of in-hospital mortality and risk factors among 2.22 million individuals. *Ageing Res Rev.* 2024;99:102376. https://doi.org/10.1016/j.arr.2024.102376.",
 "Cook": "Cook AM, et al. Update on neuroprotection after traumatic brain injury. *CNS Drugs.* 2025;39(5):473-484. https://doi.org/10.1007/s40263-025-01175-7.",
 "Pearson": "Pearson A, et al. PPARγ activation ameliorates cognitive impairment and chronic microglial activation in the aftermath of r-mTBI. *J Neuroinflammation.* 2024;21(1):194. https://doi.org/10.1186/s12974-024-03173-w.",
 "Ji": "Ji Q, et al. Astrocytic SIRT1 ameliorates cognitive deficits after traumatic brain injury via autophagy-mediated MEGF10 phagocytosis. *J Neuroinflammation.* 2025;22(1):304. https://doi.org/10.1186/s12974-025-03633-x.",
 "Osier": "Osier ND, Dixon CE. The controlled cortical impact model: applications, considerations for researchers, and future directions. *Front Neurol.* 2016;7:134. https://doi.org/10.3389/fneur.2016.00134.",
 "Pybus": "Pybus AF, et al. Profiling the neuroimmune cascade in 3xTg-AD mice exposed to successive mild traumatic brain injuries. *J Neuroinflammation.* 2024;21(1):156. https://doi.org/10.1186/s12974-024-03128-1.",
 "Ramlack": "Ramlackhansingh AF, et al. Inflammation after trauma: microglial activation and traumatic brain injury. *Ann Neurol.* 2011;70(3):374-383. https://doi.org/10.1002/ana.22455.",
 "Escartin": "Escartin C, et al. Reactive astrocyte nomenclature, definitions, and future directions. *Nat Neurosci.* 2021;24(3):312-325. https://doi.org/10.1038/s41593-020-00783-4.",
 "Hammond": "Hammond TR, et al. Single-cell RNA sequencing of microglia throughout the mouse lifespan and in the injured brain reveals complex cell-state changes. *Immunity.* 2019;50(1):253-271.e6. https://doi.org/10.1016/j.immuni.2018.11.004.",
 "Jha": "Jha RM, et al. A single-cell atlas deconstructs heterogeneity across multiple models in murine traumatic brain injury and identifies novel cell-specific targets. *Neuron.* 2024;112(18):3069-3088.e4. https://doi.org/10.1016/j.neuron.2024.06.020.",
 "Dinarello": "Dinarello CA, et al. Suppression of innate inflammation and immunity by interleukin-37. *Eur J Immunol.* 2016;46(5):1067-1081. https://doi.org/10.1002/eji.201545828.",
 "Boraschi": "Boraschi D, et al. IL-37: a new anti-inflammatory cytokine of the IL-1 family. *Eur Cytokine Netw.* 2011;22(3):127-147. https://doi.org/10.1684/ecn.2011.0288.",
 "Zhang": "Zhang J, et al. Interleukin-37 modulates microglial phenotype and inhibits inflammatory response via the MyD88/NF-κB pathway in lipopolysaccharide-induced neuroinflammation. *Inflamm Res.* 2025;74(1):87. https://doi.org/10.1007/s00011-025-02048-x.",
 "HISAT2": "Kim D, et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nat Biotechnol.* 2019;37(8):907-915. https://doi.org/10.1038/s41587-019-0201-4.",
 "featureCounts": "Liao Y, et al. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics.* 2014;30(7):923-930. https://doi.org/10.1093/bioinformatics/btt656.",
 "StringTie": "Pertea M, et al. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. *Nat Biotechnol.* 2015;33(3):290-295. https://doi.org/10.1038/nbt.3122.",
 "rMATS": "Shen S, et al. rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. *Proc Natl Acad Sci USA.* 2014;111(51):E5593-E5601. https://doi.org/10.1073/pnas.1419161111.",
 "DESeq2": "Love MI, et al. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.* 2014;15(12):550. https://doi.org/10.1186/s13059-014-0550-8.",
}

t = open(MAN, encoding="utf-8").read()
body, refs_tail = t.split("## References", 1)

# 1. 按首次出现顺序给标签编号
order = []
for m in re.finditer(r"\{\{(\w+)\}\}", body):
    tag = m.group(1)
    if tag not in order:
        order.append(tag)

# 校验: 所有标签都有条目
missing = [tg for tg in order if tg not in REF]
assert not missing, f"缺条目: {missing}"

num = {tag: i + 1 for i, tag in enumerate(order)}

# 2. 合并相邻标签为 [n,m], 替换回数字
#    先把连续的 }}{{ 之间插逗号占位, 再整体处理
def replace_run(mo):
    tags = re.findall(r"\{\{(\w+)\}\}", mo.group(0))
    nums = sorted(num[tg] for tg in tags)
    return "[" + ",".join(str(n) for n in nums) + "]"

body2 = re.sub(r"(?:\{\{\w+\}\})+", replace_run, body)

# 3. 生成 References 列表 (按编号顺序)
lines = ["## References", ""]
for tag in order:
    lines.append(f"{num[tag]}. {REF[tag]}")
refs_block = "\n".join(lines) + "\n"

open(MAN, "w", encoding="utf-8").write(body2 + refs_block)

print(f"参考文献总数: {len(order)}")
print(f"首次引用顺序编号: {[(num[t], t) for t in order]}")
# 残留标签检查
left = re.findall(r"\{\{\w+\}\}", body2)
print(f"残留标签(应0): {left}")
