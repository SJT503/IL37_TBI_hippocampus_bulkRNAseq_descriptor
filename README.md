# IL-37 TBI Hippocampus Bulk RNA-seq Data Descriptor — Analysis Code

**Companion code repository for:**  
Lin R†, Liu X†, Sun Q, Chen X, Li K, Chen W\*, Sheng J\*  
"A bulk RNA sequencing dataset of the mouse hippocampus after controlled cortical impact with AAV-mediated human interleukin-37 overexpression"  
*Scientific Data* (submitted 2026)

---

## Overview

This repository contains all analysis code used to compute the technical-validation metrics and to generate the figures and tables reported in the above Data Descriptor, starting from the released gene-level count matrix deposited in GEO (accession GSE337710).

The dataset covers nine bulk RNA-seq libraries from the ipsilateral hippocampus of C57BL/6NJ male mice across three conditions (Sham / CCI+RFP / CCI+IL-37, n = 3 per group), sequenced on an Illumina NovaSeq 6000.

---

## Repository structure

```
scripts/
  novogene_qc_stats.py        # Ground-truth QC constants (Novogene report)
  extended_stats.py           # Extended analysis constants (AS/SNV/novel transcripts)
  01_technical_validation.py  # Main pipeline: tables + all figures
  make_figures.py             # Main figure generation (Fig 1–6)
  make_fig1_design.py         # Fig 1 study design schematic
  supplementary.py            # Supplementary tables S1–S5 + figures S1–S3
  02_build_submission.py      # Build submission package (docx + figures)
  03_rebuild_references.py    # Reference renumbering utility
  04_add_refs.py              # Reference insertion utility
```

---

## How to reproduce all figures and tables

### Requirements
```
Python >= 3.9
pandas, numpy, matplotlib, scikit-learn, python-docx, openpyxl, pandoc (CLI)
```

Install:
```bash
pip install pandas numpy matplotlib scikit-learn python-docx openpyxl
```

### Data
Download the count matrix from GEO (GSE337710):
- `bulk_counts.tsv` — 41,045 genes × 9 samples (raw integer counts)

Place it at: `F:/IL-37/analysis/deconvolution/results/bulk_counts.tsv`  
(or update the `COUNTS` path in `01_technical_validation.py`)

### Run
```bash
# Generate all tables + figures
python scripts/01_technical_validation.py

# Outputs:
# results/tables/Table1_sequencing_qc.csv
# results/tables/Table2_sample_info.csv
# results/tables/Table3_marker_response.csv
# results/tables/Table4_deg_summary.csv
# results/tables/supplementary/TableS1_DEG_full_lists.csv
# results/tables/supplementary/TableS2_logCPM_all_genes.csv
# results/tables/supplementary/TableS3_marker_CPM_per_sample.csv
# results/tables/supplementary/TableS4_AS_events.csv
# results/tables/supplementary/TableS5_SNP_annotation.csv
# results/figures/Fig1_study_design.png  (requires make_fig1_design.py)
# results/figures/Fig2_seq_qc.png
# results/figures/Fig3_pca.png
# results/figures/Fig4_correlation.png
# results/figures/Fig5_marker_response.png
# results/figures/Fig6_volcano.png
# results/figures/supplementary/FigS1_gene_detection.png
# results/figures/supplementary/FigS2_AS_events.png
# results/figures/supplementary/FigS3_variant_annotation.png
```

---

## Data sources

All quantitative inputs used by the scripts are derived from:
1. **Novogene official QC report** — hardcoded in `novogene_qc_stats.py` (Q30, mapping rates, insert size, correlation, PCA)
2. **Gene-level count matrix** (`bulk_counts.tsv`) — primary processed data, deposited in GEO
3. **Novogene DEG results** — provider DESeq2 output, used for Table 4 and Supplementary Table S1
4. **Extended analysis outputs** — rMATS AS events, GATK+SnpEff SNV/Indel counts, StringTie novel transcripts; summarised in `extended_stats.py`

---

## Citation

If you use this code or the dataset, please cite:  
> Lin R, Liu X, Sun Q, Chen X, Li K, Chen W\*, Sheng J\*. A bulk RNA sequencing dataset of the mouse hippocampus after controlled cortical impact with AAV-mediated human interleukin-37 overexpression. *Scientific Data* (2026). https://doi.org/10.5281/zenodo.21214811 (code)

GEO accession: **GSE337710** (pending public release)

---

## License

MIT License. See LICENSE file.

---

## Contact

Jiangtao Sheng — jtsheng@stu.edu.cn  
Department of Microbiology and Immunology, Shantou University Medical College
