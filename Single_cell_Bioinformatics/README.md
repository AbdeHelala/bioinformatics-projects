
---

# Single-Cell Bioinformatics Projects

This repository contains two analysis reports with executable RMarkdown notebooks:

* **Project 2 – scATAC-seq** (single-cell chromatin accessibility with ArchR)
* **Project 3 – Spatial Transcriptomics** (10x Visium with Seurat)

Each project directory includes the **.Rmd** (code + narrative) and a rendered **.pdf** report.

> Authors: **Abdelsalam Helala (ID: 7056985)**, **Ahmed Lamloum (ID: 7003029)**.&#x20;

---

## Repository structure

```
Single_cell_Bioinformatics/
├── Project 2- scATACseq/
│   ├── Project 2- scATACseq.Rmd
│   └── Project 2- scATACseq.pdf
└── SCB24_25__Spatial_Transcriptomics/
    ├── Project 3- Spatial Transcriptomics.Rmd
    └── Project 3- Spatial Transcriptomics.pdf

---

## Quick start (R)

```r
# Recommended R >= 4.2
install.packages(c("tidyverse","ggplot2","patchwork","dplyr","BiocManager"))

# ArchR & friends (for scATAC-seq)
BiocManager::install(c("ArchR","ComplexHeatmap"))
install.packages("Rmagic")     # optional smoothing

# Seurat & spatial ecosystem (for Visium)
install.packages(c("Seurat","SeuratData","hdf5r"))
BiocManager::install(c("CellChat","Biobase"))
# For deconvolution (optional):
# remotes::install_github("meichendong/SCDC")

# Render either project
rmarkdown::render("Project 2- scATACseq/Project 2- scATACseq.Rmd")
rmarkdown::render("SCB24_25__Spatial_Transcriptomics/Project 3- Spatial Transcriptomics.Rmd")
```

> **Data paths:** The notebooks assume local data layouts noted inside each `.Rmd`. Adjust `project_dir`/file paths as needed.

---

## Project 2 — scATAC-seq (ArchR pipeline)

### Overview

End-to-end single-cell ATAC-seq workflow: **preprocessing & QC → dimensionality reduction → batch correction → clustering → peak calling → gene activity → smoothing (MAGIC) → TF motif activity**.&#x20;

### Data & QC (highlights)

* **Input:** fragment files per sample; Arrow files built with **1 kb** tiling matrix.
* **Initial QC thresholds:** Fragments > 500; TSS enrichment > 4; doublets flagged via `addDoubletScores`.&#x20;
* **Stricter filter (final set):** Fragments **> 3,200 & < 100,000**, TSS **> 10**, doublet score **< 50** → **6,176** high-quality cells retained (from \~9,500).&#x20;
* **Doublets removed:** **327** total (cleaner downstream signal).&#x20;

**QC interpretation:**
Expected fragment-length peaks at \~147 bp (NFR) and \~300 bp (mono-nucleosome); clear TSS enrichment around 0 bp across samples.&#x20;

### DR, integration & clustering

* **Iterative LSI** (TF-IDF) preferred over PCA for sparse binary accessibility; **UMAP** for visualization.&#x20;
* **Harmony** batch correction: reduces sample-driven structure, improves biological mixing.&#x20;
* **Louvain clustering:** **12** clusters (C1–C12); largest: C8–C10; smallest: C1, C4.&#x20;

### Peaks, markers & gene activity

* **MACS2** peak calling on **cluster-grouped** coverages (reproducible peaks per cluster).
* **Marker peaks:** Wilcoxon test (bias on TSS & log10 fragments), **FDR ≤ 0.1**, **log2FC ≥ 0.5**; **139,142** significant peaks identified.&#x20;
* **Gene activity:** `GeneScoreMatrix` computed; marker genes by Wilcoxon with the same thresholds (heatmap summaries in report).&#x20;

### MAGIC smoothing (optional)

Graph-based imputation **clarifies sparse signals** for top marker genes, improving visualization across clusters.&#x20;

### TF motif activity (CIS-BP)

CIS-BP annotations added; variable motifs projected on UMAP and profiled across clusters. Example motifs: **TAL1\_62** and **TAL2\_822** show cluster-specific activity patterns.&#x20;

---

## Project 3 — Spatial Transcriptomics (10x Visium, Seurat)

### Overview

Spatial analysis on two Visium sections: **loading images & counts → QC → SCTransform → PCA/UMAP → clustering → differential expression & spatial features → merge & integration → (planned) cell-type ID & deconvolution → (planned) cell–cell communication**.&#x20;

### Platform notes (Visium)

* **Spot** diameter: **55 µm**; **center-to-center**: **100 µm**; **4,992 spots** per capture area (78×64).
  → Each spot aggregates expression from **multiple cells**; analyses capture **tissue domains** rather than single cells.&#x20;

### Loading & QC

* Images via `Read10X_Image`; counts via `Load10X_Spatial` for **Section\_1** and **Section\_2**.&#x20;
* QC metrics: `nFeature_Spatial`, `nCount_Spatial`, `% mitochondrial`. Example thresholds used:

  * Section 1: retain **306 / 3355** spots
  * Section 2: retain **938 / 3289** spots
    Using **nFeature\_Spatial 2000–7500**, **nCount\_Spatial 15k–50k**, \`%mt < 10–15\*\*.&#x20;

### Normalization, DR & clustering

* **SCTransform** on Spatial assay; **PCA** (chosen **30 PCs** from elbow), **UMAP**, **Louvain** clustering.

  * Section 1: **6** clusters; Section 2: **16** clusters (example resolutions and outputs shown in report).&#x20;

### Differential expression & spatial features

* `FindAllMarkers` identifies top markers per cluster (e.g., **Car8, Itpr1, Calb1, Pcp4, Ppp1r17** among others).&#x20;
* Spatially variable features via `FindSpatiallyVariableFeatures(selection.method="markvariogram")` and a small helper to retrieve ranked features; top examples include **Trh, Ttr, Ccn3** with strong spatial patterns.&#x20;

### Merging sections & batch correction

* **Without integration:** merged object shows section-specific structure; cluster/sample cross-tabs quantify overlap.&#x20;
* **With integration:** Seurat integration (anchors + `IntegrateData`) improves **mixing** across sections; recommended for downstream tasks.&#x20;

### (Planned / partially implemented) downstreams

* **Cell-type ID**: label transfer from a reference (Allen cortex); **manual** curation via DE markers.
* **Deconvolution**: SCDC to estimate **cell type proportions per spot** using scRNA-seq reference.&#x20;
* **Cell–cell communication**: **CellChat** pipeline to infer ligand–receptor networks and visualize pathways (e.g., WNT).&#x20;

> Note: Some code blocks in the PDF show placeholder objects (e.g., `merged_integrated`) to illustrate intended steps; ensure the integrated object is created and stored before running label transfer, SCDC, and CellChat sections.&#x20;

---

## Reproducibility & tips

* Set seeds where provided (e.g., `set.seed(42)`), and keep R/bioc package versions consistent.&#x20;
* For scATAC-seq, ensure **MACS2** is installed and on `PATH` for peak calling.&#x20;
* Spatial data require **Visium outputs** (Space Ranger image & HDF5 matrix); update `project_dir` in the `.Rmd` to your local paths.&#x20;

---
