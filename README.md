# Bioinformatics Projects

A curated collection of my coursework and projects across several bioinformatics tracks, covering **algorithms**, **structural bioinformatics**, **single-cell omics**, and **neural networks for molecular property prediction**. Each subfolder is a self-contained project with its own code and (where relevant) report/notebook.

> Author: **Abdelsalam Helala**
> Academic year: 2024–2025

---

## Repository layout

```
bioinformatics-projects/
├── Neural-Networks/
│   └── Molecular Property Prediction/        # NN project on SMILES → lipophilicity
├── Single_cell_Bioinformatics/               # scATAC-seq + Spatial transcriptomics (RMarkdown + PDFs)
├── Structural_Bioinformatics/                # CLI tools + PyMOL plugin + MSA logo
└── python2024-main/                          # Algorithmic Bioinformatics assignments + AFM project
```

---

## Quick setup

### Python

* Recommended: **Python 3.10+**
* Common scientific stack (project-dependent):

  ```bash
  pip install numpy matplotlib pandas biopython logomaker streamlit plotly
  ```
* For the PyMOL plugin in *Structural\_Bioinformatics*: install **PyMOL** separately.

### R / Bioconductor (for *Single\_cell\_Bioinformatics*)

* Recommended: **R ≥ 4.2**
* Core packages (project-dependent):

  ```r
  install.packages(c("tidyverse","ggplot2","patchwork","dplyr","Seurat","SeuratData","hdf5r"))
  BiocManager::install(c("ArchR","ComplexHeatmap","Biobase","CellChat"))
  # optional utilities
  install.packages("Rmagic")
  # remotes::install_github("meichendong/SCDC")  # for deconvolution, optional
  ```

> Large artifacts (e.g., PDFs, models) are tracked with Git LFS where needed.

---

## Projects

### 1) Structural\_Bioinformatics

Utility scripts and a PyMOL extension for protein structure/sequence analysis.

* **`prot-stats.py`** — Summarize a PDB:

  * AA composition (counts & % via Cα), hydrophobic/hydrophilic (Kyte–Doolittle), charge classes, atomic composition, HETATM residue counts, most distant Cα pair, **radius of gyration**.
  * Run: `python prot-stats.py path/to/structure.pdb`

* **`pdb-hydropathy.py`** — Sliding-window hydropathy plot (Kyte–Doolittle):

  * Extracts a chain and saves PNG plots for windows (e.g., 3/5/7/9).
  * Run: `python pdb-hydropathy.py path/to/structure.pdb`

* **`hydrocolor.py`** — **PyMOL** plugin to color residues by hydropathy:

  * Usage inside PyMOL:

    ```
    run hydrocolor.py
    hydrocolor <selection>, <palette>, <odd window>
    # e.g., hydrocolor chain A, blue_white_red, 5
    ```

* **`msa_logo.py`** — Build a normalized PSSM and a sequence logo from an MSA:

  * Outputs: `<surname>_sequence_profile.csv` and `<surname>_sequence_logo.png`
  * Run:
    `python msa_logo.py alignment.fasta --format fasta --outdir results --surname Helala`

> Dependencies: `biopython`, `pandas`, `logomaker`, `matplotlib` (and PyMOL for the plugin).

---

### 2) Single\_cell\_Bioinformatics

Two RMarkdown projects with rendered PDFs (methods + figures).

**Project 2 — scATAC-seq (ArchR pipeline)**

* **Workflow:** preprocessing/QC → iterative LSI → Harmony batch correction → UMAP → Louvain clustering → cluster-wise **MACS2** peak calling → marker peaks/genes → gene activity (`GeneScoreMatrix`) → optional **MAGIC** smoothing → motif activity (CIS-BP).
* **QC highlights:** fragment count & TSS enrichment thresholds; doublet scoring/removal; nucleosome signal and TSS profiles checked.
* **Outputs:** UMAPs (raw/integrated), QC summaries, marker heatmaps, motif activity views.

**Project 3 — Spatial Transcriptomics (10x Visium, Seurat)**

* **Workflow:** `Load10X_Spatial` + `Read10X_Image` → QC (spot filters on features/counts/%MT) → **SCTransform** → PCA/UMAP → clustering → differential expression → **spatially variable features** → merge sections → Seurat **integration** (anchors + `IntegrateData`).
* **Planned downstreams:** label transfer (cell-type ID), **SCDC** deconvolution, **CellChat** receptor–ligand analysis.

> Render notebooks:
> `rmarkdown::render("Project 2- scATACseq/Project 2- scATACseq.Rmd")`
> `rmarkdown::render("SCB24_25__Spatial_Transcriptomics/Project 3- Spatial Transcriptomics.Rmd")`

---

### 3) Algorithmic Bioinformatics — `python2024-main`

Course assignments (01–04) with example tests and a multi-part **AFM** project.

**Assignments (with public example tests):**

* **A01:** list processing; DNA reverse complement; k-mers utilities; **k-mer base-4 encoding/decoding** via bit ops.
* **A02:** **anagram grouping**; **Hamming distance** across a file (RNA→DNA normalization & validation); **common k-mers** from FASTA; **vampire numbers** + generator.
* **A03:** parcel stack rearrangement parser/simulator; **PWM** & all possible **consensus sequences**; tiny arithmetic interpreter using `match/case`.
* **A04 (NumPy):** **Counting Bloom Filter** (insert/delete/contains/count); circular list “treasure” puzzle; **quantile normalization**.

Run tests (example):

```bash
cd python2024-main/assignments
pytest -q  # or python test_assignment_02.py
```

**AFM Project (Parts I–III):**

* **Part I:** parse AFM text, extract distance (m) & force (N), plot **f vs d** per (series, i, j).
  `python project-afm/plotafm.py --textfile sample.txt --show --plotprefix curve`
* **Part II:** estimate **slope** in the linear left region (esp. series 0), overlay fit, and print `s i j slope` lines suitable for automated grading.
  `python project-afm/plotafm.py -t sample.txt --plotprefix curve`
* **Part III:** **Streamlit** app showing **heights** (heatmap) and **per-pixel slopes** (heatmap), with an interactive curve viewer at (s,i,j).
  `streamlit run project-afm/afm.py -- afm`  *(expects `afm.heights.npy` & `afm.data.pickled`)*

---

### 4) Neural-Networks / Molecular Property Prediction

Neural models for **SMILES → lipophilicity** regression.

* **Task 1:** baseline supervised fine-tuning + **unsupervised MLM pre-finetuning** to improve downstream performance.
* **Task 2:** **influence-function** data selection (LiSSA inverse-HVP) to pick impactful external samples and boost generalization.
* **Task 3:** **PEFT** methods — **BitFit**, **LoRA**, **IA³** — to reduce trainable params while preserving accuracy.

> Environment: PyTorch, Hugging Face Transformers, scikit-learn; optional `peft` and WandB.

---

## Data notes

* Large raw datasets (e.g., full AFM text, Visium outputs) are **not committed**. Paths are parameterized in scripts/`.Rmd` files; place data locally as described there.
* Figures/reports are included as PDFs where appropriate (LFS-tracked).

---
