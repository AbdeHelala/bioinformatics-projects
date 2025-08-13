
---

## README.md

# Protein Bioinformatics Tools

This repository contains **four independent Python scripts** for analyzing protein structures, sequences, and multiple sequence alignments. Each script has been designed for a specific task, with clean CLI usage, documented inputs/outputs, and reproducible workflows.

The scripts cover:

1. **`prot-stats.py`** – Composition & geometry analysis of PDB structures
2. **`pdb-hydropathy.py`** – Sliding-window Kyte–Doolittle hydropathy profiling & plotting
3. **`hydrocolor.py`** – PyMOL plugin to color residues by hydropathy scale
4. **`msa-logo.py`** – Generate PSSM and sequence logos from multiple sequence alignments

---

## 1️⃣ `prot-stats.py` — Protein Composition & Geometry Summary

A CLI tool to analyze a PDB file and report:

* Amino acid composition (count & % based on **Cα** atoms)
* Hydrophobicity breakdown (Kyte–Doolittle scale)
* Charge breakdown (positive/negative/neutral)
* Atomic composition by element
* HETATM residue counts (excluding water)
* Most distant Cα pair and distance
* Radius of gyration (Å, unit masses)

**Usage:**

```bash
python prot-stats.py path/to/structure.pdb
```

**Example output:**

```
Amino acid composition:
  ALA    15   5.76%
  LEU    20   7.69%
  ...
Hydrophobicity:
  Hydrophobic   85  48.02%
  Hydrophilic   92  51.98%
...
Most distant Cα residues: 5 ↔ 102 = 45.67 Å
Radius of gyration: 17.35 Å
```

---

## 2️⃣ `pdb-hydropathy.py` — Kyte–Doolittle Hydropathy Profiling

Generates sliding-window hydropathy plots for a protein chain from a PDB file.

**Features:**

* Extracts residues from a specified chain (default A)
* Computes **Kyte–Doolittle hydropathy values**
* Applies user-defined sliding window sizes
* Saves publication-quality PNG plots

**Usage:**

```bash
python pdb-hydropathy.py 1abc.pdb
```

**Example output files:**

* `Helala_hydrophobicity_plot_3.png`
* `Helala_hydrophobicity_plot_5.png`

---

## 3️⃣ `hydrocolor.py` — PyMOL Hydropathy Coloring Plugin

A PyMOL extension that **colors residues by hydrophobicity** directly inside PyMOL using the Kyte–Doolittle scale, with optional sliding window smoothing.

**Features:**

* Works on any selection (`all`, `chain A`, custom selections)
* Accepts **any PyMOL palette** (`blue_red`, `blue_white_red`, `rainbow`, etc.)
* Sliding-window smoothing with odd window sizes
* Writes hydropathy values to **B-factor field** for `spectrum` coloring
* Backwards-compatible alias `hydrophobicity` for older commands

**Installation & usage inside PyMOL:**

```pymol
run hydrocolor.py
hydrocolor chain A, blue_white_red, 5
```

or

```pymol
hydrophobicity all, rainbow, 3
```

---

## 4️⃣ `msa-logo.py` — Multiple Sequence Alignment Logo Generator

Converts a multiple sequence alignment (MSA) into:

1. A normalized **position-specific scoring matrix (PSSM)** CSV
2. A **sequence logo PNG** visualizing residue probabilities

**Features:**

* Supports FASTA, Clustal, and other Biopython-supported formats
* Normalizes per-position probabilities, ignoring gaps/unknowns
* Fixed **20-AA column order** in CSV & logo
* Generates both CSV and PNG in one run

**Usage:**

```bash
python msa-logo.py alignment.fasta --format fasta --outdir results --surname Helala
```

**Example output files:**

* `Helala_sequence_profile.csv`
* `Helala_sequence_logo.png`

---

## Installation

All scripts require Python 3.9+ and the following packages:

```bash
pip install biopython pandas logomaker matplotlib
```

`hydrocolor.py` requires [PyMOL](https://pymol.org) and must be run inside PyMOL.

---

## Project Structure

```
.
├── prot-stats.py          # PDB composition & geometry summarizer
├── pdb-hydropathy.py      # Kyte–Doolittle hydropathy plotter
├── hydrocolor.py          # PyMOL hydropathy coloring plugin
├── msa-logo.py            # MSA → PSSM + sequence logo generator
└── README.md              # This documentation
```

---

## Example Workflow

A typical workflow using all four scripts:

1. **Analyze a PDB file** for composition & geometry:

   ```bash
   python prot-stats.py protein.pdb
   ```

2. **Plot hydropathy profile**:

   ```bash
   python pdb-hydropathy.py protein.pdb
   ```

3. **Color protein in PyMOL**:

   ```pymol
   run hydrocolor.py
   hydrocolor protein, blue_red, 5
   ```

4. **Generate MSA logo** from multiple homologs:

   ```bash
   python msa-logo.py homologs.aln --format clustal
   ```

---

## Author

Scripts written and refined by **\[Abdelsalam / Helala]** for protein bioinformatics coursework, integrating:

* **Structural analysis** (PDB parsing & geometry metrics)
* **Sequence-based hydropathy profiling**
* **Interactive visualization in PyMOL**
* **Multiple sequence alignment visualization & scoring**

---
