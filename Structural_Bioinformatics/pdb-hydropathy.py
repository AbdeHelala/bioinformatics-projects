#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hydrophobicity_plot.py
Compute and plot sliding-window hydrophobicity from a PDB file.

Usage:
    python hydrophobicity_plot.py path/to/file.pdb
"""

import sys
import os
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List

# Write your surname here to name the output plot file
STUDENT_SURNAME = "Helala"

# Kyte–Doolittle hydropathy scale
KYTE_DOOLITTLE = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
    'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
    'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
    'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
}

def get_residues(pdb_file: Path) -> List[str]:
    """
    Extracts the three-letter residue codes from chain A of a PDB file.
    """
    residues = []
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == "A":  # chain A only
                res = line[17:20].strip().upper()
                residues.append(res)
    return residues

def compute_avg_hydrophobicity(residues: List[str], window_size: int) -> List[float]:
    """
    Computes sliding-window average hydrophobicity.
    """
    # Filter residues to ones in the KD scale
    values = [KYTE_DOOLITTLE[res] for res in residues if res in KYTE_DOOLITTLE]
    if len(values) < window_size:
        return []
    return [
        sum(values[i:i+window_size]) / window_size
        for i in range(len(values) - window_size + 1)
    ]

def plot_hydrophobicity(hydro_values: List[float], window_size: int) -> Path:
    """
    Saves a hydrophobicity plot for the given values.
    """
    plt.figure()
    plt.plot(hydro_values, marker="o")
    plt.xlabel("AA position")
    plt.ylabel("Hydrophobicity")
    plt.title(f"Hydrophobicity (window size={window_size})")
    out_path = Path(f"{STUDENT_SURNAME}_hydrophobicity_plot_{window_size}.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    return out_path

if __name__ == "__main__":
    if not STUDENT_SURNAME:
        sys.exit("Please set STUDENT_SURNAME before running.")

    if len(sys.argv) < 2:
        sys.exit("Usage: python hydrophobicity_plot.py path/to/file.pdb")

    pdb_file = Path(sys.argv[1])
    if not pdb_file.exists():
        sys.exit(f"PDB file not found: {pdb_file}")

    residues = get_residues(pdb_file)
    print(f"Extracted {len(residues)} residues from chain A")

    # Example: generate plots for multiple window sizes
    for window_size in [3, 5, 7, 9]:
        hydro_vals = compute_avg_hydrophobicity(residues, window_size)
        if hydro_vals:
            plot_file = plot_hydrophobicity(hydro_vals, window_size)
            print(f"Saved: {plot_file}")
        else:
            print(f"Window size {window_size} too large for sequence length")

    # Example: remove a specific file if needed
    file_to_remove = Path(f"{STUDENT_SURNAME}_hydrophobicity_plot_3.png")
    if file_to_remove.exists():
        os.remove(file_to_remove)
        print(f"Removed: {file_to_remove}")
