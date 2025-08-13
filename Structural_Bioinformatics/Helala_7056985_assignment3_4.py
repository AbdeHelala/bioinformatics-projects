#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
msa-logo: build a normalized PSSM and sequence logo from a multiple sequence alignment.

Examples:
  python msa_logo.py alignment.fasta
  python msa_logo.py alignment.aln --format clustal --outdir results --surname Helala
"""

import sys
import csv
import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd
from Bio import AlignIO
from Bio.Align import AlignInfo
import logomaker

# Standard 20 AA order (columns in CSV/logo)
AA20 = list("ACDEFGHIKLMNPQRSTVWY")

def calculate_sequence_profiles(msa_path: Path, fmt: str = "fasta") -> List[Dict[str, float]]:
    """
    Read an alignment and return a list of per-position probability dicts over AA20.
    Non-AA (e.g., gaps, X, B, Z) are ignored in normalization.
    """
    alignment = AlignIO.read(msa_path, fmt)
    summary = AlignInfo.SummaryInfo(alignment)

    # Biopython's PSSM here is counts by residue per position
    raw_pssm = summary.pos_specific_score_matrix()  # List[Dict[str, float]]

    profiles: List[Dict[str, float]] = []
    for pos_dict in raw_pssm:
        # Keep only standard AAs
        counts = {aa: float(pos_dict.get(aa, 0.0)) for aa in AA20}
        total = sum(counts.values())
        if total == 0:
            # all gaps — leave zeros
            probs = {aa: 0.0 for aa in AA20}
        else:
            probs = {aa: counts[aa] / total for aa in AA20}
        profiles.append(probs)
    return profiles

def write_pssm_csv(pssm: List[Dict[str, float]], out_csv: Path) -> None:
    """Write PSSM (probabilities) to CSV with columns: Position + AA20 in fixed order."""
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Position"] + AA20)
        for i, pos in enumerate(pssm, start=1):
            writer.writerow([i] + [pos.get(aa, 0.0) for aa in AA20])

def sequence_logo_creator(pssm_csv: Path, out_png: Path) -> None:
    """Read the CSV and save a sequence logo PNG."""
    df = pd.read_csv(pssm_csv)
    # Ensure AA columns present and in correct order
    for aa in AA20:
        if aa not in df.columns:
            df[aa] = 0.0
    profile = df[AA20]  # drop Position, keep fixed AA order

    logo = logomaker.Logo(profile)
    logo.ax.set_xlabel("Position")
    logo.ax.set_ylabel("Probability")
    logo.fig.savefig(out_png, dpi=300, bbox_inches="tight")
    # no plt.close needed; logomaker manages its own fig

def main():
    ap = argparse.ArgumentParser(description="Generate a normalized PSSM and sequence logo from an MSA.")
    ap.add_argument("alignment", type=Path, help="Path to alignment file (FASTA/Clustal/etc.)")
    ap.add_argument("--format", default="fasta",
                    help="Alignment format for Biopython AlignIO (default: fasta, e.g., 'clustal')")
    ap.add_argument("--outdir", type=Path, default=Path("."),
                    help="Output directory (default: current dir)")
    ap.add_argument("--surname", default="Helala",
                    help="Surname/prefix for output filenames (default: Helala)")
    args = ap.parse_args()

    if not args.alignment.exists():
        sys.exit(f"Alignment file not found: {args.alignment}")

    # Compute profiles
    pssm = calculate_sequence_profiles(args.alignment, fmt=args.format)

    # Paths
    csv_path = args.outdir / f"{args.surname}_sequence_profile.csv"
    png_path = args.outdir / f"{args.surname}_sequence_logo.png"

    # Write CSV + logo
    write_pssm_csv(pssm, csv_path)
    sequence_logo_creator(csv_path, png_path)

    print(f"Wrote PSSM CSV:   {csv_path}")
    print(f"Wrote logo PNG:   {png_path}")

if __name__ == "__main__":
    main()
