# -*- coding: utf-8 -*-
"""
prot-stats: PDB composition & geometry 

Short description: A tiny CLI that reads a PDB file and reports amino‑acid composition (counts & %), hydrophobicity and charge breakdowns, atomic composition, heteroatoms, most distant Cα pair, and radius of gyration.
prot-stats: PDB composition & geometry 
Usage:
    python prot_stats.py path/to/structure.pdb
"""

from __future__ import annotations
import argparse
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# Kyte–Doolittle hydropathy scale (3-letter)
KYTE_DOOLITTLE = {
    'ALA': 1.8,'ARG':-4.5,'ASN':-3.5,'ASP':-3.5,'CYS': 2.5,
    'GLN':-3.5,'GLU':-3.5,'GLY':-0.4,'HIS':-3.2,'ILE': 4.5,
    'LEU': 3.8,'LYS':-3.9,'MET': 1.9,'PHE': 2.8,'PRO':-1.6,
    'SER':-0.8,'THR':-0.7,'TRP':-0.9,'TYR':-1.3,'VAL': 4.2
}

AA_CHARGE = {
    'ARG': 'Positive', 'HIS': 'Positive', 'LYS': 'Positive',
    'ASP': 'Negative', 'GLU': 'Negative',
    'ALA': 'Neutral', 'CYS': 'Neutral', 'PHE': 'Neutral', 'GLY': 'Neutral', 'ILE': 'Neutral',
    'LEU': 'Neutral', 'MET': 'Neutral', 'ASN': 'Neutral', 'PRO': 'Neutral', 'GLN': 'Neutral',
    'SER': 'Neutral', 'THR': 'Neutral', 'VAL': 'Neutral', 'TRP': 'Neutral', 'TYR': 'Neutral'
}

# Common PDB non-standard residue normalizations
AA_NORMALIZE = {
    'MSE': 'MET',  # Selenomethionine
}

def read_pdb_lines(pdb_path: Path) -> List[str]:
    """Return all PDB lines (kept as raw strings)."""
    with pdb_path.open('r', encoding='utf-8', errors='ignore') as f:
        return [ln.rstrip('\n') for ln in f]

def is_atom(line: str) -> bool:
    return line.startswith('ATOM  ') or line.startswith('HETATM')

def is_het(line: str) -> bool:
    return line.startswith('HETATM')

def get_resname(line: str) -> str:
    # columns 18-20 (1-based), Python slice [17:20]
    res = line[17:20].strip().upper()
    return AA_NORMALIZE.get(res, res)

def get_atom_name(line: str) -> str:
    # columns 13-16
    return line[12:16].strip()

def get_resseq(line: str) -> Optional[int]:
    # columns 23-26
    try:
        return int(line[22:26].strip())
    except ValueError:
        return None

def get_xyz(line: str) -> Tuple[float, float, float]:
    # columns 31-38, 39-46, 47-54
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))

def get_element(line: str) -> str:
    # Prefer columns 77–78; if blank, derive from atom name
    elt = (line[76:78] if len(line) >= 78 else '').strip().upper()
    if elt:
        return elt
    name = get_atom_name(line)
    # Heuristic: first letter(s) of atom name (e.g. "CA" -> "C", "ZN" -> "ZN")
    if len(name) >= 2 and name[0].isalpha() and name[1].isalpha():
        return name[:2].upper()
    return name[:1].upper()

def iter_atom_lines(lines: Iterable[str]) -> Iterable[str]:
    for ln in lines:
        if is_atom(ln):
            yield ln

def iter_backbone_ca(lines: Iterable[str]) -> Iterable[str]:
    for ln in lines:
        if ln.startswith('ATOM  ') and get_atom_name(ln) == 'CA':
            yield ln

def amino_acid_composition(lines: List[str]) -> Dict[str, int]:
    comp: Dict[str, int] = {}
    for ln in iter_backbone_ca(lines):
        aa = get_resname(ln)
        if len(aa) == 3:
            comp[aa] = comp.get(aa, 0) + 1
    return comp

def to_percentages(counts: Dict[str, int], total: int) -> Dict[str, float]:
    if total == 0:
        return {k: 0.0 for k in counts}
    return {k: (v / total) * 100.0 for k, v in counts.items()}

def hydrophobicity_breakdown(aa_counts: Dict[str, int]) -> Tuple[Dict[str, int], Dict[str, float]]:
    buckets = {'Hydrophobic': 0, 'Hydrophilic': 0}
    total = 0
    for aa, n in aa_counts.items():
        kd = KYTE_DOOLITTLE.get(aa)
        if kd is None:
            continue  # skip unknowns
        buckets['Hydrophobic' if kd > 0 else 'Hydrophilic'] += n
        total += n
    return buckets, to_percentages(buckets, total)

def charge_breakdown(aa_counts: Dict[str, int]) -> Dict[str, int]:
    out = {'Positive': 0, 'Negative': 0, 'Neutral': 0}
    for aa, n in aa_counts.items():
        out[AA_CHARGE.get(aa, 'Neutral')] += n
    return out

def atomic_composition(lines: List[str]) -> Tuple[Dict[str, int], Dict[str, float]]:
    atoms: Dict[str, int] = {}
    atoms_lines = list(iter_atom_lines(lines))
    for ln in atoms_lines:
        el = get_element(ln)
        if el:
            atoms[el] = atoms.get(el, 0) + 1
    return atoms, to_percentages(atoms, len(atoms_lines))

def hetero_residues(lines: List[str], exclude_water: bool = True) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for ln in lines:
        if is_het(ln):
            res = get_resname(ln)
            if exclude_water and res in {'HOH', 'WAT'}:
                continue
            counts[res] = counts.get(res, 0) + 1
    return counts

def most_distant_ca_pair(lines: List[str]) -> Optional[Tuple[Tuple[int, int], float]]:
    ca_records = [(get_resseq(ln), get_xyz(ln)) for ln in iter_backbone_ca(lines)]
    ca_records = [(r, xyz) for r, xyz in ca_records if r is not None]
    if len(ca_records) < 2:
        return None
    max_d = -1.0
    best: Optional[Tuple[int, int]] = None
    for i in range(len(ca_records)):
        ri, (xi, yi, zi) = ca_records[i]
        for j in range(i + 1, len(ca_records)):
            rj, (xj, yj, zj) = ca_records[j]
            d = math.dist((xi, yi, zi), (xj, yj, zj))
            if d > max_d:
                max_d = d
                best = (ri, rj)
    return (best, max_d) if best else None

def center_of_mass(lines: List[str]) -> Tuple[float, float, float]:
    xs = ys = zs = total = 0.0
    for ln in iter_atom_lines(lines):
        x, y, z = get_xyz(ln)
        m = 1.0  # unit mass (no element masses used)
        xs += m * x; ys += m * y; zs += m * z; total += m
    if total == 0:
        return (0.0, 0.0, 0.0)
    return (xs / total, ys / total, zs / total)

def radius_of_gyration(lines: List[str]) -> float:
    com = center_of_mass(lines)
    total = 0.0
    accum = 0.0
    for ln in iter_atom_lines(lines):
        x, y, z = get_xyz(ln)
        dx = x - com[0]; dy = y - com[1]; dz = z - com[2]
        m = 1.0
        accum += m * (dx*dx + dy*dy + dz*dz)
        total += m
    return math.sqrt(accum / total) if total else 0.0

def pretty_print_report(pdb_path: Path, lines: List[str]) -> None:
    aa_counts = amino_acid_composition(lines)
    aa_total = sum(aa_counts.values())
    aa_pct = to_percentages(aa_counts, aa_total)

    hyd_counts, hyd_pct = hydrophobicity_breakdown(aa_counts)
    charge_counts = charge_breakdown(aa_counts)

    atom_counts, atom_pct = atomic_composition(lines)
    het_counts = hetero_residues(lines)

    distant = most_distant_ca_pair(lines)
    rog = radius_of_gyration(lines)

    print(f"File: {pdb_path.name}")
    print("\nAmino acid composition:")
    for aa in sorted(aa_counts):
        print(f"  {aa:>3}  {aa_counts[aa]:4d}  {aa_pct[aa]:6.2f}%")
    print(f"  Total residues (Cα): {aa_total}")

    print("\nHydrophobicity (Kyte–Doolittle):")
    for k in ('Hydrophobic', 'Hydrophilic'):
        print(f"  {k:<11} {hyd_counts.get(k,0):4d}  {hyd_pct.get(k,0.0):6.2f}%")

    print("\nCharge breakdown:")
    for k in ('Positive','Negative','Neutral'):
        print(f"  {k:<8} {charge_counts.get(k,0):4d}")

    print("\nAtomic composition (by element):")
    for el in sorted(atom_counts):
        print(f"  {el:<2}  {atom_counts[el]:6d}  {atom_pct[el]:6.2f}%")

    print(f"\nHETATM residues (non-water): {len(het_counts)} types")
    for res, n in sorted(het_counts.items()):
        print(f"  {res:<4} {n}")

    if distant:
        (r1, r2), d = distant
        print(f"\nMost distant Cα residues: {r1} ↔ {r2} = {d:.2f} Å")
    else:
        print("\nMost distant Cα residues: N/A")

    print(f"Radius of gyration: {rog:.2f} Å")

def main():
    ap = argparse.ArgumentParser(description="Summarize composition and geometry from a PDB file.")
    ap.add_argument("pdb", type=Path, help="Path to PDB file")
    args = ap.parse_args()

    if not args.pdb.exists():
        ap.error(f"PDB file not found: {args.pdb}")

    lines = read_pdb_lines(args.pdb)
    pretty_print_report(args.pdb, lines)

if __name__ == "__main__":
    main()
