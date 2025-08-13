# hydrocolor.py
# PyMOL extension: color residues by Kyte–Doolittle hydropathy (with optional smoothing window).
#
# Usage inside PyMOL:
#   run hydrocolor.py
#   hydrocolor [selection [, palette [, window]]]
#
# Examples:
#   hydrocolor all, blue_white_red, 1
#   hydrocolor chain A, blue_red, 5
#
# Palettes to try: blue_white_red, blue_red, green_red, red_white_blue, rainbow

from pymol import cmd
from typing import List, Tuple

# Kyte–Doolittle hydropathy scale
KD = {
    'ALA': 1.8,'ARG':-4.5,'ASN':-3.5,'ASP':-3.5,'CYS': 2.5,
    'GLN':-3.5,'GLU':-3.5,'GLY':-0.4,'HIS':-3.2,'ILE': 4.5,
    'LEU': 3.8,'LYS':-3.9,'MET': 1.9,'PHE': 2.8,'PRO':-1.6,
    'SER':-0.8,'THR':-0.7,'TRP':-0.9,'TYR':-1.3,'VAL': 4.2
}

# Common normalization for non‑standard residues
AA_NORMALIZE = { 'MSE': 'MET' }

def _residue_list(selection: str) -> List[Tuple[str, str, str, str]]:
    """
    Return an ordered list of residues in the selection as tuples:
    (segi, chain, resi, resn)
    """
    model = cmd.get_model(selection)
    seen = set()
    residues = []
    for a in model.atom:
        key = (a.segi, a.chain, a.resi, a.resn.upper())
        if key not in seen:
            seen.add(key)
            residues.append(key)
    return residues

def _kd_value(resn: str) -> float:
    resn = AA_NORMALIZE.get(resn.upper(), resn.upper())
    return KD.get(resn, 0.0)  # unseen residues get 0.0

def _smooth(vals: List[float], window: int) -> List[float]:
    """Odd-length sliding average; if len<window, returns empty list."""
    if window < 1 or window % 2 == 0:
        raise ValueError("window must be a positive odd integer")
    n = len(vals)
    if n < window:
        return []
    w = window
    return [sum(vals[i:i+w]) / w for i in range(0, n - w + 1)]

@cmd.extend
def hydrocolor(selection: str = "all", palette: str = "blue_red", window: int = 1):
    """
    Color residues by Kyte–Doolittle hydropathy.
      selection: PyMOL selection (default 'all')
      palette  : palette for 'spectrum b' (e.g., blue_red, blue_white_red, green_red, rainbow)
      window   : odd sliding window size for smoothing (1 = no smoothing)
    """
    try:
        window = int(window)
    except Exception:
        print("hydrocolor: 'window' must be an integer.")
        return

    if window < 1 or window % 2 == 0:
        print("hydrocolor: please provide an ODD window size >= 1 (e.g., 1, 3, 5, 7).")
        return

    residues = _residue_list(selection)
    if not residues:
        print("hydrocolor: no residues found in selection.")
        return

    kd_vals = [_kd_value(resn) for (_, _, _, resn) in residues]

    # Optional smoothing with odd window
    if window > 1:
        smoothed = _smooth(kd_vals, window)
        if not smoothed:
            print("hydrocolor: window larger than sequence length; using no smoothing.")
            smoothed = kd_vals
            window = 1
        k = (window - 1) // 2
        per_res_vals = [None] * len(residues)
        for i, v in enumerate(smoothed):
            per_res_vals[i + k] = v
        # Fill edges with nearest available value
        first_val = next((v for v in per_res_vals if v is not None), 0.0)
        last_val = next((v for v in reversed(per_res_vals) if v is not None), first_val)
        for i in range(len(per_res_vals)):
            if per_res_vals[i] is None:
                per_res_vals[i] = first_val if i < k else last_val
    else:
        per_res_vals = kd_vals

    # Write values into B-factor (b) and color via spectrum
    for (segi, chain, resi, _), val in zip(residues, per_res_vals):
        sel_res = f"({selection}) and resi {resi}"
        if chain:
            sel_res += f" and chain {chain}"
        if segi:
            sel_res += f" and segi {segi}"
        cmd.alter(sel_res, f"b={float(val)}")

    cmd.rebuild()
    try:
        cmd.spectrum("b", palette, selection)
    except Exception:
        print(f"hydrocolor: unknown palette '{palette}'. Falling back to blue_red.")
        cmd.spectrum("b", "blue_red", selection)

    print(f"hydrocolor: colored '{selection}' by Kyte–Doolittle (window={window}, palette='{palette}').")

# Backwards-compat command alias: 'hydrophobicity'
@cmd.extend
def hydrophobicity(selection: str = "all", palette: str = "blue_red", window: int = 1):
    """Alias for hydrocolor() to keep older notes/commands working."""
    return hydrocolor(selection, palette, window)
