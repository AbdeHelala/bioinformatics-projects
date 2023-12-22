from pymol import cmd, stored, util

kd = { 'ALA': 1.8,'ARG':-4.5,'ASN':-3.5,'ASP':-3.5,'CYS': 2.5,
        'GLN':-3.5,'GLU':-3.5,'GLY':-0.4,'HIS':-3.2,'ILE': 4.5,
        'LEU': 3.8,'LYS':-3.9,'MET': 1.9,'PHE': 2.8,'PRO':-1.6,
        'SER':-0.8,'THR':-0.7,'TRP':-0.9,'TYR':-1.3,'VAL': 4.2 }
@cmd.extend
def color_by_hydrophobicity(selection="all", palette="blue_red", window=1, _self=cmd):
    # Get the list of residues in the selection
    residues = cmd.get_model(selection).get_residues()
    # Compute average hydrophobicity
    hydrophobicity = compute_avg_hydrophobicity(residues, window)
    # Create a new color for each residue

    for i, res in enumerate(residues[:-window+1]):
        cmd.set_color("hydro_{}".format(i), util.cbay(hydrophobicity[i]))
        cmd.color("hydro_{}".format(i), "resi {}".format(res.resi))
    # Color the remaining residues with the last color
    for res in residues[-window+1:]:
        cmd.color("hydro_{}".format(len(hydrophobicity)-1), "resi {}".format(res.resi))
assert window % 2 == 1
# Register the function
cmd.extend("hydrophobicity", color_by_hydrophobicity)

    
    
    
    #TODO Write a PyMOL extension that will colour each residue in a given PyMOL selection according to the Kyte-Doolittle scale. The extension has to be imported
    # into PyMOL with run surname_matriculatonnumber_assignment2_3.py from the PyMOL command line, and has to be
    # launchable with the hydrophobicity [, selection [, palette]] from inside PyMOL. Refer to the PyMOL wiki and provided code snippet for guidance

