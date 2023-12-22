
# Load the structure
pdb_id = "1GZX"
cmd.fetch(pdb_id)

# Remove water molecules
cmd.remove('resn HOH')

# Colour ligands cyan
cmd.select('ligands', 'organic')
cmd.color('cyan', 'ligands')

# Colour binding sites grey
cmd.select('binding_sites', 'byres ligands around 5')
cmd.color('grey', 'binding_sites')

# Colour Tryptophan residues magenta
cmd.select('tryptophan', 'resn TRP')
cmd.color('magenta', 'tryptophan')

# Represent the ligand and binding site as sticks
cmd.show('sticks', 'ligands')
cmd.show('sticks', 'binding_sites')

# Represent the rest of the protein as cartoon
cmd.hide('everything', 'not ligands and not binding_sites')
cmd.show('cartoon', 'not ligands and not binding_sites')

# Colour the rest of the protein (not including ligands, binding sites, and Tryptophan residues) green
cmd.color('green', 'not ligands and not binding_sites and not tryptophan')

# export the result as a png file
cmd.png('./output.png', width=10cm , dpi=300, ray=1)



