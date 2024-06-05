import sys
import matplotlib.pyplot as plt
import os

#TODO write your surname here to name file in compute_hydrophobicity function
STUDENT_SURNAME = 'Helala'

residues = []
def get_residues(pdb_file):

    with open(pdb_file, 'r') as f:
        lines = [l for l in f if l.startswith('ATOM')]
        residues = [l[17:20].strip() for l in lines if l[21] == 'A']
    return residues

kd = { 'ALA': 1.8,'ARG':-4.5,'ASN':-3.5,'ASP':-3.5,'CYS': 2.5,
        'GLN':-3.5,'GLU':-3.5,'GLY':-0.4,'HIS':-3.2,'ILE': 4.5,
        'LEU': 3.8,'LYS':-3.9,'MET': 1.9,'PHE': 2.8,'PRO':-1.6,
        'SER':-0.8,'THR':-0.7,'TRP':-0.9,'TYR':-1.3,'VAL': 4.2 }

hydrophobicity = []
def compute_hydrophobicity(residues, window_size=3):

    hydrophobicity = [kd[res] for res in residues if res in kd]
    return [sum(hydrophobicity[i:i+window_size])/window_size for i in range(len(hydrophobicity)-window_size+1)]    
  
avg = []
def compute_avg_hydrophobicity(residues, window_size=3):
    hydrophobicity_smooth = []
    for i in range(len(hydrophobicity) - window_size + 1):
        window = hydrophobicity[i:i+window_size]
        avg = sum(window) / window_size
        hydrophobicity_smooth.append(avg)
    return hydrophobicity_smooth

    #TODO calculate average value with window size and store it as a list in average_values

    plt.plot(avg)
    plt.xlabel("AA position")
    plt.ylabel("Hydrophobicity")
    plt.savefig(f"./{STUDENT_SURNAME}_hydrophobicity_plot_{window_size}.png")
    plt.clf()

# Nothing to do here
if __name__ == "__main__":
    if STUDENT_SURNAME == None:
        print(f"Helala")
        sys.exit()
    residues = get_residues(sys.argv[1])
    print(f"Residues: {residues}")
    print(f"Length of residues: {len(residues)}")
    window_sizes = [5, 9]
    for window_size in window_sizes:
        compute_hydrophobicity(residues, window_size)
        plt.clf()

    for window_size in [3, 7]:
        compute_hydrophobicity(residues, window_size)
    os.remove(f"./{STUDENT_SURNAME}_hydrophobicity_plot_3.png")

