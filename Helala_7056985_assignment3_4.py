import sys
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
import logomaker
from Bio import AlignIO
import csv
import os
import pandas as pd

#TODO write your surname here to name file in the outputs
STUDENT_SURNAME = "Helala"

def calculate_sequence_profiles(clustal_omega_file):
    """
    Function to calculate sequence profiles

    Parameters
    ----------
    consensus_file : str
        path to consensus file (clustal omega fasta output)

    Return
    ----------
    pssm : pssm file
        # Hint check imports
    """
    # Read the alignment file
    alignment = AlignIO.read(clustal_omega_file, "fasta")

    # Create a summary of the alignment
    summary_align = AlignInfo.SummaryInfo(alignment)

    # Calculate the position-specific scoring matrix (pssm)
    pssm = summary_align.pos_specific_score_matrix()

    # Normalize the pssm
    length = len(alignment)
    for position in pssm:
        for key in position:
            position[key] /= length

    return pssm


def sequence_logo_creator(sequence_profile_csv):
    """
    Function to read sequence profile csv file, plot sequence logo and saves as png

    Parameters
    ----------
    sequence_profile_csv : str
        path to sequence profile csv

    Return
    ----------
    None
    """

    # Read the csv file as a pandas dataframe
    df = pd.read_csv(sequence_profile_csv)

    # Drop the 'Position' column 
    profile_data = df.drop(columns='Position')

    # Create the sequence logo
    logo = logomaker.Logo(profile_data)

    # Save the sequence logo as a png
    logo.fig.savefig(f"./{STUDENT_SURNAME}_sequence_logo.png")


def csv_writer(pssm, output_file):
    """
    Function to write probabilities to csv file

    Parameters
    ----------
    pssm : pssm file
        returned pssm from calculate_sequence_profiles

    output_file : str
        path for output file

    Return
    ----------
    None
    """

    # Write the pssm to a csv file
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Position'] + list(pssm[0].keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

        for i, position in enumerate(pssm, start=1):
            row = {'Position': i}
            row.update(position)
            writer.writerow(row)
if __name__ == "__main__":
    if STUDENT_SURNAME == None:
        print(f"Helala")
        sys.exit()

    seq_profiles = calculate_sequence_profiles(sys.argv[1])
    csv_writer(seq_profiles, f"./{STUDENT_SURNAME}_sequence_profile.csv")


    # Bonus point:
    sequence_logo_creator(f"./{STUDENT_SURNAME}_sequence_profile.csv")

