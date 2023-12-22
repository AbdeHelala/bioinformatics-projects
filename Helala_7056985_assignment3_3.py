import sys
import Bio

def calculate_consensus_sequence(clustal_omega_file):
    """
    Function to calculate consensus sequence

    Parameters
    ----------
    clustal_omega_file : str
        path to clustal omega output file

    Return
    ----------
    None
    """
    # Read the alignment file
    alignment = AlignIO.read(clustal_omega_file, "fasta")

    # Create a summary of the alignment
    summary_align = AlignInfo.SummaryInfo(alignment)

    # Calculate the consensus sequence
    consensus = summary_align.dumb_consensus()

    # Print the consensus sequence
    print(consensus)

if __name__ == "__main__":
    calculate_consensus_sequence(sys.argv[1])


