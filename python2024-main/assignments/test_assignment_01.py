from assignment_01 import *

# Task 1
def test_process_list_1():
	assert process_list([2, 3, 4, 4, 7, 7, 8, 9, 10, 1, 1]) == 25

# Task 2
def test_DNA_complement():
	assert DNA_complement("ACGATCGATCGATTC") == "GAATCGATCGATCGT"
	

# Task 3.1
def test_list_kmers():
	assert list_kmers("ACGATCGATC", 5) == ["ACGAT", "CGATC", "GATCG", "ATCGA", "TCGAT", "CGATC"]

# Task 3.2
def test_unique_kmers():
	assert number_of_unique("ACGATCGATC", 5) == 4

# Task 4.1	
def test_integer_encoding_ACTG():
	assert kmer_code("ACTG") == 30
	
# Task 4.2
def test_integer_decoding_ACTG():
	assert kmer_decode(30, 4) == "ACTG"