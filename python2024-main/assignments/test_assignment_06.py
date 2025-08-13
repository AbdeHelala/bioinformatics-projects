from assignments_06 import *

def test_possible_sequences():
    assert possible_sequences(9, 4) == 40

def test_sequence_complexity():
    assert sequence_complexity('ATTTGGATT', 4) == 0.875

def test_align_score():
    assert align('GCATGCT', 'GATTACA', (1, -1, -1))[0] == 0

def test_align_bonus():
    assert align('GCATGCT', 'GATTACA', (1, -1, -1))[1] == [('GCA-TGCT', 'G-ATTACA'), ('GCAT-GCT', 'G-ATTACA'), ('GCATG-CT', 'G-ATTACA')]

def test_convert_read_1():
    assert convert_read('test.fastq', 0.0005) == "CAAAGAGAGAAAGAAAAGTCAATGATTTTATAGCCAGGCAAAATGACTTTCAAGTAAAAAATATAAAGCACCTTACAAACTAGTATCAAAATGCATTTCT"

def test_convert_read_2():
    assert convert_read('test.fastq', 0.0001) == "NNNNNNNNNNNNNAAAAGNCNATGATTTTATAGCCAGGCAAAATGACTTTCAAGTAAAAAATATAAAGCACNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"