from assignments_02 import *

def test_anagrams():
    assert anagrams(['no', 'peels', 'sleep', 'on', 'leeps', 'fruit']) == {'peels', 'sleep', 'leeps'}

def test_hamming_distance():
    assert hamming_distance('ass_02.2.txt') == (1,2)

def test_common_kmers():
    assert common_kmers('ass_02.3.fasta', 8) == {'TGGGCTCA', 'GGCTGGTC', 'CACTGTGC', 'TCAGCATC', 'CTGGCTGC', 'CAGCCTCA'}

def test_no_common_kmers():
    assert common_kmers('ass_02.3.fasta', 21) == set()

def test_vampire_numbers():
    assert is_vampire(1260)
    assert not is_vampire(126000)

def test_vampire_generator():
    VG = vampire_generator()
    assert [next(VG) for _ in range(5)] == [1260, 1395, 1435, 1530, 1827]
