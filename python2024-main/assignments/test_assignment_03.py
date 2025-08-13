from assignments_03 import *

def test_read_config():
    assert read_config('ass_03.1.txt') == (dict({1:['A', 'B'], 2:['C', 'D', 'E'], 3: ['F']}), [(1,1,3), (2,2,1)])

def test_rearrange_parcels():
    assert rearrange_parcels(dict({1:['A', 'B'], 2:['C', 'D', 'E'], 3: ['F']}), [(1,1,3), (2,2,1)]) == 'DCB'

def test_pwm():
    assert pwm('ass_03.2.txt') == [[2,0,0,0,0,3,0], [0,4,3,0,1,0,4], [0,0,1,0,3,0,0], [2,0,0,4,0,1,0]]

def test_consensus_sequence():
    assert consensus_sequence([[1,0,0,0,0,3,0], [0,4,3,0,1,0,4], [0,0,1,0,3,0,0], [3,0,0,4,0,1,0]]) == ['TCCTGAC']

def test_consensus_sequence():
    assert consensus_sequence([[2,0,0,0,0,3,0], [0,4,3,0,1,0,4], [0,0,1,0,3,0,0], [2,0,0,4,0,1,0]]) == ['ACCTGAC', 'TCCTGAC']

def test_solve():
    assert solve(10, "add 6 mul 2 div 4 mod 5 mul 3") == 9