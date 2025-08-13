from assignments_05 import *
import numpy as np
import signal

def test_move_mean_python_1():
    assert move_mean_python([2,4,6,8,12,14,16,18,20], 5) == [6.4,8.8,11.2,13.6,16]

def test_move_mean_numba_1():
    assert (move_mean_numba(np.array([2,4,6,8,12,14,16,18,20]), 5) == np.array([6.4,8.8,11.2,13.6,16])).all()


def handle_timeout(sig, frame):
    raise TimeoutError('took too long')

def test_move_mean_python_2():
    rng = np.random.default_rng()
    lst = list(rng.integers(10**6, size=10**7))
    signal.signal(signal.SIGALRM, handle_timeout)
    signal.alarm(10)
    move_mean_python(lst, 20)
    signal.alarm(0)
    
def test_move_mean_numba_2():
    rng = np.random.default_rng()
    lst = rng.integers(10**6, size=10**7)
    signal.signal(signal.SIGALRM, handle_timeout)
    signal.alarm(1)
    move_mean_numba(lst, 20)
    signal.alarm(0)

def test_largest_k_2():
    assert largest_k(2) == 2

def test_largest_k_4():
    assert largest_k(4) == 3

def test_largest_k_8():
    assert largest_k(8) == 4

def test_largest_k_12():
    assert largest_k(12) == 6