from assignments_04 import *
import numpy as np

def test_counting_bloom_filter():
    rng = np.random.default_rng()
    numbers = rng.integers(low=0, high=10000, size=100) 

    filter = CountingBloomFilter(100, 3, [6441, 1117, 7577])
    for n in numbers:
        assert filter.contains(n) == False
        assert filter.count(n) == 0

    for n in numbers:
        filter.insert(n)
    
    for n in numbers:
        assert filter.contains(n) == True


def test_counting_bloom_filter_2():
    filter = CountingBloomFilter(100, 2, [6441, 1117])
    numbers = [10, 20, 65, 10, 10, 99, 72, 99]
    for n in numbers:
        filter.insert(n)
    
    assert filter.count(10) == 3
    assert filter.count(20) == 1
    assert filter.count(65) == 1
    assert filter.count(99) == 2
    assert filter.count(72) == 1

    filter.delete(10)
    assert filter.count(10) == 2

    filter.delete(72)
    assert filter.count(72) == 0


def test_treasure_hunt():
    assert get_treasure('ass_04.2.txt') == 3

def test_quantile_normalize():
    assert (quantile_normalize(np.array([[2,4,4,5], [5,14,4,7], [4,8,6,9], [3,8,5,8], [3,9,3,5]], dtype=np.float64)) == np.array([[3.5,3.5,5,3.5],[8.5,8.5,5.5,5.5],[6.5,5,8.5,8.5],[5,5.5,6.5,6.5],[5.5,6.5,3.5,5]])).all()