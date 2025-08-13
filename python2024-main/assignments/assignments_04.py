import numpy as np

class CountingBloomFilter:
    
    def __init__(self, size, num_hashes, hash_factors):
        self.size = size
        self.num_hashes = num_hashes
        self.hash_factors = hash_factors
        self.counters = np.zeros(size, dtype=np.uint8)

    def _hash(self, item, seed):
        return (item * self.hash_factors[seed]) % self.size

    def insert(self, item):
        for seed in range(self.num_hashes):
            index = self._hash(item, seed)
            self.counters[index] += 1

    def delete(self, item):
        for seed in range(self.num_hashes):
            index = self._hash(item, seed)
            if self.counters[index] > 0:
                self.counters[index] -= 1

    def contains(self, item):
        for seed in range(self.num_hashes):
            index = self._hash(item, seed)
            if self.counters[index] == 0:
                return False
        return True

    def count(self, item):
        min_count = float('inf')
        for seed in range(self.num_hashes):
            index = self._hash(item, seed)
            min_count = min(min_count, self.counters[index])
        return min_count if min_count != float('inf') else 0


def get_treasure(filename):
    with open(filename, 'r') as file:
        numbers = [int(line.strip()) for line in file]

    n = len(numbers)
    indexed_numbers = list(enumerate(numbers))
    
    def move_element(arr, index, value):
        elem = arr.pop(index)
        new_index = (index + value) % len(arr)
        arr.insert(new_index, elem)

    for original_index, value in enumerate(numbers):
        current_index = next(i for i, (idx, val) in enumerate(indexed_numbers) if idx == original_index)
        move_element(indexed_numbers, current_index, value)

    zero_index = next(i for i, (idx, val) in enumerate(indexed_numbers) if val == 0)

    result = 0
    for offset in [1000, 2000, 3000]:
        result += indexed_numbers[(zero_index + offset) % n][1]

    return result


def quantile_normalize(matrix):
    sorted_matrix = np.sort(matrix, axis=0)
    
    rank_means = np.mean(sorted_matrix, axis=1)
    
    ranked_matrix = np.argsort(np.argsort(matrix, axis=0), axis=0)
    
    quantile_normalized_matrix = np.zeros_like(matrix)
    for col in range(matrix.shape[1]):
        for row in range(matrix.shape[0]):
            quantile_normalized_matrix[row, col] = rank_means[ranked_matrix[row, col]]
    
    return quantile_normalized_matrix    
