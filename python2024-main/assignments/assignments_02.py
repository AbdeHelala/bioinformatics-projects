import itertools
from collections import defaultdict
# Task 1


def anagrams(words):
    anagrams = defaultdict(set)

    if not words:
        return defaultdict(set)
    
    for word in words:
        sortedword = ''.join(sorted(word))
        anagrams[sortedword].add(word)
    
    largestsubset = max(anagrams.values(), key=len)
    
    return largestsubset




# Task 2
def hamming_distance(filename):
 
    try:
        with open(filename, 'r') as file:
            sequences = [line.strip().replace('U', 'T') for line in file]
    except FileNotFoundError:
        return -1, -1  
    
    length = len(sequences[0])
    if not all(len(seq) == length for seq in sequences):
        return -1, -1  

    def count_hamming_distance(str1, str2):
        return sum(c1.upper() != c2.upper() for c1, c2 in zip(str1, str2))

    mn = length  
    mx = 0       

    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            distance = count_hamming_distance(sequences[i], sequences[j])
            mn = min(mn, distance)
            mx = max(mx, distance)

    return mn, mx








# Task 3
def common_kmers(fasta, k):
    
    def extract_sequences(file):
       
        sequences = []
        with open(file, 'r') as file:
            sequence = ''
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if sequence:
                        sequences.append(sequence)
                        sequence = ''
                else:
                    sequence += line
            if sequence:
                sequences.append(sequence)
        return sequences

    def generate_kmers(sequence, k):
        
        kmers = set()
        n = len(sequence)
        for i in range( n - k + 1):
            kmers.add(sequence[i:i+k])
        return kmers

    sequences = extract_sequences(fasta)
    if not sequences:
        return set()  

    common = generate_kmers(sequences[0], k)

    for seq in sequences[1:]:
        kmers = generate_kmers(seq, k)
        common.intersection_update(kmers)

    return common


# Task 4
def is_vampire(n):

    number_of_digits = [int(d) for d in str(n)]
    if len(number_of_digits) % 2 != 0:
        return False
    
    fangs = itertools.permutations(number_of_digits, len(number_of_digits)//2)
    
    for fang_pairs in fangs:
        first_fang = int(''.join(map(str, fang_pairs)))
        second_fang_digits = [digit for digit in number_of_digits if digit not in fang_pairs]
        if first_fang != 0 and second_fang_digits:
            second_fang = int(''.join(map(str, second_fang_digits)))
            if first_fang * second_fang == n:
                return True
    return False

def vampire_generator(n = 1259):
    
    while True:
        if is_vampire(n):
            yield n
        n += 1