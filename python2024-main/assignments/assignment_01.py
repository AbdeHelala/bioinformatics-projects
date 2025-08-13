''' 
    Task 1: Write a function process_list (2 points)
    Write a function process_list that modifies each
    element in a list and returns the maximum value of the 
    modified list.
    Specification:
    Name: process_list
    Parameter: list of integers
    Return value: largest number in the list after modification
    Modification:
    1- If an elements is odd, multiply the element
      by 2; even, divide the element by 2.
    2- If the position in the list is a multiple of 7,
       add the position to the value.
'''
def process_list(nums):
    if not nums:  
        return None

    mx = nums[0]

    for i, num in enumerate(nums): 
        if num % 2 == 0:  
            nums[i] //= 2
        else:  
            nums[i] *= 2

        if i % 7 == 0:  
            nums[i] += i

        if(mx < nums[i]):
            mx = nums[i]

    return mx


'''
    Task 2: DNA Reverse Complement (2 points)
    Write a function DNA_complement that takes a string as input 
    and computes its reverse, DNA complement (A ↔ T, C ↔ G).
    For example, ACGATCGATCGATTC ↔ GAATCGATCGATCGT.
    You can assume that the input string only contains the 
    upper case letters A,C,G,T.
'''
def DNA_complement(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = ''
    for i in seq[::-1]:
        reverse_seq += complement_dict[i]
    return reverse_seq

'''
    Task 3.1: Extracting k-mers (substrings of length k) from a DNA sequence (3 points)
    1- Write a function list_kmers that, given a DNA sequence s and an integer k,
    returns a list of all k-mers (strings of length k, from left to right, 
    overlapping). Handle edge cases correctly (k > |s| ?!).
    
'''
def list_kmers(seq, k):
    kmers = []
    n = len(seq)
    
    if k > n:
        return kmers  
    
    for i in range(n - k + 1):
        kmers.append(seq[i:i+k])
    
    return kmers
'''
    Task3.2: Write a function number_of_unique_kmers that, given s and k,
    returns the number of unique k-mers in s (k-mers that appear only once).
'''

def number_of_unique(s, k):
    kmer_count = {}
    n = len(s)
    for i in range( n - k + 1):
        kmer = s[i:i+k]
        kmer_count[kmer] = kmer_count.get(kmer, 0) + 1
    unique_count = 0
    for count in kmer_count.values():
        if count == 1:
            unique_count += 1                   
    return unique_count





'''
    Task 4.1: Integer encoding of k-mers (3 points)
    For a given k, there are 4k different possible k-mers.
    Let’s encode every k-mer
    bijectively as an integer 0 ≤ c < 4k .
    Use this encoding: Convert A >→ 0, C >→ 1, G >→ 2, T >→ 3.
    Interpret the sequence of digits as an integer with base 4.
    Example: ACGT >→ (0, 1, 2, 3)^4 = 0·4^3 + 1·4^2 + 2·4^1 + 3·4^0 = 16 + 8 + 3 = 27.
    Write a function kmer_code(kmer) that returns the integer encoding of the string
    kmer using bit operations.
    Your are not allowed to use any packages, like numpy, in this exercise.
'''


def kmer_code(kmer):
    encoding_lst = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    code = 0
    for nucleotide in kmer:
        code <<= 2  
        code += encoding_lst.get(nucleotide)
    return code

'''
    Task4.2: Write a function kmer_decode(code, k) that returns the correct string for a
    given integer encoding and k-mer length k using bit operations.
    Your are not allowed to use any packages, like numpy, in this exercise.
'''

def kmer_decode(code, k):
    decoding_lst = ['A', 'C', 'G', 'T']
    kmer = ''
    for i in range(k):
        index = code & 3  
        kmer = decoding_lst[index] + kmer
        code >>= 2  
    return kmer

