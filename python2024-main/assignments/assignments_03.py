def read_config(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    stacks = {}
    moves = []
    
    stack_lines = []
    move_lines = []
    parsing_moves = False

    for line in lines:
        if line.startswith('m'):
            parsing_moves = True
        if parsing_moves:
            move_lines.append(line.strip())
        else:
            stack_lines.append(line.rstrip())

    last_stack_line = stack_lines.pop().strip()
    num_stacks = len(last_stack_line.split())

    for i in range(1, num_stacks + 1):
        stacks[i] = []

    for line in stack_lines:
        parcels = [line[i:i+4].strip() for i in range(0, len(line), 4)]
        for i, parcel in enumerate(parcels):
            if parcel != '':
                stacks[i + 1].insert(0, parcel.strip('[]'))


    for line in move_lines:
        parts = line.split()
        if len(parts) == 6:
            num_parcels = int(parts[1])
            current_stack = int(parts[3])
            new_stack = int(parts[5])
            moves.append((num_parcels, current_stack, new_stack))


    return stacks, moves

def rearrange_parcels(stacks, moves):
    for numParcels, currentStack, newStack in moves:
        for _ in range(numParcels):
            parcel = stacks[currentStack].pop()
            stacks[newStack].append(parcel)
    
    top_parcels = ''.join(stacks[key][-1] for key in sorted(stacks.keys()))
    return top_parcels

# Task 2

def pwm(filename):
    with open(filename, 'r') as file:
        sequences = [line.strip() for line in file.readlines()]
    
    if not sequences:
        return []

    n = len(sequences[0])
    
    pwm = [[0] * n for _ in range(4)]
    
    nucleotide_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    for seq in sequences:
        for i, nucleotide in enumerate(seq):
            pwm[nucleotide_to_index[nucleotide]][i] += 1
    
    return pwm

def consensus_sequence(pwm):
    from itertools import product

    if not pwm:
        return []

    n = len(pwm[0])
    
    index_to_nucleotide = ['A', 'C', 'G', 'T']
    
    consensus_chars = []
    
    for j in range(n):
        counts = [pwm[i][j] for i in range(4)]
        max_count = max(counts)
        
        consensus_bases = [index_to_nucleotide[i] for i in range(4) if counts[i] == max_count]
        
        consensus_chars.append(consensus_bases)
    
    all_consensus_sequences = [''.join(seq) for seq in product(*consensus_chars)]
    
    return sorted(all_consensus_sequences)



#task3
def solve(n, operations):
    current_value = n
    ops = operations.split()
    
    i = 0
    while i < len(ops):
        operation = ops[i]
        value = int(ops[i + 1])
        
        match operation:
            case 'add':
                current_value += value
            case 'mul':
                current_value *= value
            case 'div':
                current_value //= value
            case 'mod':
                current_value %= value
            
        
        i += 2
    
    return current_value