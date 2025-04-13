# Viterbi algorithm for exon-intron prediction in DNA sequences
# The emission probabilities are known
# The transition probabilities are known
# Find genes in the sequence and find the exons in the genes
# The sequence is in a fasta file
# The output is a list of exon positions in the sequence

import random
import re
from math import log
import os
import sys

# Function to read the sequence from a fasta file
def read_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

# Viterbi algorithm to find exon-intron boundaries
def viterbi_algorithm(sequence):
    states = ['exon', 'intron']
    emission_probs = {
        'exon': {base: 0.25 for base in 'ACGT'},
        'intron': {'A': 0.4, 'C': 0.1, 'G': 0.1, 'T': 0.4}
    }
    # use the log of the probabilities to avoid underflow
    log_emission_probs = {
        state: {base: log(prob) for base, prob in probs.items()}
        for state, probs in emission_probs.items()
    }
    transition_probs = {
        'exon': {'exon': 0.9933, 'intron': 0.0067},
        'intron': {'exon': 0.0067, 'intron': 0.9933}
    }
    # use the log of the probabilities to avoid underflow
    log_transition_probs = {
        state: {next_state: log(prob) for next_state, prob in next_states.items()}
        for state, next_states in transition_probs.items()
    }
    start_probs = {'exon': 1.0, 'intron': 0.0}

    # initialize the Viterbi table and backpointer
    n = len(sequence)
    viterbi_table = [[float('-inf')] * len(states) for _ in range(n)]
    backpointer = [[0] * len(states) for _ in range(n)]
    # initialize the first column of the Viterbi table
    for s, state in enumerate(states):
        if start_probs[state] > 0:
            viterbi_table[0][s] = log(start_probs[state]) + log_emission_probs[state].get(sequence[0], float('-inf'))
        else:
            viterbi_table[0][s] = float('-inf')
        backpointer[0][s] = 0

    for t in range(1, n):
        for s, state in enumerate(states):
            max_prob = float('-inf')
            max_state = 0
            # best path to all the states at the given time t
            for s_prev, prev_state in enumerate(states):
                prob = viterbi_table[t - 1][s_prev] + log_transition_probs[prev_state].get(state, float('-inf'))
                if prob > max_prob:
                    max_prob = prob
                    max_state = s_prev
            # given the state what is the probability of observing the current base
            viterbi_table[t][s] = max_prob + log_emission_probs[state].get(sequence[t], float('-inf'))
            backpointer[t][s] = max_state

    best_path = []
    best_last_state = max(range(len(states)), key=lambda s: viterbi_table[n - 1][s])
    best_path.append(best_last_state)
    for t in range(n - 1, 0, -1):
        best_path.append(backpointer[t][best_path[-1]])
    best_path.reverse()

    best_path_states = [states[state] for state in best_path]

    # Extract exon positions from the best path
    exon_positions = []
    current_exon_start = None
    for i, state in enumerate(best_path_states):
        if state == 'exon':
            if current_exon_start is None:
                current_exon_start = i
        else:
            if current_exon_start is not None:
                exon_positions.append((current_exon_start, i - 1))
                current_exon_start = None
    if current_exon_start is not None:
        exon_positions.append((current_exon_start, n - 1))

    return exon_positions

# Function to find all exons in the sequence

# Function to find all gene regions in the sequence
def find_all_genes(sequence):
    genes = []
    i = 0
    while i < len(sequence) - 2:
        if sequence[i:i+3] == 'ATG':
            for j in range(i+3, len(sequence)-2, 3):
                codon = sequence[j:j+3]
                if codon in ['TAA', 'TAG', 'TGA']:
                    genes.append((i, j+3))  # include stop codon
                    i = j + 3  # move past this gene
                    break
            else:
                i += 3  # no stop codon found, keep looking
        else:
            i += 1
    return genes

# Modified function to find exons in all gene regions
def find_exons_in_sequence(sequence):
    gene_regions = find_all_genes(sequence)
    if not gene_regions:
        print("No genes found in the sequence.")
        return []

    all_exon_positions = []
    for idx, (gene_start, gene_end) in enumerate(gene_regions):
        gene_sequence = sequence[gene_start:gene_end]
        print(f"Gene {idx+1} found from {gene_start} to {gene_end}, length = {gene_end - gene_start}")
        exon_positions = viterbi_algorithm(gene_sequence)
        # Reject exons that are < 100 and > 200 bp
        exon_positions = [(start, end) for start, end in exon_positions if 100 <= (end - start + 1) <= 200]
        exon_positions_in_sequence = [
            (gene_start + start, gene_start + end)
            for start, end in exon_positions
        ]
        all_exon_positions.extend(exon_positions_in_sequence)

    return all_exon_positions

# Function to write the exon positions to a file
def write_exons_to_file(exon_positions, output_file):
    with open(output_file, 'w') as file:
        file.write('Start\tEnd\n')
        for start, end in exon_positions:
            file.write(f"{start}\t{end}\n")
    print(f"Exons written to {output_file}")

# === Main Driver Code ===
if __name__ == "__main__":
    fasta_file = 'simulated_dna_sequence.fasta'
    output_file = 'exons.txt'
    sequence = read_fasta(fasta_file)
    exon_positions = find_exons_in_sequence(sequence)
    write_exons_to_file(exon_positions, output_file)
