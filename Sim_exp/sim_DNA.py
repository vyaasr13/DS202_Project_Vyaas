# Creating the simulated DNA sequence
# Rules: Consists of introns, exons, and intergenic regions
# exons are 100-200 bp long, introns are 1000-2000 bp long, intergenic regions are 5000-10000 bp long
# The sequence is 100000 bp long
# The sequence starts with an intergenic region and ends with an intergenic region
# Then gene should start with an exon and end with an exon
# the gene should be followed by an intergenic region
# The gene shoulf start with a start codon and end with a stop codon
# The start codon is ATG and the stop codon is TAA, TAG, or TGA
# The sequence should be in the format of a FASTA file

# Create a function to generate a random DNA sequence of a given length
import random
import string
import sys
import os
import time
import argparse
import re

# in an exon, the probability of all codons is equal
# in an intron, A 0.4, C 0.2, G 0.2, T 0.4
# in intergenic region, A 0.3, C 0.2, G 0.2, T 0.3

# store all the possible codons in a list
codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT',
          'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATT',
          'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
          'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
          'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
          'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
          'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
          'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG']

def generate_dna_sequence(region_type):
    if region_type == 'exon':
        length = random.randint(100, 200)
    elif region_type == 'intron':
        length = random.randint(1000, 2000)
    elif region_type == 'intergenic':
        length = random.randint(5000, 10000)

    if region_type == 'exon':
        sequence = ''
        while len(sequence) < length:   
            if len(sequence) < 3:
                base = random.choices(['A', 'C', 'G', 'T'], weights=[0.25, 0.25, 0.25, 0.25], k=1)[0]
                sequence += base
            # check if the last 3 bases are a stop codon or start codon
            if len(sequence) >= 3:
                last_codon = sequence[-3:]
                if last_codon in {'ATG', 'TAA', 'TAG', 'TGA'}:
                    # remove the last 1 bases
                    sequence = sequence[:-1]
                # add a random base to the end of the sequence
                sequence += random.choices(['A', 'C', 'G', 'T'], weights=[0.25, 0.25, 0.25, 0.25], k=1)[0]
        return sequence
    elif region_type == 'intron':
        sequence = ''
        while len(sequence) < length:   
            if len(sequence) < 3:
                base = random.choices(['A', 'C', 'G', 'T'], weights=[0.4, 0.1, 0.1, 0.4], k=1)[0]
                sequence += base
            # check if the last 3 bases are a stop codon or start codon
            if len(sequence) >= 3:
                last_codon = sequence[-3:]
                if last_codon in {'ATG', 'TAA', 'TAG', 'TGA'}:
                    # remove the last 3 bases
                    sequence = sequence[:-3]
                # add a random base to the end of the sequence
                sequence += random.choices(['A', 'C', 'G', 'T'], weights=[0.4, 0.1, 0.1, 0.4], k=1)[0]
        return sequence
    elif region_type == 'intergenic':
        sequence = ''
        while len(sequence) < length:
            if len(sequence) < 3:
                base = random.choices(['A', 'C', 'G', 'T'], weights=[0.3, 0.2, 0.2, 0.3], k=1)[0]
                sequence += base
            # check if the last 3 bases are a stop codon or start codon
            if len(sequence) >= 3:
                last_codon = sequence[-3:]
                if last_codon in {'ATG', 'TAA', 'TAG', 'TGA'}:
                    # remove the last 3 bases
                    sequence = sequence[:-3]
                # add a random base to the end of the sequence
                sequence += random.choices(['A', 'C', 'G', 'T'], weights=[0.3, 0.2, 0.2, 0.3], k=1)[0]
        return sequence
    else:
        raise ValueError("Invalid region type. Must be 'exon', 'intron', or 'intergenic'.")

# Modified gene generator to track exon positions
def generate_gene(start_pos):
    gene = ''
    exon_positions = []
    num_exons = random.randint(2, 5)
    current_pos = start_pos
    print(f"current_pos: {current_pos}")

    # First exon
    exon_seq = generate_dna_sequence("exon")
    gene += 'ATG' + exon_seq
    exon_start = current_pos
    exon_end = exon_start + len('ATG') + len(exon_seq) - 1
    exon_positions.append((exon_start, exon_end))
    current_pos = exon_end + 1

    for i in range(num_exons):
        # Intron
        if i < num_exons:
            # Generate intron sequence
            intron_seq = generate_dna_sequence("intron")
            gene += intron_seq
            current_pos += len(intron_seq)

        # Next exon
        exon_seq = generate_dna_sequence("exon")
        exon_start = current_pos
        gene += exon_seq
        exon_end = exon_start + len(exon_seq) - 1
        exon_positions.append((exon_start, exon_end))
        current_pos = exon_end + 1

    # Append stop codon
    stop_codon = random.choice(['TAA', 'TAG', 'TGA'])
    gene += stop_codon
    current_pos += len(stop_codon)
    print(f"current_pos: {current_pos}")

    return gene, exon_positions, current_pos

def generate_dna_sequence_with_genes(total_length):
    sequence = ''
    current_length = 0
    all_exon_positions = []

    # Start with an intergenic region
    intergenic_seq = generate_dna_sequence('intergenic')
    sequence += intergenic_seq
    current_length += len(intergenic_seq)

    num_genes = 0

    while current_length < total_length:
        num_genes += 1
        gene_start = current_length
        gene, exon_positions, gene_end = generate_gene(gene_start)
        sequence += gene
        all_exon_positions.extend(exon_positions)
        current_length = gene_end

        if current_length < total_length:
            intergenic_seq = generate_dna_sequence('intergenic')
            sequence += intergenic_seq
            current_length += len(intergenic_seq)

    return sequence, all_exon_positions, num_genes

def write_fasta_file(sequence, filename):
    with open(filename, 'w') as fasta_file:
        fasta_file.write('>simulated_dna_sequence\n')
        for i in range(0, len(sequence), 80):
            fasta_file.write(sequence[i:i+80] + '\n')

def write_exon_positions(exon_positions, filename):
    with open(filename, 'w') as f:
        f.write('Start\tEnd\n')
        for start, end in exon_positions:
            f.write(f'{start}\t{end}\n')

dna_sequence, exon_positions, num_genes = generate_dna_sequence_with_genes(100000)
write_fasta_file(dna_sequence, 'simulated_dna_sequence.fasta')
write_exon_positions(exon_positions, 'exon_positions.txt')
print("Simulated DNA sequence and exon positions have been written to files.") 
print(f"number of genes: {num_genes}") 

