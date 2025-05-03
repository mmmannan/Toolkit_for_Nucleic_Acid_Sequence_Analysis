# Code testing file
from DNA_toolkit import *

import random
# Creating a random DNA sequence for testing
randDNAStr = ''.join([random.choice(DNA_Nucleotides)
                      for nuc in range(50)])

DNAStr = validate_Seq(randDNAStr)

print(f'\nSequence: {DNAStr}')
print(f'[1] + Sequence Length: {len(DNAStr)}')
print(f'[2] + Nucleotide Frequency: {countNucFrequency(DNAStr)}\n')
print(f'[3] + DNA -> RNA Transcription: {transcription(DNAStr)}\n')
print(f"[4] + DNA String + Complement + Reverse Complement:\n5' {DNAStr} 3'")
print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
print(f"3' {complement(DNAStr)} 5'  [Complement]")
print(f"5' {reverse_complement(DNAStr)} 3'  [Reverse Complement]\n")

print(f'\n[5] + GC Content: {gc_content(DNAStr)}%')
print(f'[6] + GC Content of Subsection: {gc_content_subsec(DNAStr)}')
print(f'[6] + GC Content of Subsection with length 5: {gc_content_subsec(DNAStr, 5)}\n')

print(f'\n[7] + Aminoacid Sequence from DNA: {translate_DNA_sequence(DNAStr)}')
print(f'[7] + Aminoacid Sequence from DNA starting at index 1: {translate_DNA_sequence(DNAStr, 1)}')
print(f'[7] + Aminoacid Sequence from DNA starting at index 2: {translate_DNA_sequence(DNAStr, 2)}')
print(f'[7] + Aminoacid Sequence from DNA starting at index 3: {translate_DNA_sequence(DNAStr, 3)}')
print(f'[7] + Aminoacid Sequence from DNA starting at index 4: {translate_DNA_sequence(DNAStr, 4)}')
print(f'[7] + Aminoacid Sequence from DNA starting at index 5: {translate_DNA_sequence(DNAStr, 5)}\n\n')

print(f'[8] + Codon frequency of amino acid (L): {codon_usage(DNAStr, "L")}\n\n')

print(f'[9] + Reading frames:')
for frame in gen_reading_frames(DNAStr):
    print(frame)


test_aa_seq = ['K', 'L', 'V', 'M', 'F', 'G', 'D', 'W', 'E', 'I', 'R', 'A', 'P', 'S', 'T', '_', 'N']   # example
print(f'\n\n{proteins_from_rf(test_aa_seq)}\n')


print('\n[10] + All proteins in 6 open reading frames:')
for prot in all_proteins_from_orfs(DNAStr, 0, 0, True):
    print(f'{prot}')

