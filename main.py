# DNA Toolset/Code testing file
from nucleic_acid_toolkit import *

import random
# Creating a random DNA sequence for testing:
randDNAStr = ''.join([random.choice(DNA_Nucleotides)
                      for nuc in range(50)])

DNAStr = validate_Seq(randDNAStr)

print(f'\nSequence: {DNAStr}\n')
print(f'[1] + Sequence Length: {len(DNAStr)}\n')
print(f'[2] + Nucleotide Frequency: {countNucFrequency(DNAStr)}\n')
print(f'[3] + DNA/RNA Transcription: {transcription(DNAStr)}\n')
print(f"[4] + DNA String + Complement + Reverse Complement:\n5' {DNAStr} 3'")
print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
print(f"3' {complement(DNAStr)} 5'  [Complement]\n")
print(f"5' {reverse_complement(DNAStr)} 3'  [Reverse Complement]\n")

print(f'[5] + GC Content: {gc_content(DNAStr)}%\n')
print(f'[6] + GC Content in Subsection: {gc_content_subsec(DNAStr)}\n')
print(f'[6] + GC Content in Subsection k=5: {gc_content_subsec(DNAStr, 5)}\n')

print(f'[7] + Aminoacids Sequence from DNA: {translate_DNA_sequence(DNAStr)}\n')
print(f'[7] + Aminoacids Sequence from DNA: {translate_DNA_sequence(DNAStr, 1)}\n')
print(f'[7] + Aminoacids Sequence from DNA: {translate_DNA_sequence(DNAStr, 2)}\n')
print(f'[7] + Aminoacids Sequence from DNA: {translate_DNA_sequence(DNAStr, 3)}\n')
print(f'[7] + Aminoacids Sequence from DNA: {translate_DNA_sequence(DNAStr, 4)}\n')
print(f'[7] + Aminoacids Sequence from DNA: {translate_DNA_sequence(DNAStr, 5)}\n')

print(f'[8] + Codon frequency of amino acid (L): {codon_usage(DNAStr, "L")}\n')

print(f'[9] + Reading frames:')
for frame in gen_reading_frames(DNAStr):
    print(frame)


test_aa_seq = ['K', 'L', 'V', 'M', 'F', 'G', 'D', 'W', 'E', 'I', 'R', 'A', 'P', 'S', 'T', '_', 'N']
print(f'\n\n{proteins_from_rf(test_aa_seq)}\n\n')


print('\n[10] + All prots in 6 open reading frames:')
for prot in all_proteins_from_orfs(DNAStr, 0, 0, True):
    print(f'{prot}')

