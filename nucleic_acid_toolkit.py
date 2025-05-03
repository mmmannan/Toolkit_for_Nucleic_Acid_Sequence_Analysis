from bio_structures import *
from collections import Counter


## docstring '''  ''' is describes what a function, class, or module does.

# Check the sequence to make sure it is a DNA string
def validate_Seq(dna_seq):
    '''Check the sequence to make sure it is a valid DNA string.'''
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in DNA_Nucleotides:
            return False
    return tmpseq


# Counting nucleotides
def countNucFrequency(seq):
    '''Counting nucleotides in a given sequence. Return a dictionary.'''
    tmpFreqDict = {"A":0, "C":0, "G":0, "T":0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict
### OR,
# import collections
# def countNucFrequency(seq):
#    return dict(collections.Counter(seq))
#### OR,
# def countNucFrequency(seq):
#    return dict(Counter(seq))


# Transcription of coding DNA strand
def transcription(seq):
    '''Transcribe a DNA strand into RNA sequence.'''
    return seq.replace("T", "U")


DNA_ReverseComplement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
# Reverse complement
def reverse_complement(seq):
    '''Return 5' to 3' DNA reverse complement strand.'''
    return ''.join([DNA_ReverseComplement[nuc] 
                    for nuc in seq])[::-1]

# Pythonic approach. A little bit faster solution.
# def reverse_complement(seq):
#    mapping = str.maketrans('ATCG', 'TAGC')
#    return seq.translate(mapping)[::-1]


DNA_Complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
# Complement, not reverse
def complement(seq):
    '''Return 3' to 5' DNA complement strand'''
    return ''.join([DNA_Complement[nuc]
                    for nuc in seq])

# {reverse_complement(randDNAStr)}[::-1]


# GC content in a DNA or RNA sequence
def gc_content(seq):
    '''Return GC content in a DNA or RNA sequence.'''
    return round(((seq.count('C') + seq.count('G')) / len(seq)) * 100, 6)


# GC Content in a DNA or RNA subsequence of length k. k=20 by default.
# Window sliding where window size = k, not overlapping
def gc_content_subsec(seq, k=20):
    '''Calculate GC Content in a DNA or RNA subsequence of length k. k=20 by default. Return a list.'''
    res = []
    for i in range(0, len(seq)-k+1, k):         # jump or gap = k = window size
        subseq = seq[i:i+k]
        res.append(gc_content(subseq))
    return res



# Translation of a coding DNA sequence into amino acid sequence
def translate_DNA_sequence(seq, init_pos = 0):
    '''Translate a coding DNA sequence into amino acid sequence.'''
    return [DNA_Codons[seq[pos : pos + 3]]
            for pos in range(init_pos, len(seq) - 2, 3)]


# Provides the frequency of each codon encoding a given amino acid in a DNA sequence
def codon_usage(seq, aminoacid):
    '''Provides the frequency of each codon encoding a given amino acid in a DNA sequence.'''
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i : i+3]] == aminoacid:
            tmpList.append(seq[i : i+3])

    freqDict = dict(Counter(tmpList))
    totalWight = sum(freqDict.values())
    for codon in freqDict:
        freqDict[codon] = round(freqDict[codon] / totalWight, 2)
    return freqDict


# Generate the 6 reading frames of a DNA sequence, including reverse complement
def gen_reading_frames(seq):
    '''Generate the 6 reading frames of a DNA sequence, including reverse complement'''
    frames = []
    frames.append(translate_DNA_sequence(seq, 0))
    frames.append(translate_DNA_sequence(seq, 1))
    frames.append(translate_DNA_sequence(seq, 2))
    frames.append(translate_DNA_sequence(reverse_complement(seq), 0))
    frames.append(translate_DNA_sequence(reverse_complement(seq), 1))
    frames.append(translate_DNA_sequence(reverse_complement(seq), 2))
    return frames


# Compute all possible proteins in an amino acid sequence and return a list of possible proteins
def proteins_from_rf(aa_seq):
    '''Compute all possible proteins in an amino acid sequence and return a list of possible proteins'''
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            # STOP accumulating amino acids if _ -STOP was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            # START accumulating amino acids if M - START was found
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins


def all_proteins_from_orfs(seq, startReadPos=0, endReadPos=0, ordered=False):
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startRead : endRead])
    else:
        rfs = gen_reading_frames(seq)

    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)
    
    if ordered:
        return sorted(res, key=len, reverse=True)
    
    return res


