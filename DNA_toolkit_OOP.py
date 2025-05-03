from bio_structures import DNA_Nucleotides, DNA_Codons
from collections import Counter
import random


class bio_seq:
    '''DNA sequence class. Default value: ATCG, DNA, No label.'''

    def __init__(self, seq="ATCG", seq_type="DNA", label='No Label'):
        '''Sequence initialization and validation.'''
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.validate()
        assert self.is_valid, f"Provided data does not seem to be correct {self.seq_type} sequence"


    def validate(self):
        '''Check the sequence to make sure it is a valid DNA string.'''
        return set(DNA_Nucleotides).issuperset(self.seq)


    def get_seq_biotype(self):
        '''Return sequence type.'''
        return self.seq_type


    def get_seq_info(self):
        '''Return 4 strings. Full sequence information.'''
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"


    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        '''Generate a random DNA sequence, provided the length.'''
        seq = ''.join([random.choice(DNA_Nucleotides)
                       for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")


    def countNucFrequency(self):
        '''Counting nucleotides in a given sequence. Return a dictionary.'''
        return dict(Counter(self.seq))


    def transcription(self):
        '''Transcribe a DNA strand into RNA sequence.'''
        return self.seq.replace("T", "U")


    def reverse_complement(self):
        '''Return 5' to 3' DNA reverse complement strand.'''
        mapping = str.maketrans('ATCG', 'TAGC')
        return self.seq.translate(mapping)[::-1]


    def gc_content(self):
        '''Return GC content of a DNA sequence.'''
        return round(((self.seq.count('C') + self.seq.count('G')) / len(self.seq)) * 100, 2)


    def gc_content_subsec(self, k=20):
        '''Calculate GC Content of a DNA subsequence of length k. k=20 by default. Return a list.'''
        res = []
        for i in range(0, len(self.seq)-k+1, k):         # jump or gap = k = window size
            subseq = self.seq[i:i+k]
            res.append(
                round(((subseq.count('C') + subseq.count('G')) / len(subseq)) * 100))
        return res


    def translate_DNA_sequence(self, init_pos = 0):
        '''Translate a coding DNA sequence into amino acid sequence.'''
        return [DNA_Codons[self.seq[pos : pos + 3]]
                for pos in range(init_pos, len(self.seq) - 2, 3)]


    def codon_usage(self, aminoacid):
        '''Provide the frequency of each codon encoding a given amino acid in a DNA sequence.'''
        tmpList = []
        for i in range(0, len(self.seq) - 2, 3):
            if DNA_Codons[self.seq[i : i+3]] == aminoacid:
                tmpList.append(self.seq[i : i+3])

        freqDict = dict(Counter(tmpList))
        totalWight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWight, 2)
        return freqDict


    def gen_reading_frames(self):
        '''Generate the 6 reading frames of a DNA sequence, including reverse complement.'''
        frames = []
        frames.append(self.translate_DNA_sequence(0))
        frames.append(self.translate_DNA_sequence(1))
        frames.append(self.translate_DNA_sequence(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_DNA_sequence(0))
        frames.append(tmp_seq.translate_DNA_sequence(1))
        frames.append(tmp_seq.translate_DNA_sequence(2))
        del tmp_seq
        return frames


    def proteins_from_rf(self, aa_seq):
        '''Compute all possible proteins in an amino acid sequence and return a list of possible proteins.'''
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


    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(
                self.seq[startReadPos : endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()
        
        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)
        
        if ordered:
            return sorted(res, key=len, reverse=True)
        
        return res


