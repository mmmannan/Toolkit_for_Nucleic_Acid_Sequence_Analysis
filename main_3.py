# Nucleic_Acid_toolkit_OOP code testing file

from Nucleic_Acid_toolkit_OOP import bio_seq

# Create an instance test_Str of bio_seq class
# to have all of the functionality from bio_seq class to be accessible through test_Str instance
test_Str = bio_seq()
test_Str.generate_rnd_seq(92, "RNA")
print(test_Str.get_seq_info())
print(f'\nNucleotide Frequency: {test_Str.countNucFrequency()}')
print(f"\nDNA -> RNA Transcription: {test_Str.transcription()}")
print(f"\n5' -> 3' Reverse Complement: {test_Str.reverse_complement()}")
print(f"\nGC Content: {test_Str.gc_content()}%")
print(f'GC Content of Subsection: {test_Str.gc_content_subsec()}')
print(f'GC Content of Subsection with length 10: {test_Str.gc_content_subsec(10)}\n')
print(f"Aminoacid Sequence from DNA: {test_Str.translate_sequence()}")
print(f"Aminoacid Sequence from DNA starting at index 17: {test_Str.translate_sequence(17)}")
print(f'\nCodon frequency of amino acid (R): {test_Str.codon_usage("R")}\n')

print("Reading frames:")
for rf in test_Str.gen_reading_frames():
    print(rf)

print(f"\n{test_Str.proteins_from_rf(
    ['C', 'G', 'R', 'M', 'E', 'R', 'S', 'W', 'L', 'G', '_', 'K', 'T'])}")   # example

print(f"\nAll proteins in 6 open reading frames:\n{test_Str.all_proteins_from_orfs()}")

