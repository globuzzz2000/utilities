from Bio import SeqIO
import re

# Configuration
fasta_file = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Arenosa2.fna"
min_len = 100
max_chromosomes = 8  # only process first 8 sequences

print("Chromosome\tStart_Pos\tLength")

# Parse sequences
for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
    if i >= max_chromosomes:
        break  # Skip any chromosomes beyond the first 8

    seq = str(record.seq).upper()
    chrom = record.id
    print(f"Processing {chrom} ({i + 1}) with length {len(seq)} bases...")

    # Find all N-stretches >= min_len
    stretches = [(m.start(), len(m.group())) for m in re.finditer(r'N{' + str(min_len) + ',}', seq)]
    if not stretches:
        continue

    for start, length in stretches:
        print(f"{chrom}\t{start}\t{length}")