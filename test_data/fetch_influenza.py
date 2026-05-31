from Bio import Entrez, SeqIO

# Always tell NCBI who you are (use your student email)
Entrez.email = "04282213040@student.qau.edu.pk"

# The 8 segments of Influenza A (H1N1) PR8 strain
influenza_segments = [
    "NC_007373.1", "NC_007372.1", "NC_007371.1", "NC_007366.1", 
    "NC_007369.1", "NC_007368.1", "NC_007370.1", "NC_007367.1"
]

print("Downloading 8 Influenza segments from NCBI...")
records = []

# Fetch each segment from the NCBI API
for acc in influenza_segments:
    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    records.append(record)
    print(f"Downloaded: {record.description}")

from pathlib import Path

# Write all 8 segments into ONE single FASTA file
output_file = Path(__file__).parent / "04_segmented_virus" / "test_04_influenza_A.fasta"
SeqIO.write(records, str(output_file), "fasta")
print(f"\nSuccess! Saved all segments to {output_file}")