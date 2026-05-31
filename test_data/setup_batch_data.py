import os
from pathlib import Path
from Bio import Entrez, SeqIO

# Always identify yourself to NCBI
Entrez.email = "04282213040@student.qau.edu.pk"

# Create/target the pre-existing large batch directory
batch_dir = Path(__file__).parent / "20_large_batch_dataset"
os.makedirs(batch_dir, exist_ok=True)

print("Searching NCBI for Dengue Virus Genomes...")
# Search NCBI nucleotide database for Dengue Virus (taxid:12637)
handle = Entrez.esearch(db="nucleotide", term="txid12637[Organism]", retmax=20)
record = Entrez.read(handle)
id_list = record["IdList"]

print(f"Found {len(id_list)} genomes. Downloading...")

# Download and save each genome as a separate FASTA file
for i, seq_id in enumerate(id_list):
    fetch_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
    seq_record = SeqIO.read(fetch_handle, "fasta")
    
    # Save to file
    file_path = os.path.join(batch_dir, f"dengue_genome_{i+1}_{seq_id}.fasta")
    SeqIO.write(seq_record, file_path, "fasta")
    print(f"Saved: {file_path}")

print("Batch data generation complete!")