#!/usr/bin/env python
"""
setup_test_data.py

This script initializes the directory structure for validating the "PathoScope AI"
viral functional genomics pipeline. It automatically creates 10 test folders and
populates them with synthetic genomics data designed to test various biological
edge cases, standard inputs, quality control thresholds, and error handlers.

Author: Senior Bioinformatics QA Engineer
Date: May 2026
"""

import os
import sys
import gzip
from pathlib import Path

def create_directory_structure(base_dir: Path) -> dict:
    """
    Creates the 10 standard directories required for functional validation.
    
    Args:
        base_dir (Path): The root directory where test_data will be created.
        
    Returns:
        dict: A mapping of folder index/name to its absolute path.
    """
    folders = [
        "01_small_bacteriophage",
        "02_ssRNA_virus",
        "03_dsDNA_virus",
        "04_segmented_virus",
        "05_large_virus",
        "06_fastq_raw_reads",
        "07_invalid_input",
        "08_ambiguous_sequences",
        "09_duplicate_sequences",
        "20_large_batch_dataset"
    ]
    
    paths = {}
    print("=== Initializing PathoScope AI Test Directory Structure ===")
    
    for folder in folders:
        folder_path = base_dir / folder
        try:
            folder_path.mkdir(parents=True, exist_ok=True)
            paths[folder] = folder_path
            print(f"[SUCCESS] Created directory: {folder_path.relative_to(base_dir.parent)}")
            
            # Place a .keep file in folders that are not explicitly populated with synthetic data
            # to preserve directory structures in version control systems.
            keep_file = folder_path / ".keep"
            if not keep_file.exists():
                keep_file.touch()
                
        except Exception as e:
            print(f"[ERROR] Failed to create directory {folder_path}: {e}", file=sys.stderr)
            
    return paths

def generate_valid_phage_fasta(folder_path: Path):
    """
    Generates a valid 500bp synthetic bacteriophage genome in FASTA format.
    Includes a biologically plausible structure:
      - Flanking non-coding sequences (5' UTR)
      - A valid open reading frame (ORF) beginning with 'ATG' and ending with 'TAA'
      - Flanking 3' UTR sequence
      - Standard 80-column line wrapping
    """
    fasta_file = folder_path / "small_phage.fasta"
    print(f"\nGenerating valid synthetic bacteriophage genome: {fasta_file.name}")
    
    # 50bp 5' UTR
    utr_5 = "AGTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
    
    # 393bp Open Reading Frame (131 codons)
    # Starts with ATG, ends with TAA
    orf_body = (
        "ATGGCCCGTGCACGTGCACGTGGTTTTACACGTGCACGTCGTGGTGCAGTTACACGTGCACGTG"
        "GTGGTGTTACACGTGCACGTCGTGGTGGTTTTACACGTGCACGTGGTGGTGCAGTTACACGTGC"
        "ACGTCGTGGTGGTGTTACACGTGCACGTCGTGGTGGTTTTACACGTGCACGTGGTGGTGCAGTT"
        "ACACGTGCACGTCGTGGTGGTGTTACACGTGCACGTCGTGGTGGTTTTACACGTGCACGTGGTG"
        "GTGCAGTTACACGTGCACGTCGTGGTGGTGTTACACGTGCACGTCGTGGTGGTTTTACACGTGC"
        "ACGTGGTGGTGCAGTTACACGTGCACGTCGTGGTGGTGTTACACGTGCACGTCGTGGTGGTTTT"
        "ACACGTTAA"
    )
    
    # 57bp 3' UTR
    utr_3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT"
    
    full_sequence = utr_5 + orf_body + utr_3 # Exactly 500bp
    
    assert len(full_sequence) == 500, f"Error: Generated sequence is {len(full_sequence)}bp instead of 500bp."
    
    try:
        with open(fasta_file, "w", encoding="utf-8") as f:
            f.write(">NC_000001.1 Synthetic Escherichia phage Lambda-mock, complete genome\n")
            # Wrap sequence lines at 80 characters
            for i in range(0, len(full_sequence), 80):
                f.write(full_sequence[i:i+80] + "\n")
        print(f"[SUCCESS] Wrote valid 500bp FASTA to {fasta_file}")
    except Exception as e:
        print(f"[ERROR] Failed to write phage FASTA: {e}", file=sys.stderr)

def generate_mock_fastq(folder_path: Path):
    """
    Generates a high-fidelity synthetic FASTQ file containing standard reads
    with Phred+33 quality scores. This includes:
      - 3 High-quality reads (average Q-score > 35)
      - 2 Low-quality reads (average Q-score < 15, to test sliding window and mean QC trims)
      - Standard Illumina header formats
    """
    fastq_file = folder_path / "mock_reads.fastq"
    print(f"\nGenerating high-fidelity synthetic raw reads: {fastq_file.name}")
    
    # 4-line FASTQ records
    reads = [
        # Read 1: High Quality, 50bp
        "@SRR123456.1 Illumina_mock_read_1 length=50\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n", # Phred 40 for all bases
        
        # Read 2: High Quality, 50bp
        "@SRR123456.2 Illumina_mock_read_2 length=50\n"
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\n"
        "+\n"
        "JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n", # Phred 41 for all bases
        
        # Read 3: High Quality, 50bp
        "@SRR123456.3 Illumina_mock_read_3 length=50\n"
        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        "+\n"
        "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n", # Phred 39 for all bases
        
        # Read 4: Low Quality on 3' end, 50bp (tests sliding window trimming)
        "@SRR123456.4 Illumina_mock_read_4_low_3prime length=50\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGNNNNNNNNNN\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!!!!!!!!!!\n", # Drops to Phred 0 on the last 10 bases
        
        # Read 5: Low Quality Overall, 50bp (tests mean Q-score filtering)
        "@SRR123456.5 Illumina_mock_read_5_poor_overall length=50\n"
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
        "+\n"
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" # Phred 0 for the whole read (discard candidate)
    ]
    
    try:
        with open(fastq_file, "w", encoding="utf-8") as f:
            for read in reads:
                f.write(read)
        print(f"[SUCCESS] Wrote synthetic raw reads FASTQ to {fastq_file}")
        print("[INFO] Real raw FASTQ sequencing files can also be downloaded and placed here.")
    except Exception as e:
        print(f"[ERROR] Failed to write mock FASTQ: {e}", file=sys.stderr)

def generate_invalid_characters_fasta(folder_path: Path):
    """
    Generates an invalid FASTA file containing non-IUPAC characters like 'X' and 'Z'
    to test the pipeline preprocessor's input validation constraints.
    """
    fasta_file = folder_path / "invalid_characters.fasta"
    print(f"\nGenerating invalid input dataset: {fasta_file.name}")
    
    # Sequence containing invalid non-IUPAC characters 'X' and 'Z'
    bad_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGXATCGATCGZATCG"
    
    try:
        with open(fasta_file, "w", encoding="utf-8") as f:
            f.write(">seq_invalid_chars FASTA containing invalid X and Z characters\n")
            for i in range(0, len(bad_sequence), 80):
                f.write(bad_sequence[i:i+80] + "\n")
        print(f"[SUCCESS] Wrote invalid FASTA to {fasta_file}")
    except Exception as e:
        print(f"[ERROR] Failed to write invalid FASTA: {e}", file=sys.stderr)

def generate_ambiguous_fasta(folder_path: Path):
    """
    Generates a FASTA sequence where exactly 35% of the bases are ambiguous ('N').
    This is used to validate the pipeline's 'max_ambiguous_pct' filtering logic
    (which typically drops sequences exceeding a 5% threshold).
    """
    fasta_file = folder_path / "ambiguous_35percent.fasta"
    print(f"\nGenerating ambiguous input dataset: {fasta_file.name}")
    
    # 35 'N' bases and 65 valid bases = exactly 35% ambiguity
    n_bases = "N" * 35
    valid_bases = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA" # 97bp
    # Trim valid bases to 65bp
    valid_bases = valid_bases[:65]
    
    full_seq = n_bases + valid_bases
    assert len(full_seq) == 100, f"Error: Ambiguous sequence length is {len(full_seq)}bp instead of 100bp."
    
    try:
        with open(fasta_file, "w", encoding="utf-8") as f:
            f.write(">seq_ambiguous_35pct Sequence containing exactly 35% N bases\n")
            for i in range(0, len(full_seq), 80):
                f.write(full_seq[i:i+80] + "\n")
        print(f"[SUCCESS] Wrote ambiguous FASTA to {fasta_file}")
    except Exception as e:
        print(f"[ERROR] Failed to write ambiguous FASTA: {e}", file=sys.stderr)

def generate_duplicate_sequences_fasta(folder_path: Path):
    """
    Generates a FASTA sequence containing 3 identical sequence records with
    distinct headers to test the pipeline preprocessor's sequence collapsing/deduplication logic.
    """
    fasta_file = folder_path / "duplicate_sequences.fasta"
    print(f"\nGenerating duplicate sequences dataset: {fasta_file.name}")
    
    identical_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" # 120bp
    
    try:
        with open(fasta_file, "w", encoding="utf-8") as f:
            f.write(">seq_duplicate_1 Identical sequence duplicate record 1\n")
            f.write(identical_sequence + "\n")
            f.write(">seq_duplicate_2 Identical sequence duplicate record 2\n")
            f.write(identical_sequence + "\n")
            f.write(">seq_duplicate_3 Identical sequence duplicate record 3\n")
            f.write(identical_sequence + "\n")
        print(f"[SUCCESS] Wrote duplicate FASTA to {fasta_file}")
    except Exception as e:
        print(f"[ERROR] Failed to write duplicate FASTA: {e}", file=sys.stderr)

def print_real_dataset_instructions(paths: dict):
    """
    Prints high-quality biological documentation and instructions detailing how the
    user can download a real raw FASTQ sequencing run and copy it to their test folders.
    """
    instructions = f"""
================================================================================
                    BIOLOGICAL SEQUENCING DATA DOWNLOAD GUIDE
================================================================================
To test the pipeline with real biological sequencing datasets, you can download 
a raw FASTQ file from the NCBI Sequence Read Archive (SRA) or European Nucleotide Archive (ENA).

Recommended Dataset: Escherichia Phage Lambda (SRR11246031)
- Platform: Illumina MiSeq paired-end sequencing
- Purpose: High-coverage, small-genome test case ideal for local development pipelines.

INSTRUCTIONS:
1. Open your terminal or a browser.
2. Download R1 (Forward) and R2 (Reverse) fastq files:
   Forward: https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR112/031/SRR11246031/SRR11246031.fastq.gz
   (Optional: You can download and run paired-end QC!)

3. Save the downloaded gzipped FASTQ file directly into the target folder:
   Path: {paths["06_fastq_raw_reads"].resolve()}

4. The pipeline coordinator and preprocessor will automatically detect the gzipped
   fastq format (.fastq.gz) and execute quality control and adapter trimming using
   either the native 'fastp' binary or the high-fidelity pure-Python fallback parser.
================================================================================
"""
    print(instructions)
    
    # Also write these instructions to a README file in the FASTQ folder
    readme_path = paths["06_fastq_raw_reads"] / "README_DOWNLOAD.txt"
    try:
        with open(readme_path, "w", encoding="utf-8") as f:
            f.write(instructions)
    except Exception:
        pass

def main():
    """Main orchestrator for generating the test folders and synthetic files."""
    # Find the current project root directory (current directory of setup_test_data.py)
    project_root = Path(os.path.dirname(os.path.abspath(__file__))) if "__file__" in locals() else Path(os.getcwd())
    
    test_data_root = project_root / "test_data"
    
    # 1. Create directory structure
    paths = create_directory_structure(test_data_root)
    
    # 2. Populate synthetic datasets
    generate_valid_phage_fasta(paths["01_small_bacteriophage"])
    generate_mock_fastq(paths["06_fastq_raw_reads"])
    generate_invalid_characters_fasta(paths["07_invalid_input"])
    generate_ambiguous_fasta(paths["08_ambiguous_sequences"])
    generate_duplicate_sequences_fasta(paths["09_duplicate_sequences"])
    
    # 3. Print biological guide for real datasets
    print_real_dataset_instructions(paths)
    
    print("\n=== PathoScope AI Test Data Setup Complete ===")
    print(f"All synthetic datasets are populated inside: {test_data_root.resolve()}")
    print("You can now run your test and validation suites using these datasets!\n")

if __name__ == "__main__":
    main()
