# pyrefly: ignore [missing-import]
import pytest
from pathlib import Path
import tempfile
import pandas as pd
from Bio import SeqIO
from pathoscope.core.preprocessor import SequenceRecord
from pathoscope.core.orf_predictor import (
    ORFRecord,
    reverse_complement_coordinates,
    scan_single_strand,
    resolve_overlaps,
    predict_orfs_in_sequence,
    process_orf_prediction,
    InvalidCoordinateError
)
from pathoscope.utils.config_loader import AppConfig

# 1. Test coordinate mappings
def test_reverse_complement_coordinates():
    # Sequence length L = 100 bp
    # RC ORF starts at 0-based index 10 of reverse complement, and ends at 0-based index 40 (inclusive).
    # Expected: forward coordinates.
    # index 10 on RC corresponds to L - 1 - 10 = 89 on forward fwd.
    # index 40 on RC corresponds to L - 1 - 40 = 59 on forward fwd.
    # 1-based coordinates: start_fwd = 59 + 1 = 60, end_fwd = 89 + 1 = 90.
    start_fwd, end_fwd = reverse_complement_coordinates(10, 40, 100)
    assert start_fwd == 60
    assert end_fwd == 90
    
    with pytest.raises(InvalidCoordinateError):
        reverse_complement_coordinates(40, 10, 100)

# 2. Test basic ORFRecord properties
def test_orf_record_properties():
    orf = ORFRecord(
        orf_id="ORF1",
        sequence_id="viral_seq",
        start=10,
        end=70,
        strand="+",
        frame=1,
        nucleotide_seq="ATG" + "A" * 54 + "TAA", # 60bp
        protein_seq="M" + "K" * 18,               # 19aa (excluding stop)
        start_codon="ATG"
    )
    assert orf.length_bp == 60
    assert orf.length_aa == 19
    
    d = orf.to_dict()
    assert d["orf_id"] == "ORF1"
    assert d["length_aa"] == 19
    assert d["strand"] == "+"
    
    gff = orf.to_gff_line()
    assert "viral_seq\tPathoScope_ORF\tCDS\t10\t70\t0.0\t+\t0\tID=ORF1;Name=ORF1;frame=1;start_codon=ATG;overlap=None;confidence_score=0.0" in gff

# 3. Test true scanning on forward strand
def test_scan_single_strand_forward():
    # Generate mock sequence: starts withATG, has 30 bases, ends with TAA (40bp total)
    # Plus frame +1: ATG (index 0) ... TAA (index 33)
    seq = "ATGGCAAAATTTGGGCCCTTTGGGCCCTTTGGGTAA" # length 36 bp -> 12 aa (11 aa without stop)
    
    orfs = scan_single_strand(
        sequence_id="test_seq",
        seq_str=seq,
        strand="+",
        start_codons=["ATG"],
        stop_codons=["TAA"],
        min_len_aa=5,
        translation_table=1
    )
    
    assert len(orfs) == 1
    assert orfs[0].start == 1
    assert orfs[0].end == 36
    assert orfs[0].strand == "+"
    assert orfs[0].frame == 1
    assert orfs[0].start_codon == "ATG"
    assert orfs[0].prot_seq.startswith("M")

# 4. Test prevention of nested ORFs
def test_scan_single_strand_nested_prevention():
    # Sequence with two ATGs in same frame ending in TAA
    # Frame 1: ATG (idx 0) ... ATG (idx 9) ... TAA (idx 18)
    seq = "ATGGCAAAAATGGCAAAATAA" # length 21 bp
    # If nested genes resolved, we should ONLY get the outmost (longest) starting at index 0.
    # Longest: ATGGCAAAAATGGCAAAATAA -> length 21 (7 aa including stop, 6 aa without stop)
    
    orfs = scan_single_strand(
        sequence_id="test_seq",
        seq_str=seq,
        strand="+",
        start_codons=["ATG"],
        stop_codons=["TAA"],
        min_len_aa=2,
        translation_table=1
    )
    
    assert len(orfs) == 1
    assert orfs[0].start == 1
    assert orfs[0].end == 21
    assert orfs[0].length_bp == 21
    assert orfs[0].prot_seq == "MAKMAK"

# 5. Test overlap resolution algorithms
def test_resolve_overlaps_longest_only():
    # Two predicted ORFs that overlap:
    # ORF A: start 10, end 100 (90bp)
    # ORF B: start 50, end 120 (70bp) - Overlap is 50 bp
    orf_a = ORFRecord("A", "seq", 10, 100, "+", 1, "N" * 90, "M" * 29, "ATG")
    orf_b = ORFRecord("B", "seq", 50, 120, "+", 2, "N" * 70, "M" * 22, "ATG")
    
    # Under longest_only, ORF B should be discarded
    resolved = resolve_overlaps([orf_a, orf_b], overlap_threshold=30, policy="longest_only")
    assert len(resolved) == 1
    assert resolved[0].id == "A"

def test_resolve_overlaps_keep_all_flag():
    # Under keep_all_flag, both overlapping ORFs should be kept but flagged
    orf_a = ORFRecord("A", "seq", 10, 100, "+", 1, "N" * 90, "M" * 29, "ATG")
    orf_b = ORFRecord("B", "seq", 50, 120, "+", 2, "N" * 70, "M" * 22, "ATG")
    
    resolved = resolve_overlaps([orf_a, orf_b], overlap_threshold=30, policy="keep_all_flag")
    assert len(resolved) == 2
    assert resolved[0].overlap_flag == "Overlap"
    assert resolved[1].overlap_flag == "Overlap"
    assert "B" in resolved[0].overlap_with
    assert "A" in resolved[1].overlap_with

# 6. Test full process_orf_prediction orchestrator
def test_process_orf_prediction_workflow():
    # Setup standard mock sequence in clean FASTA
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
        f.write(">viral_isolate\nATGGCAAAATTTGGGCCCTTTGGGCCCTTTGGGTAA\n") # 36bp -> frame 1 ORF
        temp_cleaned = Path(f.name)
        
    temp_outdir = Path(tempfile.mkdtemp())
    
    config = AppConfig()
    config.orf_prediction.min_orf_length_aa = 5
    config.orf_prediction.overlap_resolution_policy = "keep_all_flag"
    
    try:
        stats = process_orf_prediction(temp_cleaned, temp_outdir, config)
        
        # Verify predicted stats
        assert stats["counts"]["total_orfs_predicted"] == 1
        
        # Assert files created
        assert (temp_outdir / "coordinates.gff3").exists()
        assert (temp_outdir / "coordinates.csv").exists()
        assert (temp_outdir / "proteins.fasta").exists()
        assert (temp_outdir / "orf_statistics.json").exists()
        
        # Read translated protein file
        recs = list(SeqIO.parse(str(temp_outdir / "proteins.fasta"), "fasta"))
        assert len(recs) == 1
        assert recs[0].id == "viral_isolate_predicted_orf_1"
        assert recs[0].seq.startswith("M")
        
    finally:
        if temp_cleaned.exists():
            temp_cleaned.unlink()
        for p in temp_outdir.iterdir():
            p.unlink()
        temp_outdir.rmdir()


# 7. Test in silico confidence score calculation
def test_orf_confidence_score():
    from pathoscope.core.orf_predictor import calculate_confidence
    # Long ORF with standard start codon ATG
    score_atg = calculate_confidence(length_aa=120, start_codon="ATG", nuc_seq="ATG" + "A"*357 + "TAA", bg_gc_content=0.45)
    assert 0.0 <= score_atg <= 1.0
    
    # Same ORF but alternative start codon GTG
    score_gtg = calculate_confidence(length_aa=120, start_codon="GTG", nuc_seq="GTG" + "A"*357 + "TAA", bg_gc_content=0.45)
    assert score_gtg < score_atg
    
    # Very short ORF gets penalized by sigmoidal length scaling
    score_short = calculate_confidence(length_aa=15, start_codon="ATG", nuc_seq="ATG" + "A"*42 + "TAA", bg_gc_content=0.45)
    assert score_short < score_atg


# 8. Test nested ORF filtering (both resolve and flag policies)
def test_nested_orf_filtering():
    from pathoscope.core.orf_predictor import filter_nested_orfs
    # Setup standard non-overlapping ORFs and nested child ORF
    parent = ORFRecord("parent", "seq", 10, 100, "+", 1, "N"*90, "M"*29, "ATG")
    child = ORFRecord("child", "seq", 20, 80, "+", 2, "N"*60, "M"*19, "ATG")
    
    # With resolve_nested=True, child must be filtered out
    resolved = filter_nested_orfs([parent, child], resolve_nested=True)
    assert len(resolved) == 1
    assert resolved[0].id == "parent"
    
    # With resolve_nested=False, child must be kept but flagged as Nested
    flagged = filter_nested_orfs([parent, child], resolve_nested=False)
    assert len(flagged) == 2
    assert flagged[0].overlap_flag == "None"
    assert flagged[1].overlap_flag == "Nested"
    assert "parent" in flagged[1].overlap_with


# 9. Test comparative literature genome alignment
def test_literature_comparison():
    from pathoscope.core.orf_predictor import compare_predicted_to_literature
    
    # Mock predicted ORFs matching MS2 maturation protein (ref: length 393 aa, start 130)
    maturation_orf = ORFRecord("predicted_mat", "seq", 135, 1315, "+", 1, "N"*1182, "M"*393, "ATG")
    # Mock predicted ORF matching MS2 coat protein (ref: length 130 aa, start 1335)
    coat_orf = ORFRecord("predicted_coat", "seq", 1336, 1726, "+", 1, "N"*390, "M"*130, "ATG")
    
    # Sequence length in range of MS2 (approx 3569 bp)
    comp_res = compare_predicted_to_literature([maturation_orf, coat_orf], seq_len=3569, seq_id="Bacteriophage MS2")
    assert comp_res["matched_reference_genome"] == "Bacteriophage MS2 (NC_001417)"
    
    matches = {m["reference_protein_name"]: m for m in comp_res["matches"]}
    
    # Coat protein should be aligned
    assert "Coat protein" in matches
    assert matches["Coat protein"]["predicted_orf_id"] == "predicted_coat"
    assert matches["Coat protein"]["match_confidence"] > 0.8
    
    # Maturation protein should be aligned
    assert "Maturation protein (A)" in matches
    assert matches["Maturation protein (A)"]["predicted_orf_id"] == "predicted_mat"
    assert matches["Maturation protein (A)"]["match_confidence"] > 0.8
