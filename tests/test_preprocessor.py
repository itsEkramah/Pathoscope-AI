# pyrefly: ignore [missing-import]
import pytest
from pathlib import Path
import tempfile
import pandas as pd
from Bio import SeqIO
from pathoscope.core.preprocessor import (
    SequenceRecord,
    parse_input_file,
    validate_and_clean_sequence,
    calculate_n50,
    process_sequences,
    SequenceValidationError,
    InvalidFileFormatError
)
from pathoscope.utils.config_loader import AppConfig

# 1. Test basic parsing and representation
def test_sequence_record_init():
    rec = SequenceRecord("seq1", "atcgtgac")
    assert rec.id == "seq1"
    assert rec.sequence == "ATCGTGAC" # Normalized to uppercase
    assert rec.length == 8

def test_sequence_record_invalid_init():
    with pytest.raises(SequenceValidationError):
        SequenceRecord("", "ATCG")
    with pytest.raises(SequenceValidationError):
        SequenceRecord("seq1", "")

# 2. Test quality control filters
def test_validate_and_clean_valid():
    record = SequenceRecord("seq_ok", "ATCGATCGATCG")
    is_valid, reason, stats = validate_and_clean_sequence(
        record, min_len=10, max_len=100, max_ambig_pct=5.0
    )
    assert is_valid is True
    assert reason == ""
    assert stats["length"] == 12
    assert stats["gc_percent"] == 50.0
    assert stats["ambiguous_count"] == 0

def test_validate_and_clean_length_too_short():
    record = SequenceRecord("seq_short", "ATCG")
    is_valid, reason, stats = validate_and_clean_sequence(
        record, min_len=10, max_len=100, max_ambig_pct=5.0
    )
    assert is_valid is False
    assert "below minimum threshold" in reason

def test_validate_and_clean_length_too_long():
    record = SequenceRecord("seq_long", "A" * 150)
    is_valid, reason, stats = validate_and_clean_sequence(
        record, min_len=10, max_len=100, max_ambig_pct=5.0
    )
    assert is_valid is False
    assert "exceeds maximum threshold" in reason

def test_validate_and_clean_invalid_chars():
    record = SequenceRecord("seq_bad_char", "ATCGZATCG") # 'Z' is not IUPAC
    is_valid, reason, stats = validate_and_clean_sequence(
        record, min_len=5, max_len=100, max_ambig_pct=5.0
    )
    assert is_valid is False
    assert "non-IUPAC invalid characters" in reason

def test_validate_and_clean_ambiguous_bases():
    record = SequenceRecord("seq_ambig", "ATCGNNATCG") # 2 'N's out of 10 bp = 20%
    is_valid, reason, stats = validate_and_clean_sequence(
        record, min_len=5, max_len=100, max_ambig_pct=5.0
    )
    assert is_valid is False
    assert "Ambiguous bases" in reason
    
    # 'N's under threshold
    is_valid, reason, stats = validate_and_clean_sequence(
        record, min_len=5, max_len=100, max_ambig_pct=25.0
    )
    assert is_valid is True
    assert stats["ambiguous_percent"] == 20.0

# 3. Test biological N50 assembly calculation
def test_calculate_n50():
    lengths = [100, 200, 300, 400, 500] # total = 1500, half = 750
    # Sorted: 500, 400, 300, 200, 100
    # Cum sum: 500 (not 750), 900 (>= 750). N50 is 400.
    assert calculate_n50(lengths) == 400
    
    assert calculate_n50([]) == 0
    assert calculate_n50([100]) == 100

# 4. Test parsing of multiple file formats
def test_parse_fasta_format():
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
        f.write(">seq1 test sequence\nATCG\n>seq2 another sequence\nGCTA\n")
        temp_path = Path(f.name)
        
    try:
        recs = parse_input_file(temp_path)
        assert len(recs) == 2
        assert recs[0].id == "seq1"
        assert recs[0].sequence == "ATCG"
        assert recs[1].id == "seq2"
    finally:
        temp_path.unlink()

# 5. Test complete sequence processing and de-duplication workflow
def test_process_sequences_duplicate_headers_reject():
    # Setup mock sequences with a duplicate header ID
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
        f.write(
            ">good_seq\nATCGATCGATCGATCGATCG\n" # 20bp
            ">good_seq\nCGATCGATCGATCGATCGAT\n" # duplicate header, different sequence
        )
        temp_input = Path(f.name)
        
    temp_output = Path(tempfile.mktemp(suffix=".fasta"))
    
    # Configure mock configuration
    config = AppConfig()
    config.preprocessing.min_length = 10
    config.preprocessing.handle_duplicate_headers = "reject"
    
    try:
        qc_summary = process_sequences(temp_input, temp_output, config)
        
        # Verify duplicate header was rejected
        assert qc_summary["counts"]["total_processed"] == 2
        assert qc_summary["counts"]["total_kept"] == 1
        assert qc_summary["counts"]["total_discarded"] == 1
        assert qc_summary["counts"]["duplicate_headers_encountered"] == 1
        
        # Check generated files
        parent = temp_output.parent
        assert (parent / "rejected.fasta").exists()
        assert (parent / "qc_statistics.csv").exists()
        assert (parent / "preprocessing_report.json").exists()
        
        # Check rejected file headers contain reasons
        recs_rej = list(SeqIO.parse(str(parent / "rejected.fasta"), "fasta"))
        assert len(recs_rej) == 1
        assert "REJECTED: Duplicate header" in recs_rej[0].description
        
        # Check CSV contents
        df_csv = pd.read_csv(parent / "qc_statistics.csv")
        assert len(df_csv) == 2
        assert df_csv.iloc[1]["status"] == "REJECTED"
        assert "Duplicate header" in df_csv.iloc[1]["rejection_reason"]
        
    finally:
        if temp_input.exists():
            temp_input.unlink()
        if temp_output.exists():
            temp_output.unlink()
            (temp_output.parent / "rejected.fasta").unlink()
            (temp_output.parent / "qc_statistics.csv").unlink()
            (temp_output.parent / "preprocessing_report.json").unlink()


def test_process_sequences_duplicate_headers_rename():
    # Setup mock sequences with a duplicate header ID
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
        f.write(
            ">good_seq\nATCGATCGATCGATCGATCG\n" # 20bp
            ">good_seq\nCGATCGATCGATCGATCGAT\n" # duplicate header, different sequence
        )
        temp_input = Path(f.name)
        
    temp_output = Path(tempfile.mktemp(suffix=".fasta"))
    
    # Configure mock configuration to rename
    config = AppConfig()
    config.preprocessing.min_length = 10
    config.preprocessing.handle_duplicate_headers = "rename"
    
    try:
        qc_summary = process_sequences(temp_input, temp_output, config)
        
        # Verify duplicate header was renamed and KEPT
        assert qc_summary["counts"]["total_processed"] == 2
        assert qc_summary["counts"]["total_kept"] == 2
        assert qc_summary["counts"]["total_discarded"] == 0
        assert qc_summary["counts"]["duplicate_headers_encountered"] == 1
        
        # Assert renamed ID
        cleaned_recs = list(SeqIO.parse(str(temp_output), "fasta"))
        assert len(cleaned_recs) == 2
        assert cleaned_recs[0].id == "good_seq"
        assert cleaned_recs[1].id == "good_seq_dup1"
        
    finally:
        if temp_input.exists():
            temp_input.unlink()
        if temp_output.exists():
            temp_output.unlink()
            (temp_output.parent / "rejected.fasta").unlink()
            (temp_output.parent / "qc_statistics.csv").unlink()
            (temp_output.parent / "preprocessing_report.json").unlink()


def test_process_sequences_duplicate_contents_keep_flag():
    # Setup mock sequences with duplicate contents but different headers
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
        f.write(
            ">seq_first\nATCGATCGATCGATCGATCG\n" # 20bp
            ">seq_second_dup\nATCGATCGATCGATCGATCG\n" # duplicate sequence contents
        )
        temp_input = Path(f.name)
        
    temp_output = Path(tempfile.mktemp(suffix=".fasta"))
    
    # Configure to keep duplicates (remove_duplicate_sequences = False)
    config = AppConfig()
    config.preprocessing.min_length = 10
    config.preprocessing.remove_duplicate_sequences = False
    
    try:
        qc_summary = process_sequences(temp_input, temp_output, config)
        
        # Verify both kept
        assert qc_summary["counts"]["total_processed"] == 2
        assert qc_summary["counts"]["total_kept"] == 2
        assert qc_summary["counts"]["total_discarded"] == 0
        assert qc_summary["counts"]["duplicate_sequences_encountered"] == 1
        
        # Verify flag is appended in description
        recs = list(SeqIO.parse(str(temp_output), "fasta"))
        assert len(recs) == 2
        assert "Duplicate Content Match: seq_first" in recs[1].description
        
    finally:
        if temp_input.exists():
            temp_input.unlink()
        if temp_output.exists():
            temp_output.unlink()
            (temp_output.parent / "rejected.fasta").unlink()
            (temp_output.parent / "qc_statistics.csv").unlink()
            (temp_output.parent / "preprocessing_report.json").unlink()


def test_process_sequences_duplicate_contents_remove():
    # Setup mock sequences with duplicate contents
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
        f.write(
            ">seq_first\nATCGATCGATCGATCGATCG\n" # 20bp
            ">seq_second_dup\nATCGATCGATCGATCGATCG\n" # duplicate sequence contents
        )
        temp_input = Path(f.name)
        
    temp_output = Path(tempfile.mktemp(suffix=".fasta"))
    
    # Configure to remove duplicates (remove_duplicate_sequences = True)
    config = AppConfig()
    config.preprocessing.min_length = 10
    config.preprocessing.remove_duplicate_sequences = True
    
    try:
        qc_summary = process_sequences(temp_input, temp_output, config)
        
        # Verify second sequence is removed
        assert qc_summary["counts"]["total_processed"] == 2
        assert qc_summary["counts"]["total_kept"] == 1
        assert qc_summary["counts"]["total_discarded"] == 1
        assert qc_summary["counts"]["duplicate_sequences_encountered"] == 1
        
        # Check rejected file contains rejection explanation
        parent = temp_output.parent
        recs_rej = list(SeqIO.parse(str(parent / "rejected.fasta"), "fasta"))
        assert len(recs_rej) == 1
        assert "REJECTED: Duplicate sequence content" in recs_rej[0].description
        
    finally:
        if temp_input.exists():
            temp_input.unlink()
        if temp_output.exists():
            temp_output.unlink()
            (temp_output.parent / "rejected.fasta").unlink()
            (temp_output.parent / "qc_statistics.csv").unlink()
            (temp_output.parent / "preprocessing_report.json").unlink()


def test_validate_paired_end_synchronicity_mismatched_counts():
    from pathoscope.core.preprocessor import validate_paired_end_synchronicity
    
    with tempfile.NamedTemporaryFile(suffix="_R1.fastq", delete=False, mode="w") as f1, \
         tempfile.NamedTemporaryFile(suffix="_R2.fastq", delete=False, mode="w") as f2:
        f1.write("@read1\nATCG\n+\n!!!!\n")
        f2.write("@read1\nATCG\n+\n!!!!\n@read2\nCGTA\n+\n!!!!\n")
        r1_path = Path(f1.name)
        r2_path = Path(f2.name)
        
    try:
        with pytest.raises(SequenceValidationError):
            validate_paired_end_synchronicity(r1_path, r2_path)
    finally:
        r1_path.unlink()
        r2_path.unlink()


def test_validate_paired_end_synchronicity_mismatched_headers():
    from pathoscope.core.preprocessor import validate_paired_end_synchronicity
    
    with tempfile.NamedTemporaryFile(suffix="_R1.fastq", delete=False, mode="w") as f1, \
         tempfile.NamedTemporaryFile(suffix="_R2.fastq", delete=False, mode="w") as f2:
        f1.write("@read1\nATCG\n+\n!!!!\n")
        f2.write("@read2\nATCG\n+\n!!!!\n")
        r1_path = Path(f1.name)
        r2_path = Path(f2.name)
        
    try:
        with pytest.raises(SequenceValidationError):
            validate_paired_end_synchronicity(r1_path, r2_path)
    finally:
        r1_path.unlink()
        r2_path.unlink()


def test_calculate_high_fidelity_qc_distributions():
    from pathoscope.core.preprocessor import calculate_high_fidelity_qc_distributions
    
    with tempfile.NamedTemporaryFile(suffix="_raw.fastq", delete=False, mode="w") as f1, \
         tempfile.NamedTemporaryFile(suffix="_clean.fastq", delete=False, mode="w") as f2:
        f1.write("@read1\nATCGATCGATCG\n+\n!!!!!!!!!!!!\n")
        f2.write("@read1\nATCGATCGATCG\n+\n!!!!!!!!!!!!\n")
        raw_path = Path(f1.name)
        clean_path = Path(f2.name)
        
    try:
        dist = calculate_high_fidelity_qc_distributions(
            raw_r1=raw_path,
            raw_r2=None,
            clean_r1=clean_path,
            clean_r2=None,
            adapter_trimmed_count=1,
            rejection_reasons={"Low average quality...": 0}
        )
        assert len(dist["per_base_quality"]) == 12
        assert dist["per_base_quality"][0] == 0.0
        assert dist["sequence_duplication_levels"]["1"] == 1
        assert dist["read_retention_waterfall"]["raw_reads"] == 1
    finally:
        raw_path.unlink()
        clean_path.unlink()


def test_max_reads_cap_and_vectorization():
    from pathoscope.core.preprocessor import process_sequences, sliding_window_quality_trim, FastqRecord
    
    # Test vectorized sliding window quality trimming on a high-quality record
    rec_high = FastqRecord("@read_high", "ATCGATCGATCGATCG", "+", "IIIIIIIIIIIIIIII")
    trimmed, was_trimmed = sliding_window_quality_trim(rec_high, window_size=4, min_qual=20)
    assert not was_trimmed
    assert trimmed.length == 16
    
    # Test on a record with low-quality tail
    rec_low = FastqRecord("@read_low", "ATCGATCGATCGATCG", "+", "IIIIIIIIIIII!!!!")
    trimmed, was_trimmed = sliding_window_quality_trim(rec_low, window_size=4, min_qual=20)
    assert was_trimmed
    assert trimmed.length == 11  # Trimmed the low-quality tail
    
    # Test max_reads_cap sub-sampling
    with tempfile.NamedTemporaryFile(suffix="_raw.fastq", delete=False, mode="w") as f:
        # Write 15 reads (60 lines)
        for i in range(15):
            f.write(f"@read_{i}\nATCGATCGATCGATCG\n+\nIIIIIIIIIIIIIIII\n")
        temp_input = Path(f.name)
        
    temp_output = Path(tempfile.mktemp(suffix=".fasta"))
    
    config = AppConfig()
    config.preprocessing.min_length = 10
    config.preprocessing.max_reads_cap = 10  # Set cap lower than available reads
    config.preprocessing.paired_end = False
    
    try:
        summary = process_sequences(temp_input, temp_output, config)
        # Should process exactly max_reads_cap (10) reads
        assert summary["counts"]["total_processed"] == 10
        assert temp_output.exists()
    finally:
        if temp_input.exists():
            temp_input.unlink(missing_ok=True)
        if temp_output.exists():
            temp_output.unlink(missing_ok=True)
            (temp_output.parent / "rejected.fasta").unlink(missing_ok=True)
            (temp_output.parent / "qc_statistics.csv").unlink(missing_ok=True)
            (temp_output.parent / "preprocessing_report.json").unlink(missing_ok=True)
            (temp_output.parent / "qc_report.html").unlink(missing_ok=True)
            (temp_output.parent / "cleaned_R1.fastq.gz").unlink(missing_ok=True)
            (temp_output.parent / "audit_stats.csv").unlink(missing_ok=True)
