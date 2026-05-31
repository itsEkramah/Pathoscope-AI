import unittest
from pathlib import Path
import tempfile
from pathoscope.core.preprocessor import (
    parse_input_file,
    validate_and_clean_sequence,
    process_sequences
)
from pathoscope.utils.config_loader import load_config

class TestInvalidSequences(unittest.TestCase):
    """
    Bioinformatics QA Test Case: Challenges the validation parser against extreme 
    anomalies (non-IUPAC character insertions, excessive N-bases, and duplicate headers/sequence records).
    """
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_non_iupac_characters(self):
        # Folder 07: Invalid non-IUPAC characters
        bad_fasta = Path("test_data/07_invalid_input/invalid_characters.fasta")
        self.assertTrue(bad_fasta.exists(), f"Target test data file does not exist: {bad_fasta}")
        
        records = parse_input_file(bad_fasta)
        self.assertEqual(len(records), 1)
        
        is_valid, reason, stats = validate_and_clean_sequence(
            records[0], min_len=10, max_len=500000, max_ambig_pct=5.0
        )
        self.assertFalse(is_valid, "Sequence containing X and Z was not flagged as invalid.")
        self.assertIn("non-IUPAC invalid characters", reason)

    def test_ambiguity_threshold(self):
        # Folder 08: Ambiguous sequences (>5% N bases)
        ambig_fasta = Path("test_data/08_ambiguous_sequences/ambiguous_35percent.fasta")
        self.assertTrue(ambig_fasta.exists(), f"Target test data file does not exist: {ambig_fasta}")
        
        output_fasta = self.output_dir / "cleaned_ambig.fasta"
        config = load_config(Path("config/default_config.yaml"))
        config.preprocessing.min_length = 10
        config.preprocessing.max_ambiguous_pct = 5.0
        
        summary = process_sequences(ambig_fasta, output_fasta, config)
        
        # 35% ambiguity exceeds 5% limit, so it must be discarded
        self.assertEqual(summary["counts"]["total_processed"], 1)
        self.assertEqual(summary["counts"]["total_kept"], 0)
        self.assertEqual(summary["counts"]["total_discarded"], 1)

    def test_sequence_de_duplication(self):
        # Folder 09: Duplicate sequences
        dup_fasta = Path("test_data/09_duplicate_sequences/duplicate_sequences.fasta")
        self.assertTrue(dup_fasta.exists(), f"Target test data file does not exist: {dup_fasta}")
        
        output_fasta = self.output_dir / "cleaned_unique.fasta"
        config = load_config(Path("config/default_config.yaml"))
        config.preprocessing.min_length = 10
        config.preprocessing.remove_duplicate_sequences = True
        config.preprocessing.handle_duplicate_headers = "rename" # prevent header collisions from discarding
        
        summary = process_sequences(dup_fasta, output_fasta, config)
        
        # 3 processed, 1 kept (unique), 2 discarded (sequence duplicates)
        self.assertEqual(summary["counts"]["total_processed"], 3)
        self.assertEqual(summary["counts"]["total_kept"], 1)
        self.assertEqual(summary["counts"]["total_discarded"], 2)
        self.assertEqual(summary["counts"]["duplicate_sequences_encountered"], 2)
