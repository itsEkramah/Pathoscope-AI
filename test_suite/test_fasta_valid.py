import unittest
from pathlib import Path
from pathoscope.core.preprocessor import parse_input_file, validate_and_clean_sequence
from pathoscope.utils.config_loader import load_config

class TestFastaValid(unittest.TestCase):
    """
    Bioinformatics QA Test Case: Verifies that valid FASTA assemblies parse cleanly 
    and qualify through standard quality control checks.
    """
    def test_fasta_ingestion(self):
        fasta_path = Path("test_data/01_small_bacteriophage/small_phage.fasta")
        self.assertTrue(fasta_path.exists(), f"Target test data file does not exist: {fasta_path}")
        
        # 1. Parse FASTA records
        records = parse_input_file(fasta_path)
        self.assertGreater(len(records), 0, "No records parsed from small phage assembly.")
        
        rec = records[0]
        self.assertIn("NC_", rec.id) # Standard RefSeq NCBI header id
        
        # 2. Ingest configurations and validate sequence details
        config = load_config(Path("config/default_config.yaml"))
        is_valid, reason, stats = validate_and_clean_sequence(
            rec,
            min_len=config.preprocessing.min_length,
            max_len=config.preprocessing.max_length,
            max_ambig_pct=config.preprocessing.max_ambiguous_pct
        )
        self.assertTrue(is_valid, f"FASTA Sequence Validation failed unexpectedly: {reason}")
        self.assertGreater(stats["length"], 100)
        self.assertGreater(stats["gc_percent"], 0.0)
        self.assertEqual(stats["ambiguous_count"], 0)
