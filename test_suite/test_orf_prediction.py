import unittest
from pathlib import Path
import tempfile
from pathoscope.core.preprocessor import parse_input_file
from pathoscope.core.orf_predictor import predict_orfs_in_sequence, process_orf_prediction
from pathoscope.utils.config_loader import load_config

class TestOrfPrediction(unittest.TestCase):
    """
    Bioinformatics QA Test Case: Evaluates the six-frame coordinate-aware ORF Predictor,
    asserting reverse strand conversions, ribosomal outermost sweeps, and multi-segment processing.
    """
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.temp_dir.name)
        self.config = load_config(Path("config/default_config.yaml"))

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_bacteriophage_coordinate_remapping(self):
        # Folder 01: bacteriophage coordinates
        fasta_path = Path("test_data/01_small_bacteriophage/small_phage.fasta")
        self.assertTrue(fasta_path.exists(), f"Target test data file does not exist: {fasta_path}")
        
        records = parse_input_file(fasta_path)
        self.assertGreater(len(records), 0)
        
        # Enable TT11 for prokaryotes/phages
        self.config.orf_prediction.translation_table = 11
        self.config.orf_prediction.min_orf_length_aa = 30
        
        orfs = predict_orfs_in_sequence(records[0], self.config)
        self.assertGreater(len(orfs), 0, "Failed to predict any ORFs on bacteriophage lambda-mock.")
        
        # Verify ORF records contain proper 1-based coordinates
        for orf in orfs:
            self.assertGreaterEqual(orf.start, 1)
            self.assertLessEqual(orf.end, records[0].length)
            self.assertIn(orf.strand, ["+", "-"])
            self.assertGreater(orf.length_aa, 0)
            self.assertGreater(orf.confidence_score, 0.0)

    def test_segmented_influenza_multi_sequence(self):
        # Folder 04: Segmented Influenza (multi-sequence FASTA)
        fasta_path = Path("test_data/04_segmented_virus/test_04_influenza_A.fasta")
        if not fasta_path.exists():
            fasta_path = Path("test_data/04_segmented_virus/test_04_influenza_A.fasta.fasta")
            
        self.assertTrue(fasta_path.exists(), f"Target test data file does not exist: {fasta_path}")
        
        records = parse_input_file(fasta_path)
        # Should contain multiple segments (Influenza A contains 8 segments)
        self.assertGreaterEqual(len(records), 8, f"Influenza A segmented fasta should contain at least 8 segments. Found: {len(records)}")
        
        # Test full ORF prediction directory process
        orf_outdir = self.output_dir / "orfs"
        orf_stats = process_orf_prediction(fasta_path, orf_outdir, self.config)
        
        self.assertGreater(orf_stats["counts"]["total_orfs_predicted"], 0)
        self.assertTrue((orf_outdir / "coordinates.gff3").exists())
        self.assertTrue((orf_outdir / "proteins.fasta").exists())

    def test_nested_orf_resolution_sars_cov2(self):
        # Folder 02: Positive-sense RNA complete genome (SARS-CoV-2)
        fasta_path = Path("test_data/02_ssRNA_virus/NC_045512.2.fasta")
        self.assertTrue(fasta_path.exists(), f"Target test data file does not exist: {fasta_path}")
        
        records = parse_input_file(fasta_path)
        self.assertEqual(len(records), 1)
        
        # Configure nested ORF resolution
        self.config.orf_prediction.resolve_nested = True
        self.config.orf_prediction.min_orf_length_aa = 30
        
        orfs_resolved = predict_orfs_in_sequence(records[0], self.config)
        
        # If we configure nested resolution, overlapping shorter frames are removed, keeping the longest
        self.assertGreater(len(orfs_resolved), 0)
        for orf in orfs_resolved:
            self.assertNotEqual(orf.overlap_flag, "Nested")  # All resolved nested ORFs are deleted
            
    def test_large_poxvirus_stress_test(self):
        # Folder 05: Large Poxvirus Genome (stress test)
        fasta_path = Path("test_data/05_large_virus/BK066888.1 .fasta.fasta")
        self.assertTrue(fasta_path.exists(), f"Target test data file does not exist: {fasta_path}")
        
        records = parse_input_file(fasta_path)
        self.assertEqual(len(records), 1)
        
        orf_outdir = self.output_dir / "orfs"
        self.config.orf_prediction.min_orf_length_aa = 50 # slightly larger minimum to avoid spurious hits
        
        orf_stats = process_orf_prediction(fasta_path, orf_outdir, self.config)
        # Should predict many ORFs since poxviruses are large double-stranded DNA genomes (100kbp+)
        self.assertGreaterEqual(orf_stats["counts"]["total_orfs_predicted"], 30)
