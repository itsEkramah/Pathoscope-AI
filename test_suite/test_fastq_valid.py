import unittest
from pathlib import Path
import tempfile
from pathoscope.core.preprocessor import process_sequences
from pathoscope.utils.config_loader import load_config

class TestFastqValid(unittest.TestCase):
    """
    Bioinformatics QA Test Case: Asserts that raw high-throughput sequencing datasets
    in FASTQ format are correctly audited, sliding window quality cutting trims low-quality 3' ends, 
    and reads failing overall quality thresholds are discarded.
    """
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_fastq_quality_filtering(self):
        real_fastq = Path("test_data/06_fastq_raw_reads/SRR13182871_2(2).fastq")
        if not real_fastq.exists():
            real_fastq = Path("test_data/06_fastq_raw_reads/SRR13182871_2.fastq.gz")
            
        self.assertTrue(real_fastq.exists(), f"Target real FASTQ test data file does not exist: {real_fastq}")
        
        # Squeeze out the first 100 reads (400 lines) of the real FASTQ file into a temporary FASTQ file for sub-second, real-biological testing
        sliced_fastq = self.output_dir / "real_subset_reads.fastq"
        
        import gzip
        open_func = gzip.open if str(real_fastq).endswith(".gz") else open
        mode = "rt"
        
        with open_func(real_fastq, mode, encoding="utf-8") as infile, open(sliced_fastq, "w", encoding="utf-8") as outfile:
            for _ in range(400):  # 100 reads * 4 lines = 400 lines
                line = infile.readline()
                if not line:
                    break
                outfile.write(line)
                
        output_fasta = self.output_dir / "cleaned_reads.fasta"
        
        # Load configuration and customize thresholds for real test reads
        config = load_config(Path("config/default_config.yaml"))
        config.preprocessing.min_length = 30
        config.preprocessing.min_mean_qscore = 30
        config.preprocessing.paired_end = False
        config.preprocessing.fastq_r2_path = None
        
        # Run preprocessor
        summary = process_sequences(sliced_fastq, output_fasta, config)
        
        # Assertions
        # Expecting to successfully process 100 reads sliced from the real FASTQ file
        self.assertEqual(summary["counts"]["total_processed"], 100)
        self.assertTrue(output_fasta.exists())
        # Verify that we discarded some poor-quality reads (as there are several full 'N' / low Q-score reads in the first 100)
        self.assertGreater(summary["counts"]["total_discarded"], 0, "No poor-quality reads discarded from real FASTQ slice; quality filters did not work.")
        self.assertGreater(summary["counts"]["total_kept"], 0, "All reads discarded from real FASTQ slice; too strict quality filters.")
