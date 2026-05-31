"""
test_expression.py

Unit tests for the DGE and Count Normalization module.
Validates library size CPM normalization, Welch's t-test evaluations, statsmodels FDR,
and upregulation/downregulation DEG class segmentations.
"""

import unittest
import tempfile
import numpy as np
import pandas as pd
from pathlib import Path
from pathoscope.core.expression import ExpressionAnalyzer
from pathoscope.utils.config_loader import AppConfig

class TestExpressionAnalyzer(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.outdir = Path(self.temp_dir.name)
        self.analyzer = ExpressionAnalyzer()
        
        # Construct raw mock replicated RNA-Seq count matrix
        self.mock_counts = pd.DataFrame({
            "Gene": ["IL6", "IFNB1", "TP53", "GAPDH", "ACTB"],
            "Control_1": [15, 5, 120, 2000, 2500],
            "Control_2": [20, 7, 115, 2100, 2400],
            "Virus_1": [550, 300, 130, 1950, 2600],
            "Virus_2": [600, 340, 125, 2050, 2550]
        })
        
        self.counts_csv = Path(self.temp_dir.name) / "mock_counts.csv"
        self.mock_counts.to_csv(self.counts_csv, index=False)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_cpm_stabilizing_normalization(self):
        """Assert that Counts Per Million normalization and log2 transformations stabilize values."""
        df_norm = self.analyzer.cpm_normalize_and_log_transform(self.mock_counts)
        
        self.assertEqual(len(df_norm), 5)
        # Replicates columns must exist
        self.assertIn("Control_1", df_norm.columns)
        self.assertIn("Virus_1", df_norm.columns)
        
        # Normalized values should reflect log scaling (e.g. not original integers)
        self.assertAlmostEqual(df_norm.loc[df_norm["Gene"] == "GAPDH", "Control_1"].values[0], np.log2(((2000 + 1) / 4640) * 1e6), places=2)

    def test_welch_differential_sweeps(self):
        """Assert that Welch's t-test accurately computes raw fold changes and p-values."""
        df_norm = self.analyzer.cpm_normalize_and_log_transform(self.mock_counts)
        
        controls = ["Control_1", "Control_2"]
        treated = ["Virus_1", "Virus_2"]
        
        df_deg = self.analyzer.perform_welch_t_test(df_norm, controls, treated)
        
        self.assertEqual(len(df_deg), 5)
        self.assertIn("log2_fold_change", df_deg.columns)
        self.assertIn("raw_pvalue", df_deg.columns)
        
        # IL6 must be highly upregulated (fold change > 0)
        il6_fc = df_deg.loc[df_deg["Gene"] == "IL6", "log2_fold_change"].values[0]
        self.assertGreater(il6_fc, 3.0)

    def test_hybrid_expression_orchestrator(self):
        """Assert that the high-level workflow dispatcher runs end-to-end successfully."""
        config = AppConfig()
        config.statistics.log2_fc_threshold = 1.5
        config.statistics.fdr_threshold = 0.05
        
        controls = ["Control_1", "Control_2"]
        treated = ["Virus_1", "Virus_2"]
        
        summary = self.analyzer.run_hybrid_expression_analysis(
            self.counts_csv, self.outdir, config, controls, treated
        )
        
        # Verify complete run metrics
        self.assertEqual(summary["counts"]["total_genes_ingested"], 5)
        self.assertEqual(summary["counts"]["upregulated_genes_detected"], 2) # IL6 and IFNB1 should be upregulated
        self.assertEqual(summary["counts"]["downregulated_genes_detected"], 0)
        
        # Check files exist
        self.assertTrue(Path(summary["output_files"]["complete_differential_expression_csv"]).exists())
        self.assertTrue(Path(summary["output_files"]["upregulated_subset_csv"]).exists())

if __name__ == "__main__":
    unittest.main()
