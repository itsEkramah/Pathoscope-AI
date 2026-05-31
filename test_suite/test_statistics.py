import unittest
import numpy as np
import pandas as pd
from pathoscope.core.statistics import (
    validate_ora_parameters,
    BiologicalInconsistencyError,
    apply_benjamini_hochberg_correction,
    run_ssgsea
)

class TestStatistics(unittest.TestCase):
    """
    Bioinformatics QA Test Case: Exercises the statistical models of PathoScope AI,
    asserting parameter checks, Benjamini-Hochberg FDR correction vectors, and ssGSEA running sums.
    """
    def test_ora_parameter_validation(self):
        # 1. Valid case
        try:
            validate_ora_parameters(N=10000, M=100, n=50, k=5)
        except BiologicalInconsistencyError:
            self.fail("Valid ORA parameters raised BiologicalInconsistencyError unexpectedly.")
            
        # 2. Invalid cases
        with self.assertRaises(BiologicalInconsistencyError):
            validate_ora_parameters(N=100, M=200, n=50, k=5) # N < M
            
        with self.assertRaises(BiologicalInconsistencyError):
            validate_ora_parameters(N=10000, M=100, n=50, k=60) # n < k

    def test_benjamini_hochberg_fdr_adjustment(self):
        # Construct raw enrichment DataFrame
        df_enrich = pd.DataFrame([
            {"raw_pvalue": 0.001},
            {"raw_pvalue": 0.01},
            {"raw_pvalue": 0.03},
            {"raw_pvalue": 0.05},
            {"raw_pvalue": 0.1},
            {"raw_pvalue": 0.5}
        ])
        
        df_adj = apply_benjamini_hochberg_correction(df_enrich)
        
        self.assertEqual(len(df_adj), 6)
        self.assertIn("adjusted_pvalue_fdr", df_adj.columns)
        
        # FDR adjusted p-values must be greater than or equal to raw p-values
        for idx, row in df_adj.iterrows():
            self.assertGreaterEqual(row["adjusted_pvalue_fdr"], row["raw_pvalue"])
            
        # Standard BH ranking order checks
        self.assertLess(df_adj.iloc[0]["adjusted_pvalue_fdr"], df_adj.iloc[1]["adjusted_pvalue_fdr"])

    def test_ssgsea_vectorized_engine(self):
        # 5 ranked genes
        ranked = pd.Series([2.5, 1.8, 1.2, -0.5, -2.1], index=["G1", "G2", "G3", "G4", "G5"])
        # Pathway G1, G2
        gene_sets = {
            "pathway_A": ["G1", "G2"],
            "pathway_B": ["G4", "G5"]
        }
        
        df_ssgsea = run_ssgsea(ranked, gene_sets, alpha=0.25)
        
        self.assertEqual(len(df_ssgsea), 2)
        self.assertIn("ssgsea_enrichment_score_raw", df_ssgsea.columns)
        self.assertIn("ssgsea_enrichment_score_normalized", df_ssgsea.columns)
        
        # Pathway A contains highly ranked genes, so its enrichment score should be positive and greater than B
        row_a = df_ssgsea[df_ssgsea["pathway_id"] == "pathway_A"].iloc[0]
        row_b = df_ssgsea[df_ssgsea["pathway_id"] == "pathway_B"].iloc[0]
        self.assertGreater(row_a["ssgsea_enrichment_score_normalized"], row_b["ssgsea_enrichment_score_normalized"])
