"""
test_id_normalization.py

Unit tests for the newly designed Gene ID Normalization module.
Validates mixed ID normalizations (Ensembl, Entrez, official symbols), SQLite cache seeding,
and deduplication checks.
"""

import unittest
import tempfile
import sqlite3
from pathlib import Path
from pathoscope.core.id_normalizer import GeneIDNormalizer

class TestGeneIDNormalizer(unittest.TestCase):
    def setUp(self):
        # Create a temp directory for SQLite DB cache testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.db_path = Path(self.temp_dir.name) / "test_gene_registry.db"
        self.normalizer = GeneIDNormalizer(self.db_path)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_local_dict_lookup(self):
        """Verify that pre-packaged curated synonyms translate correctly."""
        # Ensembl mapping check
        self.assertEqual(self.normalizer.normalize_id("ENSG00000141510", online_fallback=False), "TP53")
        # Entrez mapping check
        self.assertEqual(self.normalizer.normalize_id("6772", online_fallback=False), "STAT1")
        # Alias synonym check
        self.assertEqual(self.normalizer.normalize_id("IL-6", online_fallback=False), "IL6")

    def test_deduplication_and_cleansing(self):
        """Verify that mixed gene arrays normalize and collapse duplicates successfully."""
        mixed_list = ["ENSG00000141510", "7157", "TP53", "6772", "STAT1", "IL-6"]
        # Expected outputs: TP53 (from Ensembl, Entrez, and Symbol), STAT1 (from Entrez and Symbol), and IL6
        normalized = self.normalizer.normalize_gene_list(mixed_list, online_fallback=False)
        
        self.assertEqual(len(normalized), 3)
        self.assertIn("TP53", normalized)
        self.assertIn("STAT1", normalized)
        self.assertIn("IL6", normalized)

    def test_sqlite_mappings_caching(self):
        """Assert that dynamically discovered HGNC API mappings write correctly to SQLite cache."""
        # Insert a custom mapping directly to test cache write
        self.normalizer.cache_mapping("ENSG00000000001", "CUSTOMGENE")
        
        # Query local registry cache
        self.assertEqual(self.normalizer.normalize_id("ENSG00000000001", online_fallback=False), "CUSTOMGENE")

if __name__ == "__main__":
    unittest.main()
