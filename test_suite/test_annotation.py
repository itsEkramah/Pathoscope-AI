import unittest
from pathlib import Path
from pathoscope.core.annotator import smith_waterman_local_align, calculate_annotation_confidence, parse_blast_tabular, AlignmentHit

class TestAnnotation(unittest.TestCase):
    """
    Bioinformatics QA Test Case: Evaluates dynamic programming homology metrics, E-values,
    query coverage %, and asserts the correctness of the Smith-Waterman BLOSUM62 alignment algorithm.
    """
    def test_smith_waterman_local_alignment(self):
        # Alignment of two similar protein sequence fragments
        seq1 = "HEEFGTEI"
        seq2 = "HEEWGTEI"
        
        score, identity, align_len, al1, al2 = smith_waterman_local_align(seq1, seq2)
        
        # Identity: 7 out of 8 residues are identical (87.5%)
        self.assertAlmostEqual(identity, 87.5, places=1)
        self.assertEqual(align_len, 8)
        self.assertEqual(al1, seq1)
        self.assertEqual(al2, seq2)
        # Score should be high (BLOSUM62 match scores are positive, substitution is small)
        self.assertGreater(score, 0.0)

    def test_smith_waterman_gapped_alignment(self):
        # Local alignment with gap insertions
        seq1 = "HEEFGTEI"
        seq2 = "HEEGTEI" # missing F
        
        score, identity, align_len, al1, al2 = smith_waterman_local_align(seq1, seq2)
        
        self.assertGreater(score, 0.0)
        self.assertIn("-", al2) # Gap inserted in sequence 2
        self.assertEqual(align_len, 8)

    def test_annotation_confidence_score(self):
        # 1. High confidence hit: E-value 0.0, 100% identity, 100% query coverage
        conf_high = calculate_annotation_confidence(evalue=0.0, identity=100.0, query_coverage=100.0)
        self.assertEqual(conf_high, 1.0)
        
        # 2. Mid confidence hit: E-value 1e-10, 50% identity, 80% query coverage
        conf_mid = calculate_annotation_confidence(evalue=1e-10, identity=50.0, query_coverage=80.0)
        # evalue_score = 10 / 100 = 0.1
        # Combined: (0.4 * 0.1) + (0.3 * 0.5) + (0.3 * 0.8) = 0.04 + 0.15 + 0.24 = 0.43
        self.assertAlmostEqual(conf_mid, 0.43, places=2)

    def test_alignment_hit_coverage(self):
        # standard tabular columns: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen
        fields = [
            "query_orf_1", "sp|P12345|VIRAL", "80.0", "100", "20", "0", "10", "109", "1", "100", "1e-25", "180.0", "200", "200"
        ]
        hit = AlignmentHit(fields)
        
        # Aligned query length = abs(109 - 10) + 1 = 100bp
        # Query coverage % = (100 / 200) * 100.0 = 50.0%
        self.assertAlmostEqual(hit.query_coverage, 50.0, places=1)
        self.assertAlmostEqual(hit.subject_coverage, 50.0, places=1)
        self.assertEqual(hit.evalue, 1e-25)
        self.assertEqual(hit.bitscore, 180.0)
