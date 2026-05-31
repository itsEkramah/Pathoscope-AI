import unittest
from pathlib import Path
import tempfile
import pandas as pd
from pathoscope.core.orf_predictor import ORFRecord
from pathoscope.reporting.reporter import run_report_generation

class TestReportGeneration(unittest.TestCase):
    """
    Bioinformatics QA Test Case: Exercises the reporter pipeline, verifying that
    GFF3 coordinate rows, CSV tables, and interactive MultiQC dashboards render cleanly.
    """
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_gff3_line_formatting(self):
        # 1-based inclusive forward strand coordinate ORF
        orf = ORFRecord(
            orf_id="ORF_001",
            sequence_id="SARS_CoV_2",
            start=266,
            end=21555,
            strand="+",
            frame=1,
            nucleotide_seq="ATG" + ("NNN" * 7095) + "TAA", # 21290bp
            protein_seq="M" + ("X" * 7095),
            start_codon="ATG",
            confidence_score=0.95
        )
        
        gff_line = orf.to_gff_line()
        fields = gff_line.split("\t")
        
        # GFF3 fields: seqid, source, type, start, end, score, strand, phase, attributes
        self.assertEqual(len(fields), 9)
        self.assertEqual(fields[0], "SARS_CoV_2")
        self.assertEqual(fields[1], "PathoScope_ORF")
        self.assertEqual(fields[2], "CDS")
        self.assertEqual(fields[3], "266")
        self.assertEqual(fields[4], "21555")
        self.assertEqual(fields[5], "0.95")
        self.assertEqual(fields[6], "+")
        self.assertEqual(fields[7], "0")
        self.assertIn("ID=ORF_001", fields[8])
        self.assertIn("confidence_score=0.95", fields[8])

    def test_report_generation_process(self):
        # Set up mock folders matching Pipeline output layout
        pre_dir = self.output_dir / "preprocessed"
        orf_dir = self.output_dir / "orfs"
        ann_dir = self.output_dir / "annotations"
        path_dir = self.output_dir / "pathways"
        enr_dir = self.output_dir / "enrichment"
        vis_dir = self.output_dir / "visualizations"
        ai_dir = self.output_dir / "ai_interpretations"
        
        for d in [pre_dir, orf_dir, ann_dir, path_dir, enr_dir, vis_dir, ai_dir]:
            d.mkdir(parents=True)
            
        # Write dummy intermediate files
        with open(pre_dir / "qc_report.json", "w") as f:
            f.write('{"metrics": {"n50": 500, "mean_gc_percent": 40.0}, "counts": {"total_processed": 5, "total_kept": 4}}')
            
        with open(orf_dir / "coordinates.gff3", "w") as f:
            f.write("SARS_CoV_2\tPathoScope_ORF\tCDS\t10\t100\t0.9\t+\t0\tID=ORF_001\n")
            
        df_annot = pd.DataFrame([{
            "protein_id": "ORF_001",
            "uniprot_id": "P12345",
            "subject_db_id": "sp|P12345|VIRAL",
            "identity_percent": 80.0,
            "evalue": 1e-25,
            "annotation_status": "Annotated",
            "annotation_confidence": 0.95
        }])
        df_annot.to_csv(ann_dir / "annotated_proteins.csv", index=False)
        
        df_paths = pd.DataFrame([{
            "protein_id": "ORF_001",
            "pathway_id": "map03010",
            "pathway_description": "Ribosome",
            "source_database": "KEGG",
            "pfam_validated": 1
        }])
        df_paths.to_csv(path_dir / "mapped_pathways.csv", index=False)
        
        df_enrich = pd.DataFrame([{
            "pathway_id": "map03010",
            "description": "Ribosome",
            "raw_pvalue": 0.001,
            "adjusted_pvalue_fdr": 0.005,
            "fold_enrichment": 5.0,
            "query_count_k": 2,
            "pathway_background_M": 50
        }])
        df_enrich.to_csv(enr_dir / "significant_pathways.csv", index=False)
        
        # Execute report builder
        report_stats = run_report_generation(self.output_dir)
        
        # Verify dashboard reports created
        self.assertTrue(Path(report_stats["html_report"]).exists())
        self.assertGreater(Path(report_stats["html_report"]).stat().st_size, 0)
