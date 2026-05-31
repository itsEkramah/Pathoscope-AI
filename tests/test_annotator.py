# pyrefly: ignore [missing-import]
import pytest
from pathlib import Path
import tempfile
import pandas as pd
from Bio import SeqIO
from pathoscope.core.preprocessor import SequenceRecord
from pathoscope.core.annotator import (
    AlignmentHit,
    parse_blast_tabular,
    filter_and_rank_hits,
    process_functional_annotation,
    DatabaseNotFoundError,
    ToolExecutionError
)
from pathoscope.utils.config_loader import AppConfig

# 1. Test AlignmentHit parsing and coverage calculations
def test_alignment_hit_calculations():
    # Fields: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen
    fields = [
        "ORF_1", "NP_10239.1", "75.5", "100", "20", "2", 
        "10", "109", "50", "149", "1.2e-12", "180.5", "200", "300"
    ]
    hit = AlignmentHit(fields)
    assert hit.query_id == "ORF_1"
    assert hit.subject_id == "NP_10239.1"
    assert hit.identity == 75.5
    assert hit.evalue == 1.2e-12
    assert hit.bitscore == 180.5
    
    # Query length 200, aligned from 10 to 109 inclusive = 100 bp. Coverage = 50%
    assert hit.query_coverage == 50.0
    # Subject length 300, aligned from 50 to 149 inclusive = 100 bp. Coverage = 33.33%
    assert round(hit.subject_coverage, 2) == 33.33
    
    d = hit.to_dict()
    assert d["identity_percent"] == 75.5
    assert d["query_coverage"] == 50.0

# 2. Test tabular parser
def test_parse_blast_tabular():
    # Write mock tabular TSV
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False, mode="w") as f:
        # qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen
        f.write("ORF_1\tNP_123\t80.0\t100\t10\t1\t1\t100\t1\t100\t1e-20\t200.0\t100\t100\n")
        f.write("# comment line should be skipped\n")
        f.write("ORF_2\tNP_456\t45.0\t80\t30\t2\t10\t89\t20\t99\t1e-5\t120.0\t100\t100\n")
        temp_path = Path(f.name)
        
    try:
        hits = parse_blast_tabular(temp_path)
        assert len(hits) == 2
        assert hits[0].query_id == "ORF_1"
        assert hits[1].query_id == "ORF_2"
        assert hits[0].identity == 80.0
    finally:
        temp_path.unlink()

# 3. Test filtering and priority-based ranking
def test_filter_and_rank_hits():
    # Setup hits for query 'ORF_1'
    # Hit A: E-value 1e-10, score 150, identity 70% (Good)
    # Hit B: E-value 1e-12, score 180, identity 80% (Better E-value, higher score - Top hit)
    # Hit C: E-value 1e-4, score 90, identity 50% (Fails standard threshold evalue=1e-5)
    # Hit D: E-value 1e-15, score 210, identity 20% (Fails identity threshold 30%)
    # Hit E: E-value 1e-15, score 210, identity 70%, qlen=200, aligned 10bp (Fails query coverage 50%)
    
    # fields: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen
    raw_hit_a = AlignmentHit(["ORF_1", "A", "70.0", "100", "0", "0", "1", "100", "1", "100", "1e-10", "150.0", "100", "100"])
    raw_hit_b = AlignmentHit(["ORF_1", "B", "80.0", "100", "0", "0", "1", "100", "1", "100", "1e-12", "180.0", "100", "100"])
    raw_hit_c = AlignmentHit(["ORF_1", "C", "50.0", "100", "0", "0", "1", "100", "1", "100", "1e-4", "90.0", "100", "100"])
    raw_hit_d = AlignmentHit(["ORF_1", "D", "20.0", "100", "0", "0", "1", "100", "1", "100", "1e-15", "210.0", "100", "100"])
    raw_hit_e = AlignmentHit(["ORF_1", "E", "70.0", "10", "0", "0", "1", "10", "1", "10", "1e-15", "210.0", "200", "200"])
    
    raw_hits = [raw_hit_a, raw_hit_b, raw_hit_c, raw_hit_d, raw_hit_e]
    
    accepted, rejected = filter_and_rank_hits(
        raw_hits, evalue_thresh=1e-5, ident_thresh=30.0, cov_thresh=50.0
    )
    
    # 1. Assert only ONE top hit is accepted for the query ORF_1
    assert len(accepted) == 1
    # 2. Assert Hit B is isolated as the top hit (because of 1e-12 evalue vs 1e-10 of Hit A)
    assert accepted[0].subject_id == "B"
    
    # 3. Assert alternative and failing matches are in the rejected list
    assert len(rejected) == 4
    
    rejections_by_id = {hit.subject_id: reason for hit, reason in rejected}
    # Hit A is rejected because it is alternative
    assert "Alternative" in rejections_by_id["A"]
    # Hit C failed evalue
    assert "E-value" in rejections_by_id["C"]
    # Hit D failed identity
    assert "Identity" in rejections_by_id["D"]
    # Hit E failed query coverage
    assert "Coverage" in rejections_by_id["E"]

# 4. Test complete orchestrator loop with local cache loading
def test_process_functional_annotation_workflow():
    # Setup query protein file (3 predicted proteins)
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
        f.write(
            ">ORF_1\nMAPKRILV\n"
            ">ORF_2\nMNKPTLVR\n"
            ">ORF_3\nMYVKILVP\n"
        )
        temp_proteins = Path(f.name)
        
    temp_outdir = Path(tempfile.mkdtemp())
    
    # Write mock cache file: ORF_1 hits NP_100, ORF_2 hits NP_200, ORF_3 has no cache alignments (Hypothetical)
    cached_tsv = temp_outdir / "alignment_results_cache.tsv"
    with open(cached_tsv, "w", encoding="utf-8") as f:
        # fields: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen
        f.write("ORF_1\tNP_100\t95.0\t8\t0\t0\t1\t8\t1\t8\t1e-15\t50.0\t8\t8\n")
        f.write("ORF_2\tNP_200\t80.0\t8\t1\t0\t1\t8\t1\t8\t1e-12\t42.0\t8\t8\n")
        # No lines for ORF_3 -> hypothetical
        
    config = AppConfig()
    config.annotation.alignment_engine = "diamond"
    config.annotation.eval_threshold = 1e-5
    config.annotation.identity_threshold = 30.0
    config.annotation.coverage_threshold = 50.0
    
    try:
        stats = process_functional_annotation(temp_proteins, temp_outdir, config)
        
        # 1. Verify alignment execution bypassed (used cache) and counts are correct
        assert stats["counts"]["total_proteins_queried"] == 3
        assert stats["counts"]["total_annotated"] == 2
        assert stats["counts"]["total_hypothetical"] == 1
        assert stats["metrics"]["annotation_rate_percent"] == 66.67
        
        # 2. Check output csv tables created
        assert (temp_outdir / "annotated_proteins.csv").exists()
        assert (temp_outdir / "top_hits.csv").exists()
        assert (temp_outdir / "rejected_hits.csv").exists()
        assert (temp_outdir / "annotation_report.json").exists()
        
        # 3. Verify complete annotations table mapping
        df_anno = pd.read_csv(temp_outdir / "annotated_proteins.csv")
        assert len(df_anno) == 3
        
        # Row 0: ORF_1 (Annotated)
        assert df_anno.iloc[0]["protein_id"] == "ORF_1"
        assert df_anno.iloc[0]["annotation_status"] == "Annotated"
        assert df_anno.iloc[0]["subject_db_id"] == "NP_100"
        
        # Row 2: ORF_3 (Graceful Hypothetical Protein mapping)
        assert df_anno.iloc[2]["protein_id"] == "ORF_3"
        assert df_anno.iloc[2]["annotation_status"] == "Hypothetical Protein"
        assert pd.isna(df_anno.iloc[2]["subject_db_id"]) or str(df_anno.iloc[2]["subject_db_id"]) == "None"
        assert df_anno.iloc[2]["description"] == "hypothetical protein"
        assert pd.isna(df_anno.iloc[2]["evalue"]) # empty/none
        
    finally:
        if temp_proteins.exists():
            temp_proteins.unlink()
        for p in temp_outdir.iterdir():
            p.unlink()
        temp_outdir.rmdir()


# 5. Test Smith-Waterman Local Alignment (Dynamic Programming)
def test_smith_waterman_local_alignment():
    from pathoscope.core.annotator import smith_waterman_local_align
    
    # Align similar protein sequences
    # "HELLOWORLD" vs "HELLWORLD" (mismatch/deletion of 'O')
    score, ident, align_len, al1, al2 = smith_waterman_local_align("HELLOWORLD", "HELLWORLD")
    
    assert score > 0
    assert ident > 80.0
    assert align_len == 10
    assert al1 == "HELLOWORLD"
    assert al2 == "HELL-WORLD"

# 6. Test Annotation Confidence Scoring
def test_annotation_confidence_scoring():
    from pathoscope.core.annotator import calculate_annotation_confidence
    
    # Highly significant hit: E-value 1e-80, 100% identity, 100% coverage
    conf_high = calculate_annotation_confidence(1e-80, 100.0, 100.0, 1e-5)
    assert 0.8 <= conf_high <= 1.0
    
    # Borderline hit: E-value 1e-5, 40% identity, 60% coverage
    conf_low = calculate_annotation_confidence(1e-5, 40.0, 60.0, 1e-5)
    assert conf_low < conf_high
    assert 0.0 <= conf_low <= 1.0

# 7. Test AlignmentHit SW Refinement & 16-field format
def test_sw_refinement_pipeline():
    from pathoscope.core.annotator import AlignmentHit, smith_waterman_local_align, calculate_annotation_confidence
    
    # Fields: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen, qseq, sseq
    fields = [
        "ORF_1", "NP_10239.1", "75.0", "10", "2", "1", 
        "1", "10", "1", "10", "1.2e-12", "180.5", "10", "10",
        "HELLOWORLD", "HELLWORLD"  # aligned query and subject sequences
    ]
    hit = AlignmentHit(fields)
    assert hit.qseq == "HELLOWORLD"
    assert hit.sseq == "HELLWORLD"
    
    # Perform refinement
    raw_q = hit.qseq.replace("-", "")
    raw_s = hit.sseq.replace("-", "")
    sw_score, sw_ident, sw_len, al1, al2 = smith_waterman_local_align(raw_q, raw_s)
    
    hit.refined_score = sw_score
    hit.refined_identity = sw_ident
    hit.refined_length = sw_len
    hit.is_refined = True
    
    assert hit.refined_identity > 80.0
    assert hit.is_refined is True

# 8. Test annotations.tsv generation in orchestrator
def test_annotations_tsv_generation():
    # Setup standard mock sequence in clean FASTA
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
        f.write(">ORF_1\nMAPKRILV\n")
        temp_proteins = Path(f.name)
        
    temp_outdir = Path(tempfile.mkdtemp())
    
    # Write mock cache file with 16 columns (supported)
    cached_tsv = temp_outdir / "alignment_results_cache.tsv"
    with open(cached_tsv, "w", encoding="utf-8") as f:
        # fields: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen, qseq, sseq
        f.write("ORF_1\tNP_100\t95.0\t8\t0\t0\t1\t8\t1\t8\t1e-15\t50.0\t8\t8\tMAPKRILV\tMAPKRILV\n")
        
    config = AppConfig()
    config.annotation.alignment_engine = "diamond"
    
    try:
        stats = process_functional_annotation(temp_proteins, temp_outdir, config)
        
        # Verify both annotated_proteins.csv and annotations.tsv are created
        csv_path = temp_outdir / "annotated_proteins.csv"
        tsv_path = temp_outdir / "annotations.tsv"
        
        assert csv_path.exists()
        assert tsv_path.exists()
        
        # Read files and verify exact match
        df_csv = pd.read_csv(csv_path)
        df_tsv = pd.read_csv(tsv_path, sep="\t")
        
        assert len(df_csv) == 1
        assert len(df_tsv) == 1
        assert df_csv.iloc[0]["protein_id"] == df_tsv.iloc[0]["protein_id"]
        assert df_csv.iloc[0]["annotation_confidence"] == df_tsv.iloc[0]["annotation_confidence"]
        assert df_csv.iloc[0]["refined_identity"] == df_tsv.iloc[0]["refined_identity"]
        
    finally:
        if temp_proteins.exists():
            temp_proteins.unlink()
        for p in temp_outdir.iterdir():
            p.unlink()
        temp_outdir.rmdir()
