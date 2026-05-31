import pytest
import json
from pathlib import Path
import tempfile
import pandas as pd
from pathoscope.reporting.reporter import (
    generate_html_report_dashboard,
    run_report_generation
)

# 1. Test full HTML dashboard generation with embedded base64 plots
def test_generate_html_report_dashboard():
    temp_workspace = Path(tempfile.mkdtemp())
    (temp_workspace / "preprocessed").mkdir()
    (temp_workspace / "orfs").mkdir()
    (temp_workspace / "annotations").mkdir()
    (temp_workspace / "pathways").mkdir()
    (temp_workspace / "enrichment").mkdir()
    (temp_workspace / "final_report").mkdir()
    (temp_workspace / "visualizations").mkdir()
    
    # 1. Write mock preprocessor QC
    qc_data = {"counts": {"total_processed": 100, "total_kept": 98}, "statistics": {"n50": 1200, "gc_content_pct": 45.5}}
    with open(temp_workspace / "preprocessed" / "qc_report.json", "w") as f:
        json.dump(qc_data, f)
        
    # 2. Write mock GFF
    gff_content = "NC_001416.1\tPathoScope_ORF\tCDS\t100\t400\t.\t+\t0\tID=ORF_1;frame=0\n"
    with open(temp_workspace / "orfs" / "coordinates.gff3", "w") as f:
        f.write(gff_content)
        
    # 3. Write mock Swiss-Prot homology CSV
    annot_data = {"protein_id": ["ORF_1"], "uniprot_id": ["P12345"], "identity_pct": [85.5], "query_coverage_pct": [90.0], "e_val": [1e-12], "bitscore": [115.2]}
    pd.DataFrame(annot_data).to_csv(temp_workspace / "annotations" / "annotated_proteins.csv", index=False)
    
    # 3.5 Write mock ICTV taxonomy classification CSV
    tax_data = {
        "protein_id": ["ORF_1"],
        "subject_db_id": ["sp|P12345|ORF1_MS2"],
        "uniprot_organism_code": ["MS2"],
        "ictv_class": ["Leleviricetes"],
        "ictv_order": ["Norzivirales"],
        "ictv_family": ["Fiersviridae"],
        "ictv_genus": ["Emanavirus"],
        "baltimore_group": ["Group IV (+ssRNA)"]
    }
    pd.DataFrame(tax_data).to_csv(temp_workspace / "pathways" / "taxonomy_classification.csv", index=False)

    # 4. Write mock enrichment results
    enrich_data = {"pathway_id": ["ko03010"], "description": ["Ribosome"], "query_count_k": [1], "fold_enrichment": [50.0], "adjusted_pvalue_fdr": [1e-4]}
    pd.DataFrame(enrich_data).to_csv(temp_workspace / "enrichment" / "enrichment_results.csv", index=False)
    
    # 5. Write mock visual graphics (mocking small text file as PNG)
    with open(temp_workspace / "visualizations" / "orf_lengths.png", "wb") as f:
        f.write(b"MOCK_PNG_DATA")
        
    # Mock metadata
    metadata = {
        "pipeline": {"version": "1.0.0"},
        "run": {"input_file": "test.fasta", "input_md5": "abc123xyz", "timestamp": "2026-05-24", "status": "SUCCESS"},
        "environment": {"platform": "Windows", "python_version": "3.11.8"},
        "config_snapshot": {
            "reporting": {"theme": "dark"},
            "statistics": {"fdr_threshold": 0.05},
            "preprocessing": {"min_length": 50}
        }
    }
    
    # Mock AI
    ai_synthesis = {
        "concise_summary": "Test concise summary overview of viral isolate.",
        "detailed_biological_interpretation": "Test detailed scientific biological discussion text.",
        "implicated_host_mechanisms": ["Host transcription hijack", "Virion packaging mapping"],
        "disease_association_summary": "Associated clinical diseases host phenotypes.",
        "pathway_significance_discussion": "Ribosome pathway hijacks are highly significant.",
        "therapeutic_relevance_summary": "Broad spectrum polymerase inhibitors represent valid choices.",
        "literature_evidence_summary": "PubMed citations details support these homologs [PMID: 16262622].",
        "limitations": "Pfam search constraints model boundaries.",
        "confidence_warnings": ["Low sequence coverage alignment warnings."],
        "retrieved_literature_citations": [
            {"pmid": "16262622", "title": "Viral RNA polymerases", "authors": "Castro C", "journal": "Virus Research"}
        ]
    }
    
    try:
        html_path = generate_html_report_dashboard(temp_workspace, metadata, ai_synthesis)
        
        assert html_path.exists()
        assert html_path.name == "report.html"
        
        # Verify rendered HTML contents
        with open(html_path, "r", encoding="utf-8") as f:
            html_content = f.read()
            assert "<!DOCTYPE html>" in html_content
            assert "PathoScope<span>AI</span>" in html_content
            assert "abc123xyz" in html_content
            assert "P12345" in html_content
            assert "Ribosome" in html_content
            assert "Host transcription hijack" in html_content
            assert "Fiersviridae" in html_content
            assert "Group IV (+ssRNA)" in html_content
            assert "Dynamic ICTV Taxonomic Lineages" in html_content
            
    finally:
        # Cleanup
        for folder in ["preprocessed", "orfs", "annotations", "pathways", "enrichment", "final_report", "visualizations"]:
            fdir = temp_workspace / folder
            for p in fdir.iterdir():
                p.unlink()
            fdir.rmdir()
        temp_workspace.rmdir()


# 2. Test full run_report_generation orchestrator
def test_run_report_generation():
    temp_workspace = Path(tempfile.mkdtemp())
    (temp_workspace / "final_report").mkdir()
    
    # Write mock metadata.json
    metadata = {
        "pipeline": {"version": "2.0.0"},
        "run": {"input_file": "mock.fasta", "input_md5": "md5hash123", "timestamp": "2026-05-24", "status": "SUCCESS"},
        "environment": {"platform": "Linux", "python_version": "3.11"},
        "config_snapshot": {
            "reporting": {"theme": "light"},
            "statistics": {"fdr_threshold": 0.05},
            "preprocessing": {"min_length": 50}
        }
    }
    with open(temp_workspace / "metadata.json", "w") as f:
        json.dump(metadata, f)
        
    # Write mock ai_synthesis.json
    ai_synthesis = {
        "concise_summary": "Mock concise.",
        "detailed_biological_interpretation": "Mock detailed.",
        "implicated_host_mechanisms": ["Mock systems"],
        "disease_association_summary": "Mock diseases.",
        "pathway_significance_discussion": "Mock pathways.",
        "therapeutic_relevance_summary": "Mock drugs.",
        "literature_evidence_summary": "Mock citations [PMID: 16262622].",
        "limitations": "Mock bounds.",
        "confidence_warnings": ["Mock warning."],
        "retrieved_literature_citations": [
            {"pmid": "16262622", "title": "Viral RNA polymerases", "authors": "Castro C", "journal": "Virus Research"}
        ]
    }
    with open(temp_workspace / "final_report" / "ai_synthesis.json", "w") as f:
        json.dump(ai_synthesis, f)
        
    try:
        report_meta = run_report_generation(temp_workspace)
        
        assert report_meta["status"] == "SUCCESS"
        assert Path(report_meta["html_report"]).exists()
        
    finally:
        # Cleanup
        for folder in ["final_report"]:
            fdir = temp_workspace / folder
            for p in fdir.iterdir():
                p.unlink()
            fdir.rmdir()
        if (temp_workspace / "metadata.json").exists():
            (temp_workspace / "metadata.json").unlink()
        temp_workspace.rmdir()
