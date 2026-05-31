import pytest
import json
from pathlib import Path
import tempfile
import pandas as pd
from pathoscope.interpretation.ai_interpreter import (
    build_deterministic_context,
    generate_structured_prompt,
    generate_deterministic_rule_based_interpretation,
    run_ai_biological_interpretation
)
from pathoscope.utils.config_loader import AppConfig

# 1. Test deterministic context compilation from upstream files
def test_build_deterministic_context():
    temp_workspace = Path(tempfile.mkdtemp())
    
    # Create required subfolders
    (temp_workspace / "preprocessed").mkdir()
    (temp_workspace / "orfs").mkdir()
    (temp_workspace / "annotations").mkdir()
    (temp_workspace / "pathways").mkdir()
    (temp_workspace / "enrichment").mkdir()
    
    # Write mock QC
    qc_data = {"counts": {"total_processed": 100, "total_kept": 98}, "statistics": {"n50": 1200, "gc_content_pct": 45.5}}
    with open(temp_workspace / "preprocessed" / "qc_report.json", "w") as f:
        json.dump(qc_data, f)
        
    # Write GFF
    gff_content = "NC_001416.1\tPathoScope_ORF\tCDS\t100\t400\t.\t+\t0\tID=ORF_1;frame=0\n"
    with open(temp_workspace / "orfs" / "coordinates.gff3", "w") as f:
        f.write(gff_content)
        
    # Write Annotations
    annot_data = {"protein_id": ["ORF_1"], "uniprot_id": ["P12345"], "identity_pct": [85.5], "query_coverage_pct": [90.0], "e_val": [1e-12], "bitscore": [115.2]}
    pd.DataFrame(annot_data).to_csv(temp_workspace / "annotations" / "annotated_proteins.csv", index=False)
    
    # Write Pathways detected domains
    dom_data = {"protein_id": ["ORF_1"], "domain_name": ["Ribosomal_L1"], "domain_accession": ["PF00825"], "domain_e_value": [1e-6]}
    pd.DataFrame(dom_data).to_csv(temp_workspace / "pathways" / "pfam_domains.csv", index=False)
    
    # Write Enrichment
    enrich_data = {"pathway_id": ["ko03010"], "description": ["Ribosome"], "query_count_k": [1], "fold_enrichment": [50.0], "adjusted_pvalue_fdr": [1e-4]}
    pd.DataFrame(enrich_data).to_csv(temp_workspace / "enrichment" / "enrichment_results.csv", index=False)
    
    try:
        context = build_deterministic_context(temp_workspace)
        
        # Verify correct context parsing
        assert context["qc_metrics"]["total_kept"] == 98
        assert context["orf_metrics"]["total_predicted"] == 1
        assert context["orf_metrics"]["length_aa_range"] == "100-100"  # (400 - 100 + 1) / 3 = 100aa
        assert len(context["high_confidence_annotations"]) == 1
        assert context["high_confidence_annotations"][0]["uniprot_id"] == "P12345"
        assert len(context["conserved_domains"]) == 1
        assert context["conserved_domains"][0]["domain_name"] == "Ribosomal_L1"
        assert len(context["enriched_pathways"]) == 1
        assert context["enriched_pathways"][0]["pathway_id"] == "ko03010"
        
    finally:
        # Cleanup
        for folder in ["preprocessed", "orfs", "annotations", "pathways", "enrichment"]:
            fdir = temp_workspace / folder
            for p in fdir.iterdir():
                p.unlink()
            fdir.rmdir()
        temp_workspace.rmdir()


# 2. Test prompt generation security parameters
def test_generate_structured_prompt():
    context = {
        "qc_metrics": {"total_processed": 50, "total_kept": 48},
        "orf_metrics": {"total_predicted": 2, "length_aa_range": "30-150"},
        "high_confidence_annotations": [
            {"protein_id": "ORF_1", "uniprot_id": "P54321", "description": "viral protein", "identity_percent": 90.0, "query_coverage_percent": 95.0, "e_value": 1e-15, "bitscore": 140.0}
        ],
        "conserved_domains": [
            {"protein_id": "ORF_1", "domain_name": "Viral_polymerase", "domain_accession": "PF00680", "e_value": 1e-10}
        ],
        "enriched_pathways": [
            {"pathway_id": "map03010", "description": "Ribosome", "query_count_k": 1, "fold_enrichment": 25.0, "adjusted_pvalue_fdr": 0.01}
        ]
    }
    mock_literature = [
        {
            "pmid": "16262622",
            "title": "Viral RNA polymerases: structure, function, and therapeutic targeting",
            "authors": "Castro C, Arnold JJ, Cameron CE",
            "journal": "Virus Research",
            "abstract": "Replication of RNA viruses requires highly conserved RNA-dependent RNA polymerases (RdRp). This study outlines conserved sequence motifs (A-E), dynamic structural borders, and the development of nucleoside analog inhibitors targeting viral replication, establishing RdRp as a primary target for broad-spectrum antiviral design."
        }
    ]
    prompt = generate_structured_prompt(context, mock_literature)
    
    assert "*** CRITICAL RULES ***" in prompt
    assert "Do NOT invent biological findings" in prompt
    assert "cite relevant literature by tracking PMIDs" in prompt
    assert "P54321" in prompt
    assert "PF00680" in prompt
    assert "map03010" in prompt


# 3. Test offline fallback rule-based summarization (deterministic biological synthesis)
def test_generate_deterministic_rule_based_interpretation():
    context = {
        "qc_metrics": {"total_processed": 5, "total_kept": 5, "gc_percent": 38.2, "n50_bp": 500},
        "orf_metrics": {"total_predicted": 2, "length_aa_range": "30-150", "average_length_aa": 90},
        "high_confidence_annotations": [
            {"protein_id": "ORF_1", "uniprot_id": "Q12345", "description": "viral protein", "identity_percent": 35.5, "query_coverage_percent": 60.0, "e_value": 1e-6, "bitscore": 75.0} # Low identity homolog
        ],
        "conserved_domains": [
            {"protein_id": "ORF_1", "domain_name": "Helicase", "domain_accession": "PF00270", "e_value": 1e-4}
        ],
        "enriched_pathways": [
            {"pathway_id": "ko03010", "description": "Ribosome", "query_count_k": 1, "fold_enrichment": 10.0, "adjusted_pvalue_fdr": 0.03}
        ]
    }
    mock_literature = [
        {
            "pmid": "16262622",
            "title": "Viral RNA polymerases: structure, function, and therapeutic targeting",
            "authors": "Castro C, Arnold JJ, Cameron CE",
            "journal": "Virus Research",
            "abstract": "Replication of RNA viruses requires highly conserved RNA-dependent RNA polymerases (RdRp). This study outlines conserved sequence motifs (A-E), dynamic structural borders, and the development of nucleoside analog inhibitors targeting viral replication, establishing RdRp as a primary target for broad-spectrum antiviral design."
        }
    ]
    report = generate_deterministic_rule_based_interpretation(context, mock_literature)
    
    # Verify that all structured keys required by the virology schema exist and are correctly populated
    assert "concise_summary" in report
    assert "detailed_biological_interpretation" in report
    assert "disease_association_summary" in report
    assert "pathway_significance_discussion" in report
    assert "therapeutic_relevance_summary" in report
    assert "limitations" in report
    assert "confidence_warnings" in report
    
    assert "38.2%" in report["detailed_biological_interpretation"]
    assert "Helicase" in report["detailed_biological_interpretation"]
    assert "Ribosome" in report["pathway_significance_discussion"]
    
    # Verify confidence warning generated for low identity alignment
    low_id_warning = [w for w in report["confidence_warnings"] if "Low sequence homology" in w]
    assert len(low_id_warning) == 1
    assert "Q12345" in low_id_warning[0]


# 4. Test run_ai_biological_interpretation orchestrator and fallback integration
def test_run_ai_biological_interpretation_orchestrator():
    temp_workspace = Path(tempfile.mkdtemp())
    (temp_workspace / "preprocessed").mkdir()
    (temp_workspace / "orfs").mkdir()
    (temp_workspace / "annotations").mkdir()
    (temp_workspace / "pathways").mkdir()
    (temp_workspace / "enrichment").mkdir()
    (temp_workspace / "final_report").mkdir()
    
    # Write mock files
    qc_data = {"counts": {"total_processed": 10, "total_kept": 10}, "statistics": {"n50": 300, "gc_content_pct": 40.0}}
    with open(temp_workspace / "preprocessed" / "qc_report.json", "w") as f:
        json.dump(qc_data, f)
        
    gff_content = "NC_001416.1\tPathoScope_ORF\tCDS\t1\t90\t.\t+\t0\tID=ORF_1;frame=0\n"
    with open(temp_workspace / "orfs" / "coordinates.gff3", "w") as f:
        f.write(gff_content)
        
    pd.DataFrame(columns=["protein_id", "uniprot_id", "identity_pct", "query_coverage_pct", "e_val", "bitscore"]).to_csv(temp_workspace / "annotations" / "annotated_proteins.csv", index=False)
    pd.DataFrame(columns=["protein_id", "domain_name", "domain_accession", "domain_e_value"]).to_csv(temp_workspace / "pathways" / "pfam_domains.csv", index=False)
    pd.DataFrame(columns=["pathway_id", "description", "query_count_k", "fold_enrichment", "adjusted_pvalue_fdr"]).to_csv(temp_workspace / "enrichment" / "enrichment_results.csv", index=False)
    
    config = AppConfig()
    config.ai_interpretation.provider = "unsupported_endpoint" # Forces rule-based fallback
    
    try:
        report_meta = run_ai_biological_interpretation(temp_workspace, config)
        
        assert report_meta["status"] == "SUCCESS"
        assert (temp_workspace / "final_report" / "ai_synthesis.json").exists()
        
        # Verify synthesized JSON content
        with open(temp_workspace / "final_report" / "ai_synthesis.json", "r") as f:
            saved_report = json.load(f)
            assert "concise_summary" in saved_report
            assert "detailed_biological_interpretation" in saved_report
            assert "limitations" in saved_report
            
    finally:
        # Cleanup
        for folder in ["preprocessed", "orfs", "annotations", "pathways", "enrichment", "final_report"]:
            fdir = temp_workspace / folder
            for p in fdir.iterdir():
                p.unlink()
            fdir.rmdir()
        temp_workspace.rmdir()


# 5. Test PubMed literature retriever offline & online patterns
from pathoscope.interpretation.ai_interpreter import PubMedLiteratureRetriever

def test_pubmed_literature_retriever():
    retriever = PubMedLiteratureRetriever()
    
    # Test fallback articles matching
    articles = retriever.get_fallback_articles(["MS2"])
    assert len(articles) > 0
    assert any("MS2" in a["title"] or "MS2" in a["abstract"] for a in articles)
    
    # Test fallback articles matching Polymerase
    articles_poly = retriever.get_fallback_articles(["polymerase"])
    assert len(articles_poly) > 0
    assert any("polymerase" in a["title"].lower() or "polymerase" in a["abstract"].lower() for a in articles_poly)
    
    # Test retrieve_pubmed_literature with empty terms
    default_articles = retriever.retrieve_pubmed_literature([])
    assert len(default_articles) == 1
    assert default_articles[0]["pmid"] == "16262622"

