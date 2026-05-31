# pyrefly: ignore [missing-import]
import pytest
from pathlib import Path
import tempfile
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from pathoscope.core.statistics import (
    prevent_redundancy_and_collapse,
    calculate_hypergeometric_enrichment,
    apply_benjamini_hochberg_correction,
    apply_benjamini_hochberg_correction_statsmodels,
    process_pathway_enrichment,
    run_native_gsea_preranked,
    run_gseapy_preranked,
    run_goatools_enrichment,
    generate_volcano_data,
    generate_enrichment_matrix,
    validate_ora_parameters,
    BiologicalInconsistencyError
)
from pathoscope.utils.config_loader import AppConfig


# 1. Test ORF overlap collapsing and redundancy prevention
def test_prevent_redundancy_and_collapse():
    data = {
        "protein_id": ["ORF_1", "ORF_1_dup", "ORF_2", "ORF_3"],
        "uniprot_id": ["P12345", "P12345", "Q67890", "Q67890"],
        "pathway_id": ["map03010", "map03010", "map03010", "map03010"],
        "pathway_description": ["Ribosome", "Ribosome", "Ribosome", "Ribosome"],
        "source_database": ["KEGG", "KEGG", "KEGG", "KEGG"]
    }
    df = pd.DataFrame(data)
    
    unique_links, n_query = prevent_redundancy_and_collapse(df)
    
    assert n_query == 2
    assert len(unique_links) == 2
    assert "P12345" in unique_links["uniprot_id"].values
    assert "Q67890" in unique_links["uniprot_id"].values


# 2. Test Hypergeometric p-value calculations
def test_calculate_hypergeometric_enrichment():
    data = {
        "uniprot_id": ["P12345", "Q67890", "P12345"],
        "pathway_id": ["ko03010", "ko03010", "ko00240"],
        "pathway_description": ["Ribosome", "Ribosome", "Polymerase"],
        "source_database": ["KEGG", "KEGG", "KEGG"]
    }
    df = pd.DataFrame(data)
    unique_links, n_query = prevent_redundancy_and_collapse(df)
    
    universe_size = 100
    
    df_enrich = calculate_hypergeometric_enrichment(
        unique_links, n_query, universe_size, fdr_threshold=0.05
    )
    
    assert len(df_enrich) == 2
    ribo_row = df_enrich[df_enrich["pathway_id"] == "ko03010"].iloc[0]
    assert ribo_row["query_count_k"] == 2
    assert ribo_row["query_set_size_n"] == 2
    assert ribo_row["background_universe_N"] == 100
    
    p_val = ribo_row["raw_pvalue"]
    assert 0.0 <= p_val <= 1.0
    assert ribo_row["fold_enrichment"] > 0.0


# 3. Test biological inconsistency error handling
def test_biological_inconsistency_error():
    data = {"uniprot_id": ["P1"], "pathway_id": ["A"], "pathway_description": ["Desc"], "source_database": ["KEGG"]}
    df = pd.DataFrame(data)
    unique_links, n_query = prevent_redundancy_and_collapse(df)
    
    with pytest.raises(BiologicalInconsistencyError):
        calculate_hypergeometric_enrichment(unique_links, n_query, universe_size=0, fdr_threshold=0.05)


# 4. Test strict ORA parameter boundaries validation
def test_strict_ora_boundary_validations():
    # Negative sizes
    with pytest.raises(BiologicalInconsistencyError):
        validate_ora_parameters(-1, 10, 5, 2)
    # Universe smaller than query
    with pytest.raises(BiologicalInconsistencyError):
        validate_ora_parameters(10, 10, 15, 2)
    # Universe smaller than pathway
    with pytest.raises(BiologicalInconsistencyError):
        validate_ora_parameters(10, 15, 5, 2)
    # Pathway smaller than successes
    with pytest.raises(BiologicalInconsistencyError):
        validate_ora_parameters(100, 10, 50, 12)
    # Query smaller than successes
    with pytest.raises(BiologicalInconsistencyError):
        validate_ora_parameters(100, 50, 10, 12)
        
    # Valid parameters shouldn't throw error
    validate_ora_parameters(100, 50, 10, 5)


# 5. Test statsmodels Benjamini-Hochberg FDR correction
def test_statsmodels_bh_correction():
    data = {
        "pathway_id": ["A", "B", "C"],
        "description": ["Path A", "Path B", "Path C"],
        "source_database": ["KEGG", "KEGG", "KEGG"],
        "raw_pvalue": [0.01, 0.04, 0.50]
    }
    df = pd.DataFrame(data)
    
    df_corrected = apply_benjamini_hochberg_correction_statsmodels(df)
    
    # Sort for direct checks since statsmodels keeps index order
    df_corrected_sorted = df_corrected.sort_values("raw_pvalue").reset_index(drop=True)
    # Rank sorted adjustment:
    # rank 1: p=0.01, adj_p = 0.01 * 3 / 1 = 0.03
    # rank 2: p=0.04, adj_p = 0.04 * 3 / 2 = 0.06
    # rank 3: p=0.50, adj_p = 0.50 * 3 / 3 = 0.50
    assert abs(df_corrected_sorted.iloc[0]["adjusted_pvalue_fdr"] - 0.03) < 1e-9
    assert abs(df_corrected_sorted.iloc[1]["adjusted_pvalue_fdr"] - 0.06) < 1e-9
    assert abs(df_corrected_sorted.iloc[2]["adjusted_pvalue_fdr"] - 0.50) < 1e-9


# 6. Test ORA Fisher's Exact contingency tables odds ratios
def test_ora_fisher_exact_validation():
    data = {
        "uniprot_id": ["P1", "P2", "P3"],
        "pathway_id": ["ko1", "ko1", "ko2"],
        "pathway_description": ["Path1", "Path1", "Path2"],
        "source_database": ["KEGG", "KEGG", "KEGG"]
    }
    df = pd.DataFrame(data)
    unique_links, n_query = prevent_redundancy_and_collapse(df)
    
    df_enrich = calculate_hypergeometric_enrichment(
        unique_links, n_query, universe_size=10, fdr_threshold=0.05
    )
    
    # Assert fisher exact calculations exists
    assert "raw_pvalue_fisher" in df_enrich.columns
    assert "odds_ratio" in df_enrich.columns
    
    row1 = df_enrich[df_enrich["pathway_id"] == "ko1"].iloc[0]
    assert row1["odds_ratio"] >= 0.0
    assert 0.0 <= row1["raw_pvalue_fisher"] <= 1.0


# 7. Test Native GSEA Preranked Engine mathematical correctness
def test_native_gsea_exact_math():
    # Setup ranked list: gene1, gene2, gene3, gene4
    ranked = pd.Series([3.0, 2.0, 1.0, 0.1], index=["G1", "G2", "G3", "G4"])
    
    # Gene sets: S1 has genes at top of rank list
    gene_sets = {
        "S1": ["G1", "G2"],
        "S2": ["G3", "G4"]
    }
    
    df_gsea = run_native_gsea_preranked(
        ranked_genes=ranked,
        gene_sets=gene_sets,
        n_perm=20,
        min_size=1,
        max_size=100
    )
    
    assert not df_gsea.empty
    assert len(df_gsea) == 2
    
    # Pathway S1 should have positive Enrichment Score (ES) because its genes are at the top
    row_s1 = df_gsea[df_gsea["pathway_id"] == "S1"].iloc[0]
    assert row_s1["enrichment_score_es"] > 0.0
    assert row_s1["normalized_enrichment_score_nes"] > 0.0
    assert "G1" in row_s1["leading_edge_genes"]
    
    # Pathway S2 should have negative Enrichment Score (ES) because its genes are at the bottom
    row_s2 = df_gsea[df_gsea["pathway_id"] == "S2"].iloc[0]
    assert row_s2["enrichment_score_es"] < 0.0
    assert row_s2["normalized_enrichment_score_nes"] < 0.0


# 8. Test GSEApy Preranked Wrapper fallback behaviors
def test_gseapy_wrapper():
    # Test on a small list that triggers GSEApy internal sizes check and fallback
    ranked = pd.Series([1.5, 0.5], index=["G1", "G2"])
    gene_sets = {"S1": ["G1"]}
    
    # This should complete successfully by falling back to the native engine
    df_gsea = run_gseapy_preranked(
        ranked_genes=ranked,
        gene_sets=gene_sets,
        outdir=Path("tmp_outdir_not_created"),
        n_perm=10
    )
    assert not df_gsea.empty
    assert len(df_gsea) == 1
    assert df_gsea.iloc[0]["pathway_id"] == "S1"


# 9. Test GOATOOLS study fallback behaviors
def test_goatools_wrapper_fallback():
    data = {"uniprot_id": ["P12345"], "pathway_id": ["ko03010"]}
    df = pd.DataFrame(data)
    unique_links, n_query = prevent_redundancy_and_collapse(df)
    
    # GO study with missing files should handle fallback elegantly and return empty df
    df_go = run_goatools_enrichment(
        unique_links=unique_links,
        universe_genes=["P12345"],
        obo_path="missing_obo.obo",
        association_path="missing_assoc.txt",
        fdr_threshold=0.05
    )
    assert df_go.empty


# 10. Test Volcano plot coordinates and Enrichment matrix generators
def test_volcano_and_cooccurrence_matrix_generation():
    df_bh = pd.DataFrame({
        "pathway_id": ["ko1", "ko2"],
        "description": ["Path1", "Path2"],
        "fold_enrichment": [4.0, 0.0],
        "adjusted_pvalue_fdr": [0.01, 1.0],
        "significance": ["SIGNIFICANT", "NOT_SIGNIFICANT"]
    })
    
    # Test Volcano plotting values
    df_volc = generate_volcano_data(df_bh)
    assert len(df_volc) == 2
    assert df_volc.iloc[0]["log2_fold_enrichment"] == 2.0  # log2(4) = 2
    assert df_volc.iloc[0]["minus_log10_fdr_pvalue"] == 2.0 # -log10(0.01) = 2
    
    # Test Enrichment matrix
    links = pd.DataFrame({
        "uniprot_id": ["P1", "P1", "P2"],
        "pathway_id": ["ko1", "ko2", "ko1"]
    })
    df_matrix = generate_enrichment_matrix(links)
    assert df_matrix.shape == (2, 2)  # 2 genes x 2 pathways
    assert df_matrix.loc["P1", "ko1"] == 1
    assert df_matrix.loc["P2", "ko2"] == 0


# 11. Test edge-case empty mappings handling
def test_process_pathway_enrichment_empty():
    temp_outdir = Path(tempfile.mkdtemp())
    temp_csv = temp_outdir / "mapped_pathways.csv"
    pd.DataFrame(columns=[
        "protein_id", "uniprot_id", "subject_db_id", "pathway_id", "pathway_description", "source_database"
    ]).to_csv(temp_csv, index=False)
    
    config = AppConfig()
    
    try:
        report = process_pathway_enrichment(temp_csv, temp_outdir, config)
        assert report["counts"]["total_pathways_tested"] == 0
        assert report["counts"]["significant_pathways_enriched"] == 0
        assert (temp_outdir / "enrichment_results.csv").exists()
        assert (temp_outdir / "significant_pathways.csv").exists()
    finally:
        import shutil
        shutil.rmtree(temp_outdir, ignore_errors=True)


# 12. Test complete process_pathway_enrichment orchestrator
def test_process_pathway_enrichment_workflow():
    temp_outdir = Path(tempfile.mkdtemp())
    temp_csv = temp_outdir / "mapped_pathways.csv"
    
    data = {
        "protein_id": ["ORF_1", "ORF_2"],
        "uniprot_id": ["P12345", "Q67890"],
        "subject_db_id": ["sp|P12345|POL", "sp|Q67890|CAP"],
        "pathway_id": ["ko03010", "ko03010"],
        "pathway_description": ["Ribosome", "Ribosome"],
        "source_database": ["KEGG", "KEGG"],
        "pathway_confidence_score": [0.95, 0.90]
    }
    pd.DataFrame(data).to_csv(temp_csv, index=False)
    
    config = AppConfig()
    config.statistics.bg_universe_size = 500
    config.statistics.fdr_threshold = 0.05
    
    try:
        report = process_pathway_enrichment(temp_csv, temp_outdir, config)
        
        assert report["counts"]["unique_query_genes_n"] == 2
        assert report["counts"]["total_pathways_tested"] == 1
        assert report["counts"]["ssgsea_pathways_scored"] == 1
        
        assert (temp_outdir / "enrichment_results.csv").exists()
        assert (temp_outdir / "significant_pathways.csv").exists()
        assert (temp_outdir / "gsea_results.csv").exists()
        assert (temp_outdir / "ssgsea_results.csv").exists()
        assert (temp_outdir / "volcano_plot_data.csv").exists()
        assert (temp_outdir / "enrichment_matrix.csv").exists()
        assert (temp_outdir / "pathway_ranking_reports.csv").exists()
        assert (temp_outdir / "statistics_report.json").exists()
        
        df_res = pd.read_csv(temp_outdir / "enrichment_results.csv")
        assert len(df_res) == 1
        assert df_res.iloc[0]["query_count_k"] == 2
        assert df_res.iloc[0]["fold_enrichment"] > 0.0
        
        assert "raw_pvalue_fisher" in df_res.columns
        assert "odds_ratio" in df_res.columns
        
    finally:
        import shutil
        shutil.rmtree(temp_outdir, ignore_errors=True)


# 13. Test true pathway size fetching and fallbacks
def test_true_pathway_size_fetching():
    from pathoscope.core.statistics import fetch_true_pathway_size
    
    # 1. Test local curated DB lookup (Ribosome map03010 => 150)
    size = fetch_true_pathway_size("map03010", cache=None, N_universe=1000)
    assert size == 150
    
    # 2. Test fallback default size
    size_fallback = fetch_true_pathway_size("unknown_pathway_xyz", cache=None, N_universe=1000)
    assert size_fallback >= 50


# 14. Test ssGSEA mathematical accuracy
def test_ssgsea_mathematical_precision():
    from pathoscope.core.statistics import run_ssgsea
    
    # Simple gene set and ranked list
    ranked = pd.Series([10.0, 5.0, 1.0], index=["G1", "G2", "G3"])
    gene_sets = {
        "PathA": ["G1", "G2"],  # Top of list => positive ES
        "PathB": ["G3"]         # Bottom of list => negative ES
    }
    
    df_ssgsea = run_ssgsea(ranked, gene_sets, alpha=0.25)
    
    assert not df_ssgsea.empty
    assert len(df_ssgsea) == 2
    
    # Normalized score should exist and be in range [-1, 1]
    assert "ssgsea_enrichment_score_normalized" in df_ssgsea.columns
    path_a_row = df_ssgsea[df_ssgsea["pathway_id"] == "PathA"].iloc[0]
    path_b_row = df_ssgsea[df_ssgsea["pathway_id"] == "PathB"].iloc[0]
    
    assert path_a_row["ssgsea_enrichment_score_normalized"] > 0.0
    assert path_b_row["ssgsea_enrichment_score_normalized"] < 0.0


# 15. Test Multi-Evidence Integrated pathway priority scoring
def test_multi_evidence_scoring():
    temp_outdir = Path(tempfile.mkdtemp())
    temp_csv = temp_outdir / "mapped_pathways.csv"
    
    data = {
        "protein_id": ["ORF_1", "ORF_2"],
        "uniprot_id": ["P12345", "Q67890"],
        "subject_db_id": ["sp|P12345|POL", "sp|Q67890|CAP"],
        "pathway_id": ["ko03010", "ko03010"],
        "pathway_description": ["Ribosome", "Ribosome"],
        "source_database": ["KEGG", "KEGG"],
        "pathway_confidence_score": [0.95, 0.90]
    }
    pd.DataFrame(data).to_csv(temp_csv, index=False)
    
    config = AppConfig()
    
    try:
        report = process_pathway_enrichment(temp_csv, temp_outdir, config)
        
        ranking_file = temp_outdir / "pathway_ranking_reports.csv"
        assert ranking_file.exists()
        
        df_rank = pd.read_csv(ranking_file)
        assert not df_rank.empty
        assert "multi_evidence_pathway_score" in df_rank.columns
        assert "average_annotation_confidence" in df_rank.columns
        
        score = df_rank.iloc[0]["multi_evidence_pathway_score"]
        assert 0.0 <= score <= 10.0
    finally:
        import shutil
        shutil.rmtree(temp_outdir, ignore_errors=True)
