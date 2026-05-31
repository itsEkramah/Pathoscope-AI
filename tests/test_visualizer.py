# pyrefly: ignore [missing-import]
import pytest
from pathlib import Path
import tempfile
import pandas as pd
import numpy as np
import shutil
from pathoscope.reporting.visualizer import (
    generate_orf_distribution_plots,
    generate_annotation_distribution_plots,
    generate_enrichment_plots,
    generate_volcano_plot,
    generate_enrichment_heatmap,
    generate_pca_plot,
    generate_protein_pathway_network,
    run_all_visualizations
)
from pathoscope.utils.config_loader import AppConfig


# 1. Test ORF coordinate distribution plots
def test_generate_orf_distribution_plots():
    temp_dir = Path(tempfile.mkdtemp())
    gff_path = temp_dir / "coordinates.gff3"
    
    gff_content = (
        "##gff-version 3\n"
        "NC_001416.1\tPathoScope_ORF\tCDS\t100\t500\t.\t+\t0\tID=ORF_1;frame=0\n"
        "NC_001416.1\tPathoScope_ORF\tCDS\t600\t1200\t.\t-\t1\tID=ORF_2;frame=1\n"
        "NC_001416.1\tPathoScope_ORF\tCDS\t1500\t1800\t.\t+\t2\tID=ORF_3;frame=2\n"
    )
    with open(gff_path, "w") as f:
        f.write(gff_content)
        
    try:
        results = generate_orf_distribution_plots(gff_path, temp_dir)
        
        # Verify raster PNG exists
        assert "lengths" in results
        assert "track" in results
        assert Path(results["lengths"]).exists()
        assert Path(results["track"]).exists()
        
        # Verify vector formats exist too
        assert Path(temp_dir / "orf_lengths.svg").exists()
        assert Path(temp_dir / "orf_lengths.pdf").exists()
        assert Path(temp_dir / "orf_genomic_track.svg").exists()
        assert Path(temp_dir / "orf_genomic_track.pdf").exists()
        
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# 2. Test Swiss-Prot sequence annotation metrics distributions
def test_generate_annotation_distribution_plots():
    temp_dir = Path(tempfile.mkdtemp())
    csv_path = temp_dir / "annotated_proteins.csv"
    
    data = {
        "protein_id": ["ORF_1", "ORF_2", "ORF_3"],
        "uniprot_id": ["P12345", "Q67890", "None"],
        "identity_pct": [65.4, 42.1, 0.0],
        "query_coverage_pct": [88.0, 92.5, 0.0],
        "bitscore": [120.5, 95.0, 0.0]
    }
    pd.DataFrame(data).to_csv(csv_path, index=False)
    
    try:
        dist_plot = generate_annotation_distribution_plots(csv_path, temp_dir)
        
        assert dist_plot is not None
        assert Path(dist_plot).exists()
        assert Path(dist_plot).name == "annotation_distributions.png"
        
        # Verify vector formats
        assert Path(temp_dir / "annotation_distributions.svg").exists()
        assert Path(temp_dir / "annotation_distributions.pdf").exists()
        
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# 3. Test pathway enrichment barplots and bubble plots
def test_generate_enrichment_plots():
    temp_dir = Path(tempfile.mkdtemp())
    csv_path = temp_dir / "enrichment_results.csv"
    
    data = {
        "pathway_id": ["ko03010", "ko00240"],
        "description": ["Ribosome", "Pyrimidine metabolism"],
        "source_database": ["KEGG", "KEGG"],
        "query_count_k": [5, 2],
        "query_set_size_n": [10, 10],
        "background_count_M": [20, 15],
        "background_universe_N": [1000, 1000],
        "fold_enrichment": [25.0, 13.33],
        "raw_pvalue": [1e-6, 0.02],
        "adjusted_pvalue_fdr": [2e-6, 0.04]
    }
    pd.DataFrame(data).to_csv(csv_path, index=False)
    
    try:
        results = generate_enrichment_plots(csv_path, temp_dir)
        
        assert "barplot" in results
        assert "bubbleplot" in results
        assert Path(results["barplot"]).exists()
        assert Path(results["bubbleplot"]).exists()
        
        # Verify interactive HTML formats exist
        assert "barplot_html" in results or (temp_dir / "pathway_enrichment_barplot_interactive.html").exists()
        assert "bubbleplot_html" in results or (temp_dir / "pathway_enrichment_bubbleplot_interactive.html").exists()
        
        # Verify vector formats
        assert Path(temp_dir / "pathway_enrichment_barplot.svg").exists()
        assert Path(temp_dir / "pathway_enrichment_barplot.pdf").exists()
        
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# 4. Test Volcano Plot generation
def test_generate_volcano_plot():
    temp_dir = Path(tempfile.mkdtemp())
    csv_path = temp_dir / "volcano_plot_data.csv"
    
    data = {
        "pathway_id": ["ko03010", "ko00240", "ko00010"],
        "description": ["Ribosome", "Pyrimidine metabolism", "Glycolysis"],
        "fold_enrichment": [2.5, 0.8, 1.2],
        "raw_pvalue": [1e-4, 0.4, 0.2],
        "adjusted_pvalue_fdr": [2e-4, 0.5, 0.3],
        "significance": ["SIGNIFICANT", "NOT_SIGNIFICANT", "NOT_SIGNIFICANT"]
    }
    pd.DataFrame(data).to_csv(csv_path, index=False)
    
    try:
        results = volcano_plots = generate_volcano_plot(csv_path, temp_dir)
        
        assert "volcano" in results
        assert Path(results["volcano"]).exists()
        assert Path(results["volcano"]).name == "pathway_volcano_plot.png"
        
        # Verify vector formats
        assert Path(temp_dir / "pathway_volcano_plot.svg").exists()
        assert Path(temp_dir / "pathway_volcano_plot.pdf").exists()
        
        # Verify interactive html exists
        assert (temp_dir / "pathway_volcano_plot_interactive.html").exists()
        
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# 5. Test Gene-Pathway Profile Heatmap
def test_generate_enrichment_heatmap():
    temp_dir = Path(tempfile.mkdtemp())
    matrix_csv = temp_dir / "enrichment_matrix.csv"
    
    # Write mock binary matrix
    data = {
        "pathway_id": ["ko03010", "ko00240"],
        "P12345": [1, 0],
        "Q67890": [1, 1]
    }
    pd.DataFrame(data).set_index("pathway_id").to_csv(matrix_csv)
    
    try:
        heatmap_path = generate_enrichment_heatmap(matrix_csv, temp_dir)
        
        assert heatmap_path is not None
        assert Path(heatmap_path).exists()
        assert Path(heatmap_path).name == "gene_pathway_heatmap.png"
        
        # Verify vector formats
        assert Path(temp_dir / "gene_pathway_heatmap.svg").exists()
        assert Path(temp_dir / "gene_pathway_heatmap.pdf").exists()
        
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# 6. Test native SVD-based PCA projection Plot
def test_generate_pca_plot():
    temp_dir = Path(tempfile.mkdtemp())
    matrix_csv = temp_dir / "enrichment_matrix.csv"
    
    # PCA needs a minimum dimension of (3 genes x 2 pathways)
    data = {
        "uniprot_id": ["P1", "P2", "P3"],
        "ko1": [1, 0, 1],
        "ko2": [0, 1, 1]
    }
    pd.DataFrame(data).set_index("uniprot_id").to_csv(matrix_csv)
    
    try:
        pca_path = generate_pca_plot(matrix_csv, temp_dir)
        
        assert pca_path is not None
        assert Path(pca_path).exists()
        assert Path(pca_path).name == "functional_pca_plot.png"
        
        # Verify vector formats
        assert Path(temp_dir / "functional_pca_plot.svg").exists()
        assert Path(temp_dir / "functional_pca_plot.pdf").exists()
        assert Path(temp_dir / "functional_pca_plot_interactive.html").exists()
        
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# 7. Test bipartite protein-pathway interaction network
def test_generate_protein_pathway_network():
    temp_dir = Path(tempfile.mkdtemp())
    mappings_csv = temp_dir / "mapped_pathways.csv"
    enrich_csv = temp_dir / "enrichment_results.csv"
    
    map_data = {
        "protein_id": ["ORF_1", "ORF_1_dup", "ORF_2"],
        "uniprot_id": ["P12345", "P12345", "Q67890"],
        "pathway_id": ["ko03010", "ko03010", "ko00240"],
        "pathway_description": ["Ribosome", "Ribosome", "Pyrimidine metabolism"]
    }
    pd.DataFrame(map_data).to_csv(mappings_csv, index=False)
    
    enrich_data = {
        "pathway_id": ["ko03010", "ko00240"],
        "description": ["Ribosome", "Pyrimidine metabolism"],
        "raw_pvalue": [1e-6, 0.02],
        "adjusted_pvalue_fdr": [2e-6, 0.04]
    }
    pd.DataFrame(enrich_data).to_csv(enrich_csv, index=False)
    
    try:
        net_plot = generate_protein_pathway_network(mappings_csv, enrich_csv, temp_dir)
        
        assert net_plot is not None
        assert Path(net_plot).exists()
        assert Path(net_plot).name == "protein_pathway_network.png"
        
        # Verify vector formats
        assert Path(temp_dir / "protein_pathway_network.svg").exists()
        assert Path(temp_dir / "protein_pathway_network.pdf").exists()
        
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# 8. Test complete run_all_visualizations orchestrator
def test_run_all_visualizations():
    temp_workspace = Path(tempfile.mkdtemp())
    
    # Create required subfolders
    (temp_workspace / "orfs").mkdir()
    (temp_workspace / "annotations").mkdir()
    (temp_workspace / "pathways").mkdir()
    (temp_workspace / "enrichment").mkdir()
    
    # GFF
    gff_content = "NC_001416.1\tPathoScope_ORF\tCDS\t10\t100\t.\t+\t0\tID=ORF_1;frame=0\n"
    with open(temp_workspace / "orfs" / "coordinates.gff3", "w") as f:
        f.write(gff_content)
        
    # Annotations
    annot_data = {"protein_id": ["ORF_1"], "uniprot_id": ["P12345"], "identity_pct": [90], "query_coverage_pct": [95], "bitscore": [150]}
    pd.DataFrame(annot_data).to_csv(temp_workspace / "annotations" / "annotated_proteins.csv", index=False)
    
    # Mappings
    map_data = {"protein_id": ["ORF_1"], "uniprot_id": ["P12345"], "pathway_id": ["ko03010"], "pathway_description": ["Ribosome"]}
    pd.DataFrame(map_data).to_csv(temp_workspace / "pathways" / "mapped_pathways.csv", index=False)
    
    # Enrichment
    enrich_data = {
        "pathway_id": ["ko03010"], "description": ["Ribosome"], "query_count_k": [1], 
        "query_set_size_n": [1], "background_count_M": [10], "background_universe_N": [1000], 
        "fold_enrichment": [100], "raw_pvalue": [1e-4], "adjusted_pvalue_fdr": [1e-4]
    }
    pd.DataFrame(enrich_data).to_csv(temp_workspace / "enrichment" / "enrichment_results.csv", index=False)
    
    # Volcano Plot
    volc_data = {
        "pathway_id": ["ko03010"], "description": ["Ribosome"], "fold_enrichment": [100.0], 
        "raw_pvalue": [1e-4], "adjusted_pvalue_fdr": [1e-4], "significance": ["SIGNIFICANT"]
    }
    pd.DataFrame(volc_data).to_csv(temp_workspace / "enrichment" / "volcano_plot_data.csv", index=False)
    
    # Enrichment Matrix
    matrix_data = {
        "pathway_id": ["ko03010", "ko1", "ko2"],
        "P12345": [1, 0, 1],
        "P67890": [0, 1, 1],
        "P99999": [1, 1, 0]
    }
    pd.DataFrame(matrix_data).set_index("pathway_id").to_csv(temp_workspace / "enrichment" / "enrichment_matrix.csv")
    
    config = AppConfig()
    
    try:
        report = run_all_visualizations(temp_workspace, config)
        
        assert report["status"] == "SUCCESS"
        assert (temp_workspace / "visualizations" / "orf_lengths.png").exists()
        assert (temp_workspace / "visualizations" / "orf_genomic_track.png").exists()
        assert (temp_workspace / "visualizations" / "annotation_distributions.png").exists()
        assert (temp_workspace / "visualizations" / "pathway_enrichment_barplot.png").exists()
        assert (temp_workspace / "visualizations" / "pathway_enrichment_bubbleplot.png").exists()
        assert (temp_workspace / "visualizations" / "protein_pathway_network.png").exists()
        assert (temp_workspace / "visualizations" / "pathway_volcano_plot.png").exists()
        assert (temp_workspace / "visualizations" / "gene_pathway_heatmap.png").exists()
        assert (temp_workspace / "visualizations" / "functional_pca_plot.png").exists()
        
    finally:
        shutil.rmtree(temp_workspace, ignore_errors=True)
