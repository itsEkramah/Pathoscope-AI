import os
import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional
import pandas as pd
import numpy as np

# Force headless matplotlib backend to prevent display errors in automated/server environments
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import networkx as nx
from loguru import logger

# Try to import Plotly for premium interactive HTML features
try:
    import plotly.graph_objects as go
    import plotly.express as px
    import plotly.io as pio
    HAS_PLOTLY = True
    logger.info("Plotly interactive visualization suite successfully imported.")
except ImportError:
    HAS_PLOTLY = False
    logger.info("Plotly is not available. Defaulting to vector static matplotlib visualization suite.")

# Okabe-Ito Colorblind-Safe Color System (Scientific Standards)
OKABE_ITO = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermilion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "black": "#000000",
    "light_grey": "#f5f5f5",
    "grey": "#777777",
    "muted_grey": "#cccccc"
}

# Color allocations
STRAND_COLORS = {
    "+": OKABE_ITO["sky_blue"],
    "-": OKABE_ITO["vermilion"]
}

FRAME_COLORS = {
    0: OKABE_ITO["blue"],
    1: OKABE_ITO["bluish_green"],
    2: OKABE_ITO["orange"],
    3: OKABE_ITO["reddish_purple"],
    4: OKABE_ITO["yellow"],
    5: OKABE_ITO["grey"]
}


def setup_matplotlib_style():
    """Applies publication-style visual standards for matplotlib."""
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["font.size"] = 10
    plt.rcParams["axes.labelsize"] = 11
    plt.rcParams["axes.titlesize"] = 12
    plt.rcParams["xtick.labelsize"] = 9
    plt.rcParams["ytick.labelsize"] = 9
    plt.rcParams["legend.fontsize"] = 9
    plt.rcParams["figure.titlesize"] = 14
    plt.rcParams["grid.color"] = "#e0e0e0"
    plt.rcParams["grid.linestyle"] = "--"
    plt.rcParams["grid.linewidth"] = 0.5


def add_reproducibility_metadata(fig, config: Any = None):
    """
    Appends a standardized academic reproducibility stamp to the bottom of the figure.
    """
    pipeline_name = getattr(config.pipeline, "name", "PathoScope AI") if config else "PathoScope AI"
    version = getattr(config.pipeline, "version", "1.0.0") if config else "1.0.0"
    date_str = datetime.date.today().strftime("%Y-%m-%d")
    stamp = f"Pipeline: {pipeline_name} v{version} | Date: {date_str} | License: Academic Free | Auditable Open Science"
    fig.text(0.5, 0.015, stamp, ha="center", va="bottom", fontsize=7.5, color=OKABE_ITO["grey"], style="italic")


def save_publication_figure(fig, path: Path):
    """
    Saves the figure in multiple publication formats: PNG (300 DPI), SVG, and PDF.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save PNG (Raster high-resolution)
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    # Save SVG (Scalable vector web)
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    # Save PDF (Manuscrit print scalable)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")


def generate_enrichment_plots(
    enrichment_csv: Path,
    outdir: Path,
    config: Any = None,
    top_n: int = 15
) -> Dict[str, Path]:
    """
    Generates upgraded Horizontal Barplot and bubble plots using colorblind-safe parameters.
    Also produces Plotly interactive equivalents if available.
    """
    setup_matplotlib_style()
    enrichment_csv = Path(enrichment_csv)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    generated_files = {}
    
    if not enrichment_csv.exists() or os.path.getsize(enrichment_csv) <= 1:
        logger.warning(f"Enrichment CSV file {enrichment_csv} is empty or missing. Skipping plots.")
        return generated_files

    try:
        df = pd.read_csv(enrichment_csv)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as e:
        logger.warning(f"Enrichment CSV file {enrichment_csv} could not be parsed: {e}. Skipping plots.")
        return generated_files

    if df.empty:
        logger.warning("Empty enrichment results. Skipping visualization.")
        return generated_files

    df = df.dropna(subset=["raw_pvalue"]).copy()
    if df.empty:
        logger.warning("No valid p-values in enrichment CSV. Skipping ORA visualizers.")
        return generated_files

    df["neg_log10_p"] = -np.log10(df["raw_pvalue"] + 1e-300)
    df["neg_log10_adj_p"] = -np.log10(df["adjusted_pvalue_fdr"] + 1e-300)
    
    df_sorted = df.sort_values("raw_pvalue").head(top_n)
    df_plot = df_sorted.iloc[::-1] # horizontal reverse order

    # --- Plot 1: Pathway Enrichment Barplot ---
    fig, ax = plt.subplots(figsize=(9, min(6, len(df_plot) * 0.4 + 2)))
    
    # Colorblind safe color scale (Okabe-Ito blending using a linear blue/vermilion transition)
    norm = plt.Normalize(df_plot["fold_enrichment"].min(), df_plot["fold_enrichment"].max())
    cmap = sns.color_palette("cividis", as_cmap=True)
    colors = [cmap(norm(val)) for val in df_plot["fold_enrichment"]]
    
    bars = ax.barh(
        df_plot["description"],
        df_plot["neg_log10_adj_p"],
        color=colors,
        edgecolor=OKABE_ITO["black"],
        linewidth=0.6,
        height=0.6
    )
    
    # Cutoff line (FDR = 0.05 => -log10(0.05) ~ 1.301)
    sig_line = -np.log10(0.05)
    ax.axvline(sig_line, color=OKABE_ITO["vermilion"], linestyle="--", linewidth=1.5, label="FDR Cutoff (0.05)")
    
    ax.set_xlabel("-log10(Adjusted p-value, FDR)")
    ax.set_ylabel("Biological Pathways")
    ax.set_title("Pathway Enrichment Analysis\n(Top Enriched Implicated Pathways)", pad=15)
    ax.grid(True, axis="x")
    
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", pad=0.02, shrink=0.7)
    cbar.set_label("Fold Enrichment", fontsize=9)
    
    ax.legend(loc="lower right")
    add_reproducibility_metadata(fig, config)
    
    barplot_path = outdir / "pathway_enrichment_barplot"
    save_publication_figure(fig, barplot_path)
    plt.close()
    
    generated_files["barplot"] = barplot_path.with_suffix(".png")
    logger.info(f"Saved horizontal enrichment barplot in multiple formats at {barplot_path}.*")

    # --- Plot 2: Pathway Enrichment Bubble Plot ---
    fig, ax = plt.subplots(figsize=(10, min(7, len(df_plot) * 0.5 + 2)))
    
    k_min = df_plot["query_count_k"].min()
    k_max = df_plot["query_count_k"].max()
    sizes = df_plot["query_count_k"].apply(lambda x: 100 + (x - k_min) / max(1, k_max - k_min) * 500 if k_max > k_min else 250)

    sc = ax.scatter(
        df_plot["fold_enrichment"],
        df_plot["description"],
        s=sizes,
        c=df_plot["neg_log10_adj_p"],
        cmap="cividis",
        alpha=0.85,
        edgecolors="black",
        linewidths=0.7
    )
    
    ax.set_xlabel("Fold Enrichment")
    ax.set_ylabel("Biological Pathways")
    ax.set_title("Pathway Enrichment Bubble Plot\n(Dot size = Query Gene Count; Color = -log10 FDR)", pad=15)
    ax.grid(True, axis="both")
    
    cbar = fig.colorbar(sc, ax=ax, pad=0.02, shrink=0.7)
    cbar.set_label("-log10(Adjusted p-value)", fontsize=9)
    
    handles = []
    vals = np.linspace(k_min, k_max, min(3, k_max - k_min + 1)).astype(int)
    for v in vals:
        sz = 100 + (v - k_min) / max(1, k_max - k_min) * 500 if k_max > k_min else 250
        handles.append(
            ax.scatter([], [], s=sz, c="grey", alpha=0.6, edgecolors="black", label=f"{v} query genes")
        )
    ax.legend(handles=handles, loc="lower right", title="Pathway Counts (k)")
    
    add_reproducibility_metadata(fig, config)
    bubbleplot_path = outdir / "pathway_enrichment_bubbleplot"
    save_publication_figure(fig, bubbleplot_path)
    plt.close()
    
    generated_files["bubbleplot"] = bubbleplot_path.with_suffix(".png")
    logger.info(f"Saved enrichment bubbleplot in multiple formats at {bubbleplot_path}.*")

    # --- Interactive Plotly Equivalents ---
    if HAS_PLOTLY:
        try:
            # Barplot
            fig_int = px.bar(
                df_sorted,
                x="neg_log10_adj_p",
                y="description",
                color="fold_enrichment",
                color_continuous_scale="cividis",
                labels={"neg_log10_adj_p": "-log10 FDR", "description": "Pathway", "fold_enrichment": "Fold Enrichment"},
                title="Pathway Enrichment Analysis (Interactive)",
                orientation="h"
            )
            fig_int.add_vline(x=sig_line, line_dash="dash", line_color="red", annotation_text="Cutoff")
            fig_int.update_layout(template="plotly_white", margin=dict(l=20, r=20, t=40, b=40))
            
            bar_html = outdir / "pathway_enrichment_barplot_interactive.html"
            pio.write_html(fig_int, str(bar_html))
            generated_files["barplot_html"] = bar_html
            
            # Bubbleplot
            fig_bub = px.scatter(
                df_sorted,
                x="fold_enrichment",
                y="description",
                size="query_count_k",
                color="neg_log10_adj_p",
                color_continuous_scale="cividis",
                labels={"fold_enrichment": "Fold Enrichment", "neg_log10_adj_p": "-log10 FDR", "query_count_k": "Query Count"},
                title="Pathway Enrichment Bubble Plot (Interactive)"
            )
            fig_bub.update_layout(template="plotly_white")
            bub_html = outdir / "pathway_enrichment_bubbleplot_interactive.html"
            pio.write_html(fig_bub, str(bub_html))
            generated_files["bubbleplot_html"] = bub_html
        except Exception as e:
            logger.warning(f"Plotly bar/bubble chart generation failed: {e}")

    return generated_files


def generate_volcano_plot(
    volcano_csv: Path,
    outdir: Path,
    config: Any = None,
    fdr_cutoff: float = 0.05,
    fold_cutoff: float = 1.5
) -> Dict[str, Path]:
    """
    Generates a publication-grade Volcano Plot plotting -log10(FDR p-value) vs log2(Fold Enrichment).
    """
    setup_matplotlib_style()
    volcano_csv = Path(volcano_csv)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    generated_files = {}
    
    if not volcano_csv.exists() or os.path.getsize(volcano_csv) <= 1:
        logger.warning(f"Volcano CSV file {volcano_csv} is empty or missing. Skipping volcano plot.")
        return generated_files

    try:
        df = pd.read_csv(volcano_csv)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as e:
        logger.warning(f"Volcano CSV file {volcano_csv} could not be parsed: {e}. Skipping volcano plot.")
        return generated_files

    if df.empty:
        logger.warning("Empty volcano data. Skipping volcano plot.")
        return generated_files
        
    df = df.dropna(subset=["raw_pvalue"]).copy()
    if df.empty:
        logger.warning("No valid p-values in volcano CSV. Skipping.")
        return generated_files

    # Log coordinates
    df["log2_fold_enrichment"] = df["fold_enrichment"].apply(lambda x: np.log2(x) if x > 0 else -10.0)
    df["minus_log10_fdr_pvalue"] = df["adjusted_pvalue_fdr"].apply(lambda x: -np.log10(x) if x > 0 else 300.0)

    # Re-classify significance based on thresholds
    df["sig_type"] = "NOT_SIGNIFICANT"
    df.loc[
        (df["adjusted_pvalue_fdr"] <= fdr_cutoff) & 
        (df["fold_enrichment"] >= fold_cutoff),
        "sig_type"
    ] = "SIGNIFICANT"

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Split by significance
    df_sig = df[df["sig_type"] == "SIGNIFICANT"]
    df_nonsig = df[df["sig_type"] == "NOT_SIGNIFICANT"]

    # Draw scatters
    ax.scatter(
        df_nonsig["log2_fold_enrichment"],
        df_nonsig["minus_log10_fdr_pvalue"],
        c=OKABE_ITO["muted_grey"],
        alpha=0.6,
        edgecolors="black",
        linewidths=0.5,
        s=40,
        label="Not Significant"
    )
    ax.scatter(
        df_sig["log2_fold_enrichment"],
        df_sig["minus_log10_fdr_pvalue"],
        c=OKABE_ITO["vermilion"],
        alpha=0.9,
        edgecolors="black",
        linewidths=0.6,
        s=65,
        label=f"Significant Enriched (FDR <= {fdr_cutoff}, FE >= {fold_cutoff})"
    )

    # Cutoff threshold lines
    ax.axhline(-np.log10(fdr_cutoff), color="black", linestyle=":", linewidth=1.0)
    ax.axvline(np.log2(fold_cutoff), color="black", linestyle=":", linewidth=1.0)

    # Label top 8 significant pathways
    df_top_labels = df_sig.sort_values("raw_pvalue").head(8)
    for idx, row in df_top_labels.iterrows():
        ax.annotate(
            row["description"],
            xy=(row["log2_fold_enrichment"], row["minus_log10_fdr_pvalue"]),
            xytext=(6, 5),
            textcoords="offset points",
            fontsize=8,
            fontweight="semibold",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="black", alpha=0.7, lw=0.4)
        )

    ax.set_xlabel("log2(Fold Enrichment)")
    ax.set_ylabel("-log10(Adjusted p-value, FDR)")
    ax.set_title("Volcano Distribution of Pathway Enrichment", pad=15)
    ax.grid(True)
    ax.legend(loc="upper left")
    
    add_reproducibility_metadata(fig, config)
    volcano_path = outdir / "pathway_volcano_plot"
    save_publication_figure(fig, volcano_path)
    plt.close()
    
    generated_files["volcano"] = volcano_path.with_suffix(".png")
    logger.info(f"Saved volcano plot in multiple formats at {volcano_path}.*")

    # Plotly interactive equivalent
    if HAS_PLOTLY:
        try:
            hover_cols = ["pathway_id", "fold_enrichment", "adjusted_pvalue_fdr"]
            if "query_count_k" in df.columns:
                hover_cols.append("query_count_k")

            fig_int = px.scatter(
                df,
                x="log2_fold_enrichment",
                y="minus_log10_fdr_pvalue",
                color="sig_type",
                color_discrete_map={"NOT_SIGNIFICANT": "grey", "SIGNIFICANT": "red"},
                hover_name="description",
                hover_data=hover_cols,
                labels={"log2_fold_enrichment": "log2(Fold Enrichment)", "minus_log10_fdr_pvalue": "-log10 FDR", "sig_type": "Significance"},
                title="Volcano Distribution of Pathway Enrichment (Interactive)"
            )
            fig_int.add_vline(x=np.log2(fold_cutoff), line_dash="dot", line_color="black")
            fig_int.add_hline(y=-np.log10(fdr_cutoff), line_dash="dot", line_color="black")
            fig_int.update_layout(template="plotly_white")
            
            volc_html = outdir / "pathway_volcano_plot_interactive.html"
            pio.write_html(fig_int, str(volc_html))
            generated_files["volcano_html"] = volc_html
        except Exception as e:
            logger.warning(f"Plotly volcano generation failed: {e}")

    return generated_files


def generate_enrichment_heatmap(
    enrichment_matrix_csv: Path,
    outdir: Path,
    config: Any = None
) -> Optional[Path]:
    """
    Plots a scientific cluster heatmap mapping unique viral accessions to host pathways.
    """
    setup_matplotlib_style()
    enrichment_matrix_csv = Path(enrichment_matrix_csv)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    if not enrichment_matrix_csv.exists() or os.path.getsize(enrichment_matrix_csv) <= 1:
        logger.warning(f"Enrichment matrix {enrichment_matrix_csv} missing or empty. Skipping heatmap.")
        return None

    try:
        df = pd.read_csv(enrichment_matrix_csv, index_col=0)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as e:
        logger.warning(f"Enrichment matrix {enrichment_matrix_csv} could not be parsed: {e}. Skipping heatmap.")
        return None

    if df.empty or df.shape[0] < 1 or df.shape[1] < 1:
        logger.warning("Empty matrix cells. Skipping heatmap.")
        return None

    # Cap visual limits to prevent visual clutter
    if df.shape[0] > 40:
        df = df.head(40) # top 40 genes
    if df.shape[1] > 30:
        df = df.iloc[:, :30] # top 30 pathways

    # Draw Heatmap
    fig, ax = plt.subplots(figsize=(max(8, df.shape[1] * 0.35 + 2), max(6, df.shape[0] * 0.25 + 2)))
    
    sns.heatmap(
        df,
        cmap="cividis",
        linewidths=0.5,
        linecolor="#dddddd",
        cbar_kws={"label": "Link Active (1=Yes, 0=No)", "shrink": 0.7},
        ax=ax
    )
    
    ax.set_title("Gene-Pathway Representation Profile Heatmap\n(Binary co-occurrence coordinates mapping query genes to pathways)", pad=15)
    ax.set_ylabel("Viral Gene Accessions (collapsed)")
    ax.set_xlabel("Biological Pathways")
    
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=8)
    
    add_reproducibility_metadata(fig, config)
    
    heatmap_path = outdir / "gene_pathway_heatmap"
    save_publication_figure(fig, heatmap_path)
    plt.close()
    
    logger.info(f"Saved gene-pathway profile cluster heatmap in multiple formats at {heatmap_path}.*")
    return heatmap_path.with_suffix(".png")


def generate_pca_plot(
    enrichment_matrix_csv: Path,
    outdir: Path,
    config: Any = None
) -> Optional[Path]:
    """
    Performs custom SVD Principal Component Analysis (PCA) on gene-pathway functional profiles
    and plots PCA components PC1/PC2 along with explained variance.
    """
    setup_matplotlib_style()
    enrichment_matrix_csv = Path(enrichment_matrix_csv)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    if not enrichment_matrix_csv.exists() or os.path.getsize(enrichment_matrix_csv) <= 1:
        logger.warning(f"Enrichment matrix {enrichment_matrix_csv} missing or empty. Skipping PCA plot.")
        return None

    try:
        df = pd.read_csv(enrichment_matrix_csv, index_col=0)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as e:
        logger.warning(f"Enrichment matrix {enrichment_matrix_csv} could not be parsed: {e}. Skipping PCA plot.")
        return None

    if df.empty or df.shape[0] < 3 or df.shape[1] < 2:
        logger.info("Enrichment matrix is too small for PCA dimensional projection (needs >=3 genes, >=2 pathways). Skipping.")
        return None

    # Perform self-contained numpy SVD PCA
    X = df.values.astype(float)
    # Center the matrix
    X_centered = X - X.mean(axis=0)
    
    try:
        # Singular Value Decomposition
        U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
        
        # Project data
        PC = X_centered @ Vt.T
        
        # Calculate variance explained
        var_explained = (S ** 2) / np.sum(S ** 2)
        pc1_var = var_explained[0] * 100
        pc2_var = var_explained[1] * 100 if len(var_explained) > 1 else 0.0
        
        # Make plot
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Scatter projected coordinates
        ax.scatter(
            PC[:, 0],
            PC[:, 1] if PC.shape[1] > 1 else np.zeros_like(PC[:, 0]),
            c=OKABE_ITO["blue"],
            alpha=0.85,
            edgecolors="black",
            linewidths=0.6,
            s=80,
            zorder=3
        )
        
        # Draw centroid cross
        ax.axhline(0, color="grey", linestyle="--", linewidth=0.8, zorder=1)
        ax.axvline(0, color="grey", linestyle="--", linewidth=0.8, zorder=1)
        
        # Annotate points with gene IDs
        for i, gene_id in enumerate(df.index):
            ax.annotate(
                gene_id,
                xy=(PC[i, 0], PC[i, 1] if PC.shape[1] > 1 else 0.0),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="grey", alpha=0.6, lw=0.3)
            )
            
        ax.set_xlabel(f"PC1 ({pc1_var:.2f}% Explained Variance)")
        ax.set_ylabel(f"PC2 ({pc2_var:.2f}% Explained Variance)" if PC.shape[1] > 1 else "PC2 (0.00%)")
        ax.set_title("PCA Dimensional Projection of Gene Functional Profiles\n(PCA of pathway representation binary matrix)", pad=15)
        ax.grid(True)
        
        add_reproducibility_metadata(fig, config)
        
        pca_path = outdir / "functional_pca_plot"
        save_publication_figure(fig, pca_path)
        plt.close()
        
        logger.info(f"Saved SVD functional profile PCA projection in multiple formats at {pca_path}.*")
        
        # Interactive Plotly PCA plot
        if HAS_PLOTLY:
            try:
                df_pca = pd.DataFrame({
                    "PC1": PC[:, 0],
                    "PC2": PC[:, 1] if PC.shape[1] > 1 else np.zeros_like(PC[:, 0]),
                    "Gene": df.index
                })
                fig_int = px.scatter(
                    df_pca,
                    x="PC1",
                    y="PC2",
                    text="Gene",
                    labels={"PC1": f"PC1 ({pc1_var:.1f}%)", "PC2": f"PC2 ({pc2_var:.1f}%)"},
                    title="PCA Dimensional Projection of Gene Profiles (Interactive)"
                )
                fig_int.update_traces(marker=dict(size=12, line=dict(width=1, color="DarkSlateGrey")), textposition='top center')
                fig_int.update_layout(template="plotly_white")
                
                pca_html = outdir / "functional_pca_plot_interactive.html"
                pio.write_html(fig_int, str(pca_html))
            except Exception as e:
                logger.warning(f"Plotly PCA projection generation failed: {e}")
                
        return pca_path.with_suffix(".png")
        
    except Exception as e:
        logger.error(f"SVD projection or calculation failed in PCA visualizer: {e}")
        return None


def generate_protein_pathway_network(
    mapped_pathways_csv: Path,
    enrichment_csv: Path,
    outdir: Path,
    max_pathways: int = 5
) -> Optional[Path]:
    """
    Creates an bipartite protein-pathway interaction network graph.
    Connects unique, collapsed viral UniProt IDs to their associated biological pathways.
    """
    setup_matplotlib_style()
    mapped_pathways_csv = Path(mapped_pathways_csv)
    enrichment_csv = Path(enrichment_csv)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    if not mapped_pathways_csv.exists() or not enrichment_csv.exists():
        logger.warning("Mapped pathways or enrichment csv missing. Skipping network visualization.")
        return None
        
    df_mappings = pd.read_csv(mapped_pathways_csv)
    df_enrich = pd.read_csv(enrichment_csv)
    
    if df_mappings.empty or df_enrich.empty:
        logger.warning("Empty mapping or enrichment data. Network skipped.")
        return None
        
    df_mappings = df_mappings.drop_duplicates(subset=["uniprot_id", "pathway_id"]).copy()
    
    df_sig_pathways = df_enrich.sort_values("raw_pvalue")
    sig_paths = df_sig_pathways[df_sig_pathways["adjusted_pvalue_fdr"] <= 0.05]["pathway_id"].tolist()
    if not sig_paths:
        sig_paths = df_sig_pathways.head(max_pathways)["pathway_id"].tolist()
    else:
        sig_paths = sig_paths[:max_pathways]
        
    df_filtered = df_mappings[df_mappings["pathway_id"].isin(sig_paths)].copy()
    if df_filtered.empty:
        logger.warning("No mappings found matching significant pathways. Skipping network.")
        return None
        
    G = nx.Graph()
    
    proteins = list(df_filtered["uniprot_id"].unique())
    pathways = list(df_filtered["pathway_id"].unique())
    
    pathway_names = {}
    for idx, row in df_enrich.iterrows():
        pathway_names[row["pathway_id"]] = row["description"]
        
    G.add_nodes_from(proteins, bipartite=0)
    pathway_labels = {pid: pathway_names.get(pid, pid) for pid in pathways}
    G.add_nodes_from(pathways, bipartite=1)
    
    for idx, row in df_filtered.iterrows():
        G.add_edge(row["uniprot_id"], row["pathway_id"])
        
    # Stunning dark slate theme matching app portal aesthetics
    fig, ax = plt.subplots(figsize=(10, 8), facecolor="#0b0f19")
    ax.set_facecolor("#0b0f19")
    
    # Advanced layout with balanced node distribution forces
    pos = nx.spring_layout(G, k=0.45, iterations=100, seed=42)
    
    # 1. Draw elegant edges (semi-transparent indigo/slate)
    nx.draw_networkx_edges(
        G, pos,
        width=1.5,
        edge_color="#475569",
        alpha=0.45,
        ax=ax
    )
    
    # 2. Draw Translucent Glowing Halos (Simulates light glow around nodes)
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=proteins,
        node_shape="o",
        node_color="#38bdf8",
        node_size=600,
        alpha=0.18,
        ax=ax
    )
    
    # Degree-weighted sizes for pathways (higher degrees = larger diamonds)
    pathway_sizes = [550 + G.degree(node) * 120 for node in pathways]
    pathway_halos = [s * 1.5 for s in pathway_sizes]
    
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=pathways,
        node_shape="d",
        node_color="#c084fc",
        node_size=pathway_halos,
        alpha=0.15,
        ax=ax
    )
    
    # 3. Draw Core High-Fidelity Nodes
    # Viral proteins: Neon Cyan circles
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=proteins,
        node_shape="o",
        node_color="#38bdf8",
        node_size=320,
        edgecolors="#0f172a",
        linewidths=1.2,
        alpha=0.95,
        ax=ax
    )
    
    # Host pathways: Neon Violet diamonds
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=pathways,
        node_shape="d",
        node_color="#c084fc",
        node_size=pathway_sizes,
        edgecolors="#0f172a",
        linewidths=1.2,
        alpha=0.95,
        ax=ax
    )
    
    # 4. Render clean white/slate labels above nodes
    labels = {}
    for node in G.nodes():
        if node in proteins:
            labels[node] = node
        else:
            labels[node] = pathway_labels.get(node, node)
            
    pos_labels = {k: [v[0], v[1] + 0.055] for k, v in pos.items()}
    
    nx.draw_networkx_labels(
        G, pos_labels,
        labels=labels,
        font_size=7.5,
        font_weight="bold",
        font_color="#e2e8f0",
        ax=ax
    )
    
    # Curated modern legends
    patch_prot = patches.Patch(color="#38bdf8", label="Viral Protein Homology Matches")
    patch_path = patches.Patch(color="#c084fc", label="Implicated Host Pathways (Degree Hijack)")
    
    leg = plt.legend(
        handles=[patch_prot, patch_path],
        loc="upper right",
        frameon=True,
        facecolor="#0f172a",
        edgecolor="#1e293b"
    )
    plt.setp(leg.get_texts(), color="#cbd5e1", fontsize=9)
    
    plt.title(
        "Protein-Pathway Bipartite Interaction Network Map\n(Links between predicted viral proteins and host molecular pathways)",
        color="#f1f5f9",
        fontsize=12,
        fontweight="bold",
        pad=20
    )
    plt.axis("off")
    
    network_path = outdir / "protein_pathway_network"
    save_publication_figure(fig, network_path)
    plt.close()
    
    logger.info(f"Saved bipartite protein-pathway network map to {network_path}.*")
    return network_path.with_suffix(".png")


def generate_annotation_distribution_plots(
    annotated_proteins_csv: Path,
    outdir: Path,
    config: Any = None
) -> Optional[Path]:
    """
    Analyzes sequence similarity annotation statistics and generates
    a multi-panel grid plotting metrics in multiple formats.
    """
    setup_matplotlib_style()
    annotated_proteins_csv = Path(annotated_proteins_csv)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    if not annotated_proteins_csv.exists():
        logger.warning(f"Annotated proteins CSV {annotated_proteins_csv} is missing. Skipping distributions.")
        return None
        
    df = pd.read_csv(annotated_proteins_csv)
    if df.empty:
        logger.warning("Empty annotated proteins. Skipping distributions.")
        return None
        
    df_hits = df[
        (df["uniprot_id"].notna()) & 
        (df["uniprot_id"].astype(str).str.upper() != "NONE") &
        (df["uniprot_id"].astype(str).str.upper() != "HYPOTHETICAL PROTEIN")
    ].copy()
    
    if df_hits.empty:
        logger.info("Zero annotated proteins mapping with alignments. Skipping sequence alignment distributions.")
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No sequence similarity alignments\nexceeding confidence thresholds.",
                horizontalalignment="center", verticalalignment="center", fontsize=12, color="grey")
        ax.axis("off")
        info_path = outdir / "annotation_distributions"
        save_publication_figure(fig, info_path)
        plt.close()
        return info_path.with_suffix(".png")

    df_hits["identity_pct"] = pd.to_numeric(df_hits["identity_pct"], errors="coerce")
    df_hits["query_coverage_pct"] = pd.to_numeric(df_hits["query_coverage_pct"], errors="coerce")
    df_hits["bitscore"] = pd.to_numeric(df_hits["bitscore"], errors="coerce")
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    
    # 1. Percent Identity Histogram
    ax1 = axes[0]
    sns.histplot(df_hits["identity_pct"].dropna(), bins=15, kde=True, ax=ax1, color=OKABE_ITO["blue"], edgecolor="black")
    ax1.set_title("Sequence Identity Distribution")
    ax1.set_xlabel("Percent Identity (%)")
    ax1.set_ylabel("Aligned ORFs Count")
    ax1.grid(True)
    
    # 2. Query Coverage Histogram
    ax2 = axes[1]
    sns.histplot(df_hits["query_coverage_pct"].dropna(), bins=15, kde=True, ax=ax2, color=OKABE_ITO["bluish_green"], edgecolor="black")
    ax2.set_title("Query Coverage Distribution")
    ax2.set_xlabel("Query Coverage (%)")
    ax2.set_ylabel("Aligned ORFs Count")
    ax2.grid(True)
    
    # 3. Bitscore Distribution
    ax3 = axes[2]
    sns.histplot(df_hits["bitscore"].dropna(), bins=15, kde=True, ax=ax3, color=OKABE_ITO["orange"], edgecolor="black")
    ax3.set_title("Alignment Bitscore Distribution")
    ax3.set_xlabel("Bitscore")
    ax3.set_ylabel("Aligned ORFs Count")
    ax3.grid(True)
    
    plt.suptitle("Protein Functional Annotation Similarity Metrics\n(QC filtered Swiss-Prot database alignments)", y=1.02)
    plt.tight_layout()
    
    add_reproducibility_metadata(fig, config)
    
    dist_path = outdir / "annotation_distributions"
    save_publication_figure(fig, dist_path)
    plt.close()
    
    logger.info(f"Saved functional annotation similarity metric distributions to {dist_path}.*")
    return dist_path.with_suffix(".png")


def generate_orf_distribution_plots(
    orf_coordinates_gff: Path,
    outdir: Path,
    config: Any = None
) -> Dict[str, Path]:
    """
    Parses predicted ORF coordinates and lengths from GFF3 and generates
    length distributions and coordinates tracks mapping to PNG, SVG, PDF.
    """
    setup_matplotlib_style()
    orf_coordinates_gff = Path(orf_coordinates_gff)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    generated_files = {}
    
    if not orf_coordinates_gff.exists():
        logger.warning(f"GFF3 coordinates file {orf_coordinates_gff} is missing. Skipping ORF distributions.")
        return generated_files
        
    records = []
    genome_size = 0
    
    with open(orf_coordinates_gff, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            
            attrs = {}
            for item in parts[8].split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attrs[k.strip()] = v.strip()
            
            frame = int(attrs.get("frame", 0))
            orf_id = attrs.get("ID", f"ORF_{start}_{end}")
            
            length_bp = end - start + 1
            length_aa = length_bp // 3
            
            records.append({
                "orf_id": orf_id,
                "start": start,
                "end": end,
                "strand": strand,
                "frame": frame,
                "length_bp": length_bp,
                "length_aa": length_aa
            })
            
            genome_size = max(genome_size, end)
            
    if not records:
        logger.warning("No ORF records parsed from GFF3. Skipping plots.")
        return generated_files
        
    df = pd.DataFrame(records)
    
    # --- Plot 1: ORF Length Distribution Histogram ---
    fig, ax = plt.subplots(figsize=(7, 4.5))
    sns.histplot(df["length_aa"], bins=15, kde=True, color=OKABE_ITO["sky_blue"], edgecolor="black", ax=ax)
    ax.set_title("Predicted Open Reading Frames (ORFs) Length Distribution")
    ax.set_xlabel("ORF Length (amino acids)")
    ax.set_ylabel("ORF Count")
    ax.grid(True)
    
    add_reproducibility_metadata(fig, config)
    len_path = outdir / "orf_lengths"
    save_publication_figure(fig, len_path)
    plt.close()
    generated_files["lengths"] = len_path.with_suffix(".png")
    logger.info(f"Saved ORF length distribution chart to {len_path}.*")

    # --- Plot 2: Genomic Track Map (Okabe-Ito Strand Coloring) ---
    fig, ax = plt.subplots(figsize=(12, 5))
    
    ax.axhline(0, color="black", linewidth=2.0, zorder=1)
    
    ticks_step = max(100, int(genome_size / 100) * 10)
    for tick in range(0, genome_size + 1, ticks_step):
        ax.axvline(tick, color="#cccccc", linestyle=":", linewidth=0.8, zorder=0)
        
    y_offset = {
        "+": 0.4,
        "-": -0.8
    }
    
    for idx, row in df.iterrows():
        start = row["start"]
        end = row["end"]
        strand = row["strand"]
        frame = row["frame"]
        length = end - start + 1
        
        color = STRAND_COLORS.get(strand, "#888888")
        
        lane = frame % 3
        y_pos = y_offset[strand] + (lane * 0.12 if strand == "+" else -lane * 0.12)
        
        direction = 1 if strand == "+" else -1
        arrow_head = min(50, length * 0.2)
        
        rect = patches.FancyArrow(
            start if direction == 1 else end,
            y_pos,
            direction * (length - arrow_head),
            0,
            width=0.08,
            head_width=0.15,
            head_length=arrow_head,
            length_includes_head=True,
            facecolor=color,
            edgecolor="black",
            linewidth=0.5,
            alpha=0.85,
            zorder=2
        )
        ax.add_patch(rect)
        
        if length > genome_size * 0.05:
            ax.text(
                start + length / 2,
                y_pos + 0.11 if strand == "+" else y_pos - 0.15,
                row["orf_id"].split("_")[-1],
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=7.5,
                color="#333333",
                zorder=3
            )
            
    ax.set_ylim(-1.5, 1.0)
    ax.set_xlim(-100, genome_size + 100)
    ax.set_xlabel("Genomic Coordinate (bp)", fontsize=11)
    
    ax.set_yticks([0.5, -0.9])
    ax.set_yticklabels(["Positive Strand (+)", "Negative Strand (-)"], fontweight="bold")
    
    for spine in ["top", "left", "right"]:
        ax.spines[spine].set_visible(False)
        
    patch_pos = patches.Patch(color=STRAND_COLORS["+"], label="Forward Strand ORFs (+)")
    patch_neg = patches.Patch(color=STRAND_COLORS["-"], label="Reverse Strand ORFs (-)")
    ax.legend(handles=[patch_pos, patch_neg], loc="upper right")
    
    plt.title(f"PathoScope AI Predicted Viral Genome Coordinate Track Map\n(Total predicted genome sequence length = {genome_size:,} bp)", pad=15)
    
    add_reproducibility_metadata(fig, config)
    track_path = outdir / "orf_genomic_track"
    save_publication_figure(fig, track_path)
    plt.close()
    
    generated_files["track"] = track_path.with_suffix(".png")
    logger.info(f"Saved ORF genomic coordinate tracks to {track_path}.*")
    
    return generated_files


def generate_preprocessing_qc_plots(
    qc_data: Dict[str, Any],
    outdir: Path,
    config: Any = None
) -> Dict[str, Path]:
    """
    Generates high-fidelity, journal-grade visualizations for NGS preprocessing:
    - Per-Base Sequence Quality Distribution
    - Per-Read GC Content Distribution
    - Sequence Duplication Level histogram
    - Read Retention statistics waterfall
    Also generates interactive Plotly dashboards if available.
    """
    setup_matplotlib_style()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    generated_files = {}
    
    qc_plots = qc_data.get("qc_plots", {})
    if not qc_plots:
        logger.warning("No qc_plots dictionary found in preprocessing report. Skipping preprocessing QC plots.")
        return generated_files

    # 1. Per-Base Quality Score Plot
    per_base_quality = qc_plots.get("per_base_quality", [])
    if per_base_quality:
        try:
            fig, ax = plt.subplots(figsize=(8, 4.5))
            positions = list(range(1, len(per_base_quality) + 1))
            ax.plot(positions, per_base_quality, color=OKABE_ITO["blue"], linewidth=2.0, label="Mean Q-score")
            
            # Draw Q-score threshold bands
            ax.axhspan(28, 41, color="#34d399", alpha=0.15, zorder=0) # Excellent green
            ax.axhspan(20, 28, color="#f59e0b", alpha=0.12, zorder=0) # Moderate orange
            ax.axhspan(0, 20, color="#f87171", alpha=0.12, zorder=0)  # Poor red
            
            ax.set_ylim(0, 41)
            ax.set_xlabel("Read Base Position (bp)", fontsize=11)
            ax.set_ylabel("Quality Score (Phred Q)", fontsize=11)
            ax.set_title("Per-Base Sequence Quality Distribution\n(Line shows mean Phred Q scores across query reads)", pad=15)
            ax.grid(True)
            ax.legend(loc="lower left")
            
            add_reproducibility_metadata(fig, config)
            q_path = outdir / "per_base_quality"
            save_publication_figure(fig, q_path)
            plt.close()
            
            generated_files["per_base_quality"] = q_path.with_suffix(".png")
            logger.info(f"Saved per-base quality distribution chart to {q_path}.*")
            
            if HAS_PLOTLY:
                fig_int = go.Figure()
                fig_int.add_trace(go.Scatter(
                    x=positions, y=per_base_quality,
                    mode='lines',
                    line=dict(color=OKABE_ITO["blue"], width=2.5),
                    name='Mean Q-score'
                ))
                # Add background color bands
                fig_int.add_hrect(y0=28, y1=41, fillcolor="rgba(52, 211, 153, 0.15)", line_width=0)
                fig_int.add_hrect(y0=20, y1=28, fillcolor="rgba(245, 158, 11, 0.12)", line_width=0)
                fig_int.add_hrect(y0=0, y1=20, fillcolor="rgba(248, 113, 113, 0.12)", line_width=0)
                fig_int.update_layout(
                    title="Per-Base Sequence Quality (Interactive)",
                    xaxis_title="Read Base Position (bp)",
                    yaxis_title="Quality Score (Phred Q)",
                    yaxis_range=[0, 41],
                    template="plotly_white"
                )
                q_html = outdir / "per_base_quality_interactive.html"
                pio.write_html(fig_int, str(q_html))
        except Exception as e:
            logger.warning(f"Failed to generate per-base quality plots: {e}")

    # 2. Per-Read GC Content Distribution
    gc_dist = qc_plots.get("gc_content_distribution", {})
    if gc_dist:
        try:
            x_vals = [int(k) for k in gc_dist.keys()]
            y_vals = [int(v) for v in gc_dist.values()]
            xy = sorted(zip(x_vals, y_vals))
            x_sorted, y_sorted = zip(*xy)
            
            fig, ax = plt.subplots(figsize=(8, 4.5))
            ax.plot(x_sorted, y_sorted, color=OKABE_ITO["bluish_green"], linewidth=2.0, label="GC density")
            ax.fill_between(x_sorted, y_sorted, color=OKABE_ITO["bluish_green"], alpha=0.18)
            ax.set_xlabel("GC Content (%)", fontsize=11)
            ax.set_ylabel("Reads Count", fontsize=11)
            ax.set_title("Per-Read GC Content Distribution\n(Plot shows density distribution of nucleotide G-C ratios)", pad=15)
            ax.grid(True)
            
            add_reproducibility_metadata(fig, config)
            gc_path = outdir / "gc_distribution"
            save_publication_figure(fig, gc_path)
            plt.close()
            
            generated_files["gc_distribution"] = gc_path.with_suffix(".png")
            logger.info(f"Saved GC content distribution chart to {gc_path}.*")
            
            if HAS_PLOTLY:
                fig_int = px.line(
                    x=x_sorted, y=y_sorted,
                    labels={"x": "GC Content (%)", "y": "Reads Count"},
                    title="Per-Read GC Content Distribution (Interactive)"
                )
                fig_int.update_traces(line_color=OKABE_ITO["bluish_green"], fill='tozeroy')
                fig_int.update_layout(template="plotly_white")
                gc_html = outdir / "gc_distribution_interactive.html"
                pio.write_html(fig_int, str(gc_html))
        except Exception as e:
            logger.warning(f"Failed to generate GC distribution plots: {e}")

    # 3. Sequence Duplication Levels
    dup_levels = qc_plots.get("sequence_duplication_levels", {})
    if dup_levels:
        try:
            labels = list(dup_levels.keys())
            counts = [int(v) for v in dup_levels.values()]
            
            fig, ax = plt.subplots(figsize=(8, 4.5))
            ax.bar(labels, counts, color=OKABE_ITO["orange"], edgecolor="black", linewidth=0.6, width=0.45)
            ax.set_xlabel("Sequence Duplication Level", fontsize=11)
            ax.set_ylabel("Reads Count", fontsize=11)
            ax.set_title("Sequence Duplication Levels\n(Bar shows distribution of read content duplication)", pad=15)
            ax.grid(True, axis="y")
            
            add_reproducibility_metadata(fig, config)
            dup_path = outdir / "sequence_duplication"
            save_publication_figure(fig, dup_path)
            plt.close()
            
            generated_files["sequence_duplication"] = dup_path.with_suffix(".png")
            logger.info(f"Saved sequence duplication levels chart to {dup_path}.*")
            
            if HAS_PLOTLY:
                fig_int = px.bar(
                    x=labels, y=counts,
                    labels={"x": "Sequence Duplication Level", "y": "Reads Count"},
                    title="Sequence Duplication Levels (Interactive)"
                )
                fig_int.update_traces(marker_color=OKABE_ITO["orange"])
                fig_int.update_layout(template="plotly_white")
                dup_html = outdir / "sequence_duplication_interactive.html"
                pio.write_html(fig_int, str(dup_html))
        except Exception as e:
            logger.warning(f"Failed to generate duplication level plots: {e}")

    # 4. Read Retention Waterfall
    waterfall = qc_plots.get("read_retention_waterfall", {})
    if waterfall:
        try:
            categories = ["Raw Reads", "Clipped Adapters", "Short Length", "Low Quality", "Passed QC"]
            values = [
                int(waterfall.get("raw_reads", 0)),
                int(waterfall.get("adapter_trimmed", 0)),
                int(waterfall.get("length_filtered", 0)),
                int(waterfall.get("quality_filtered", 0)),
                int(waterfall.get("kept_reads", 0))
            ]
            
            colors = [
                OKABE_ITO["black"],
                OKABE_ITO["sky_blue"],
                OKABE_ITO["yellow"],
                OKABE_ITO["vermilion"],
                OKABE_ITO["bluish_green"]
            ]
            
            fig, ax = plt.subplots(figsize=(8, 4.5))
            ax.bar(categories, values, color=colors, edgecolor="black", linewidth=0.6, width=0.45)
            ax.set_ylabel("Reads Count", fontsize=11)
            ax.set_title("Sequence Preprocessing & Read Retention Waterfall\n(Waterfall tracks sequence filtration milestones)", pad=15)
            ax.grid(True, axis="y")
            
            add_reproducibility_metadata(fig, config)
            wat_path = outdir / "read_retention_waterfall"
            save_publication_figure(fig, wat_path)
            plt.close()
            
            generated_files["read_retention_waterfall"] = wat_path.with_suffix(".png")
            logger.info(f"Saved read retention waterfall chart to {wat_path}.*")
            
            if HAS_PLOTLY:
                fig_int = px.bar(
                    x=categories, y=values,
                    color=categories,
                    color_discrete_sequence=colors,
                    labels={"x": "Filter Milestones", "y": "Reads Count"},
                    title="Read Retention Statistics Waterfall (Interactive)"
                )
                fig_int.update_layout(template="plotly_white", showlegend=False)
                wat_html = outdir / "read_retention_waterfall_interactive.html"
                pio.write_html(fig_int, str(wat_html))
        except Exception as e:
            logger.warning(f"Failed to generate read retention plots: {e}")
            
    return generated_files


def generate_sspa_plots(
    ranking_csv: Path,
    outdir: Path,
    config: Any = None,
    top_n: int = 15
) -> Dict[str, Path]:
    """
    Generates publication-grade visualizations for Single-Sample Pathway Analysis (ssPA/ssGSEA)
    and Multi-Evidence Integrated Pathway Scores.
    """
    setup_matplotlib_style()
    ranking_csv = Path(ranking_csv)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    generated_files = {}
    
    if not ranking_csv.exists() or os.path.getsize(ranking_csv) <= 1:
        logger.warning(f"Ranking report CSV file {ranking_csv} is empty or missing. Skipping ssPA plots.")
        return generated_files

    try:
        df = pd.read_csv(ranking_csv)
    except Exception as e:
        logger.warning(f"Ranking report CSV file {ranking_csv} could not be parsed: {e}. Skipping ssPA plots.")
        return generated_files

    if df.empty:
        logger.warning("Empty ranking data. Skipping ssPA visualization.")
        return generated_files

    # Check if ssGSEA columns exist
    if "ssgsea_enrichment_score_normalized" not in df.columns or "multi_evidence_pathway_score" not in df.columns:
        logger.warning("ssGSEA or Multi-Evidence Score columns missing from ranking report. Skipping ssPA plots.")
        return generated_files

    df = df.dropna(subset=["ssgsea_enrichment_score_normalized"]).copy()
    if df.empty:
        return generated_files

    df_sorted = df.sort_values("multi_evidence_pathway_score", ascending=False).head(top_n)
    df_plot = df_sorted.iloc[::-1]  # Horizontal reverse order

    # --- Plot 1: ssGSEA Pathway Enrichment Scores Barplot ---
    fig, ax = plt.subplots(figsize=(9, min(6, len(df_plot) * 0.4 + 2)))
    
    norm = plt.Normalize(0.0, 10.0)
    cmap = sns.color_palette("viridis", as_cmap=True)
    colors = [cmap(norm(val)) for val in df_plot["multi_evidence_pathway_score"]]
    
    bars = ax.barh(
        df_plot["description"],
        df_plot["ssgsea_enrichment_score_normalized"],
        color=colors,
        edgecolor=OKABE_ITO["black"],
        linewidth=0.6,
        height=0.6
    )
    
    ax.set_xlabel("Normalized ssGSEA Enrichment Score")
    ax.set_ylabel("Biological Pathways")
    ax.set_title("Single-Sample Pathway Analysis (ssPA)\n(Top Pathways by ssGSEA, colored by Multi-Evidence Score)", pad=15)
    ax.grid(True, axis="x")
    
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", pad=0.02, shrink=0.7)
    cbar.set_label("Multi-Evidence Score (out of 10.0)", fontsize=9)
    
    add_reproducibility_metadata(fig, config)
    
    barplot_path = outdir / "sspa_ssgsea_barplot"
    save_publication_figure(fig, barplot_path)
    plt.close()
    
    generated_files["sspa_barplot"] = barplot_path.with_suffix(".png")
    logger.info(f"Saved horizontal ssPA barplot at {barplot_path}.*")

    # --- Plot 2: Multi-Evidence Pathway Score Bubble Plot ---
    fig, ax = plt.subplots(figsize=(10, min(7, len(df_plot) * 0.5 + 2)))
    
    k_counts = df_plot["query_count_k"].fillna(2.0)
    k_min, k_max = k_counts.min(), k_counts.max()
    sizes = k_counts.apply(lambda x: 100 + (x - k_min) / max(1, k_max - k_min) * 500 if k_max > k_min else 250)
    
    nes_vals = df_plot["normalized_enrichment_score_nes"].fillna(0.0)
    sc = ax.scatter(
        df_plot["multi_evidence_pathway_score"],
        df_plot["description"],
        s=sizes,
        c=nes_vals,
        cmap="coolwarm",
        alpha=0.85,
        edgecolors="black",
        linewidths=0.7
    )
    
    ax.set_xlabel("Multi-Evidence Integrated Pathway Score (0.0 to 10.0)")
    ax.set_ylabel("Biological Pathways")
    ax.set_title("Multi-Evidence Integrated Pathway Priority\n(Dot size = Query Gene Count; Color = GSEA NES)", pad=15)
    ax.grid(True, axis="both")
    
    cbar = fig.colorbar(sc, ax=ax, pad=0.02, shrink=0.7)
    cbar.set_label("GSEA Normalized Enrichment Score (NES)", fontsize=9)
    
    handles = []
    vals = np.linspace(k_min, k_max, min(3, int(k_max - k_min + 1))).astype(int)
    for v in vals:
        sz = 100 + (v - k_min) / max(1, k_max - k_min) * 500 if k_max > k_min else 250
        handles.append(
            ax.scatter([], [], s=sz, c="grey", alpha=0.6, edgecolors="black", label=f"{v} query genes")
        )
    ax.legend(handles=handles, loc="lower right", title="Pathway Counts (k)")
    
    add_reproducibility_metadata(fig, config)
    bubbleplot_path = outdir / "sspa_multi_evidence_bubbleplot"
    save_publication_figure(fig, bubbleplot_path)
    plt.close()
    
    generated_files["sspa_bubbleplot"] = bubbleplot_path.with_suffix(".png")
    logger.info(f"Saved multi-evidence bubbleplot at {bubbleplot_path}.*")

    # --- Interactive Plotly Equivalents ---
    if HAS_PLOTLY:
        try:
            fig_int = px.bar(
                df_sorted,
                x="ssgsea_enrichment_score_normalized",
                y="description",
                color="multi_evidence_pathway_score",
                color_continuous_scale="viridis",
                labels={"ssgsea_enrichment_score_normalized": "Normalized ssGSEA Score", "description": "Pathway", "multi_evidence_pathway_score": "Multi-Evidence Score"},
                title="Single-Sample Pathway Analysis (Interactive)",
                orientation="h"
            )
            fig_int.update_layout(template="plotly_white", margin=dict(l=20, r=20, t=40, b=40))
            bar_html = outdir / "sspa_ssgsea_barplot_interactive.html"
            pio.write_html(fig_int, str(bar_html))
            generated_files["sspa_barplot_html"] = bar_html
            
            fig_bub = px.scatter(
                df_sorted,
                x="multi_evidence_pathway_score",
                y="description",
                size="query_count_k",
                color="normalized_enrichment_score_nes",
                color_continuous_scale="coolwarm",
                labels={"multi_evidence_pathway_score": "Multi-Evidence Score", "normalized_enrichment_score_nes": "GSEA NES", "query_count_k": "Query Count"},
                title="Multi-Evidence Integrated Pathway Priority (Interactive)"
            )
            fig_bub.update_layout(template="plotly_white")
            bub_html = outdir / "sspa_multi_evidence_bubbleplot_interactive.html"
            pio.write_html(fig_bub, str(bub_html))
            generated_files["sspa_bubbleplot_html"] = bub_html
        except Exception as e:
            logger.warning(f"Plotly ssPA charts generation failed: {e}")

    return generated_files


def run_all_visualizations(
    output_dir: Path,
    config: Any
) -> Dict[str, Any]:
    """
    Orchestrates execution of the entire upgraded visualizer stage.
    """
    output_dir = Path(output_dir)
    vis_dir = output_dir / "visualizations"
    vis_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("=== Starting PathoScope AI Publication-Grade Scientific Visualization Stage ===")
    
    results = {}
    
    # Files paths
    qc_file = output_dir / "preprocessed" / "qc_report.json"
    gff_file = output_dir / "orfs" / "coordinates.gff3"
    annot_file = output_dir / "annotations" / "annotated_proteins.csv"
    mappings_file = output_dir / "pathways" / "mapped_pathways.csv"
    enrich_file = output_dir / "enrichment" / "enrichment_results.csv"
    volc_file = output_dir / "enrichment" / "volcano_plot_data.csv"
    matrix_file = output_dir / "enrichment" / "enrichment_matrix.csv"
    ranking_file = output_dir / "enrichment" / "pathway_ranking_reports.csv"
    
    # 0. Run Preprocessing QC plots
    import json
    if qc_file.exists():
        logger.info("Generating sequence quality and read retention QC charts...")
        try:
            with open(qc_file, "r", encoding="utf-8") as f:
                qc_data = json.load(f)
            qc_plots = generate_preprocessing_qc_plots(qc_data, vis_dir, config)
            results.update(qc_plots)
        except Exception as qcp_err:
            logger.warning(f"Failed to generate preprocessing QC plots: {qcp_err}")
            
    # 1. Run ORF distribution plots
    if gff_file.exists():
        logger.info("Generating ORF distribution and coordinate track plots...")
        orf_plots = generate_orf_distribution_plots(gff_file, vis_dir, config)
        results.update(orf_plots)
        
    # 2. Run annotation distributions
    if annot_file.exists():
        logger.info("Generating sequence alignment metrics distributions...")
        dist_plot = generate_annotation_distribution_plots(annot_file, vis_dir, config)
        if dist_plot:
            results["annotation_distributions"] = dist_plot
            
    # 3. Run pathway enrichment ORA charts (horizontal bar and bubble plots)
    if enrich_file.exists():
        logger.info("Generating pathway ORA enrichment bar & bubble charts...")
        enrich_plots = generate_enrichment_plots(enrich_file, vis_dir, config)
        results.update(enrich_plots)
        
    # 3.5 Run ssPA plots (single-sample pathway barplot and multi-evidence bubbleplot)
    if ranking_file.exists():
        logger.info("Generating ssPA single-sample pathway and multi-evidence priority charts...")
        try:
            sspa_plots = generate_sspa_plots(ranking_file, vis_dir, config)
            results.update(sspa_plots)
        except Exception as sspa_err:
            logger.warning(f"Failed to generate ssPA plots: {sspa_err}")

    # 4. Run pathway networks
    if mappings_file.exists() and enrich_file.exists():
        logger.info("Generating protein-pathway interaction network...")
        net_plot = generate_protein_pathway_network(mappings_file, enrich_file, vis_dir)
        if net_plot:
            results["protein_pathway_network"] = net_plot
            
    # 5. Run Volcano Plot
    if volc_file.exists():
        logger.info("Generating volcano plot of pathways significance...")
        volc_plots = generate_volcano_plot(volc_file, vis_dir, config)
        results.update(volc_plots)
        
    # 6. Run Gene-Pathway Heatmap
    if matrix_file.exists():
        logger.info("Generating gene-pathway profile heatmaps...")
        h_plot = generate_enrichment_heatmap(matrix_file, vis_dir, config)
        if h_plot:
            results["gene_pathway_heatmap"] = h_plot
            
    # 7. Run SVD-based PCA projection Plot
    if matrix_file.exists():
        logger.info("Generating SVD functional profile PCA projection maps...")
        p_plot = generate_pca_plot(matrix_file, vis_dir, config)
        if p_plot:
            results["functional_pca_plot"] = p_plot
            
    logger.info("=== Scientific Visualization Stage Completed Successfully ===")
    return {
        "status": "SUCCESS",
        "output_directory": str(vis_dir),
        "generated_files": {k: str(v) for k, v in results.items()}
    }
