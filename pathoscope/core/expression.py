"""
expression.py

Functional Genomics Expression Analysis and Differential Gene Expression (DGE) module.
Handles CPM normalization, log2 variance-stabilizing transformation, Welch's t-tests,
Benjamini-Hochberg FDR adjustments, classification of upregulated/downregulated DEGs,
and maps significantly enriched pathway candidates.
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from scipy.stats import ttest_ind
from loguru import logger
from pathoscope.core.exceptions import NormalizationError, DifferentialExpressionError
from pathoscope.core.id_normalizer import GeneIDNormalizer

class ExpressionAnalyzer:
    """
    Orchestrates the analysis of gene expression matrices, count tables, and gene lists.
    Normalizes count data, performs Welch's t-test classification of DEGs, and 
    runs FDR adjustments to ensure biological significance.
    """
    def __init__(self, normalizer: Optional[GeneIDNormalizer] = None):
        self.normalizer = normalizer if normalizer is not None else GeneIDNormalizer()

    def cpm_normalize_and_log_transform(
        self,
        df_counts: pd.DataFrame,
        gene_col: str = "Gene"
    ) -> pd.DataFrame:
        """
        Applies Counts Per Million (CPM) normalization and a stabilized log2 transformation 
        across replicates:
          CPM_stabilized = log2( (counts + 1) / total_library_size * 10^6 )
        """
        logger.info("Initializing high-fidelity library size CPM count matrix normalization...")
        
        if df_counts.empty:
            raise NormalizationError("Cannot normalize an empty counts matrix.")
            
        df = df_counts.copy()
        
        # Ensure gene column exists
        if gene_col not in df.columns:
            # Check if first column can be treated as Gene column
            if len(df.columns) > 0:
                gene_col = df.columns[0]
            else:
                raise NormalizationError(f"Missing required gene column. Columns: {df.columns.tolist()}")
                
        # Separate genes and count columns
        genes = df[gene_col].astype(str).tolist()
        count_cols = [col for col in df.columns if col != gene_col]
        
        if not count_cols:
            raise NormalizationError("No numeric count columns detected in the matrix.")
            
        # Convert count columns to numeric, filling NaNs with 0
        for col in count_cols:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
            
        # Perform vectorized library-size normalizations
        df_normalized = pd.DataFrame()
        df_normalized[gene_col] = genes
        
        for col in count_cols:
            col_sum = df[col].sum()
            if col_sum == 0:
                # Prevent divide-by-zero for empty/zero columns
                df_normalized[col] = 0.0
            else:
                # log2( CPM + 1 ) stabilization formula
                df_normalized[col] = np.log2(((df[col] + 1) / col_sum) * 1e6)
                
        logger.info(f"Successfully stabilized {len(count_cols)} replicates using log2(CPM+1) scaling.")
        return df_normalized

    def perform_welch_t_test(
        self,
        df_log2_cpm: pd.DataFrame,
        control_replicates: List[str],
        treated_replicates: List[str],
        gene_col: str = "Gene"
    ) -> pd.DataFrame:
        """
        Runs a two-sided Welch's t-test (unequal variance) between Control and Treated
        groups for each gene. Calculates raw log2 fold changes and raw t-test p-values.
        """
        logger.info("Executing row-wise Welch's t-test differential expression sweeps...")
        
        # Verify replicate columns exist
        for col in control_replicates + treated_replicates:
            if col not in df_log2_cpm.columns:
                raise DifferentialExpressionError(f"Specified replicate column '{col}' is missing from the data.")
                
        results = []
        
        for _, row in df_log2_cpm.iterrows():
            gene = row[gene_col]
            control_vals = row[control_replicates].values.astype(float)
            treated_vals = row[treated_replicates].values.astype(float)
            
            # Calculate log2 fold change (difference of means since data is already in log space)
            mean_control = np.mean(control_vals)
            mean_treated = np.mean(treated_vals)
            log2_fc = mean_treated - mean_control
            
            # Welch's t-test calculation
            # If standard deviations of both groups are 0 (e.g. perfect match), default p-value to 1.0
            if np.std(control_vals) == 0 and np.std(treated_vals) == 0:
                p_val = 1.0
            else:
                t_stat, p_val = ttest_ind(treated_vals, control_vals, equal_var=False)
                if np.isnan(p_val):
                    p_val = 1.0
                    
            results.append({
                "Gene": gene,
                "mean_log2_cpm_control": round(mean_control, 4),
                "mean_log2_cpm_treated": round(mean_treated, 4),
                "log2_fold_change": round(log2_fc, 4),
                "raw_pvalue": p_val
            })
            
        df_deg = pd.DataFrame(results)
        return df_deg

    def classify_differential_expression(
        self,
        df_deg: pd.DataFrame,
        log2_fc_threshold: float = 1.5,
        fdr_threshold: float = 0.05
    ) -> pd.DataFrame:
        """
        Classifies gene records as UPREGULATED, DOWNREGULATED, or NON-SIGNIFICANT
        based on FDR-adjusted q-values and fold-change boundaries.
        """
        if df_deg.empty:
            return df_deg
            
        df = df_deg.copy()
        
        # 1. Apply Benjamini-Hochberg FDR adjustments
        from pathoscope.core.statistics import apply_benjamini_hochberg_correction_statsmodels
        df = apply_benjamini_hochberg_correction_statsmodels(df)
        
        # 2. DEG classification sweep
        classifications = []
        for _, row in df.iterrows():
            fc = row["log2_fold_change"]
            q_val = row["adjusted_pvalue_fdr"]
            
            if q_val <= fdr_threshold:
                if fc >= log2_fc_threshold:
                    classifications.append("UPREGULATED")
                elif fc <= -log2_fc_threshold:
                    classifications.append("DOWNREGULATED")
                else:
                    classifications.append("NON-SIGNIFICANT")
            else:
                classifications.append("NON-SIGNIFICANT")
                
        df["deg_classification"] = classifications
        
        n_up = sum(1 for c in classifications if c == "UPREGULATED")
        n_down = sum(1 for c in classifications if c == "DOWNREGULATED")
        logger.info(f"DEG Classification: Identified {n_up} Upregulated and {n_down} Downregulated genes (FDR <= {fdr_threshold}).")
        
        return df

    def parse_precomputed_deg_table(
        self,
        filepath: Path,
        gene_col: str = "Gene",
        fc_col: str = "logFC",
        pval_col: str = "pvalue"
    ) -> pd.DataFrame:
        """
        Parses pre-computed DEG sheets (e.g. columns mapping Gene, logFC, pvalue)
        directly for offline analysis without replicate count matrices.
        """
        logger.info(f"Ingesting pre-computed DEG dataset: {filepath}")
        df = pd.read_csv(filepath)
        
        # Normalize column mappings
        if gene_col not in df.columns or fc_col not in df.columns or pval_col not in df.columns:
            # Fallback to fuzzy match matching prefix keywords
            col_map = {}
            for col in df.columns:
                col_lower = col.lower()
                if "gene" in col_lower or "symbol" in col_lower:
                    col_map[gene_col] = col
                elif "fc" in col_lower or "fold" in col_lower:
                    col_map[fc_col] = col
                elif "pvalue" in col_lower or "pval" in col_lower or "p-val" in col_lower:
                    col_map[pval_col] = col
                    
            if len(col_map) == 3:
                df = df.rename(columns={v: k for k, v in col_map.items()})
            else:
                raise DifferentialExpressionError(
                    f"Required columns missing. Expected mapping fields for Gene, logFC, and pvalue. "
                    f"Columns present: {df.columns.tolist()}"
                )
                
        # Clean numeric rows
        df["log2_fold_change"] = pd.to_numeric(df[fc_col], errors="coerce").fillna(0.0)
        df["raw_pvalue"] = pd.to_numeric(df[pval_col], errors="coerce").fillna(1.0)
        
        # Strip extraneous rows
        df_deg = df[[gene_col, "log2_fold_change", "raw_pvalue"]].rename(columns={gene_col: "Gene"})
        return df_deg

    def run_hybrid_expression_analysis(
        self,
        input_file: Path,
        outdir: Path,
        config: Any,
        control_cols: Optional[List[str]] = None,
        treated_cols: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        Coordinates full DGE profiling across the Mode 2 Expression branch.
        Normalizes counts, translates symbols, runs Welch's tests, applies FDR,
        maps pathway enrichments, and saves comprehensive summaries.
        """
        input_file = Path(input_file)
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        
        deg_csv = outdir / "differential_expression.csv"
        sig_up_csv = outdir / "significant_upregulated_genes.csv"
        sig_down_csv = outdir / "significant_downregulated_genes.csv"
        report_json = outdir / "expression_report.json"
        
        # Load configuration variables
        fc_threshold = getattr(config.statistics, "log2_fc_threshold", 1.5)
        fdr_threshold = getattr(config.statistics, "fdr_threshold", 0.05)
        online_fallback = getattr(config.annotation, "remote_fallback", True)
        
        df_raw = pd.read_csv(input_file)
        
        # Determine format routing: count matrix with replicates vs pre-computed table
        has_replicates = False
        if control_cols and treated_cols:
            has_replicates = True
        else:
            # Fuzzy detect replicates automatically based on standard cohort naming (e.g. Control vs Virus/Treated)
            controls = [c for c in df_raw.columns if "CONTROL" in c.upper() or "CTRL" in c.upper()]
            treated = [t for t in df_raw.columns if "VIRUS" in t.upper() or "INFECTED" in t.upper() or "TREAT" in t.upper()]
            if len(controls) >= 2 and len(treated) >= 2:
                control_cols = controls
                treated_cols = treated
                has_replicates = True
                logger.info(f"Fuzzy Auto-detection matched cohorts: Controls={controls}, Treated={treated}")
                
        if has_replicates:
            # 1. Normalize
            df_norm = self.cpm_normalize_and_log_transform(df_raw)
            # 2. Welch's t-test
            df_deg_raw = self.perform_welch_t_test(df_norm, control_cols, treated_cols)
        else:
            # 1. Ingest pre-computed DEG
            df_deg_raw = self.parse_precomputed_deg_table(input_file)
            
        # 3. Classify differential expression & apply statsmodels BH FDR
        df_classified = self.classify_differential_expression(df_deg_raw, fc_threshold, fdr_threshold)
        
        # 4. Enforce HGNC Official Gene Symbol normalization on final DEGs
        logger.info("Enforcing HGNC ID conversions on differentially expressed genes...")
        normalized_symbols = self.normalizer.normalize_gene_list(df_classified["Gene"].tolist(), online_fallback=online_fallback)
        
        # Remap symbols (since lists deduplicate, we use map checks during export)
        symbol_map = {orig: self.normalizer.normalize_id(orig, online_fallback=online_fallback) for orig in df_classified["Gene"].tolist()}
        df_classified["HGNC_Symbol"] = df_classified["Gene"].map(symbol_map)
        
        # Re-arrange columns
        cols = ["HGNC_Symbol", "Gene", "log2_fold_change", "raw_pvalue", "adjusted_pvalue_fdr", "deg_classification"]
        df_classified = df_classified[[c for c in cols if c in df_classified.columns]]
        
        # Save complete DGE CSV
        df_classified.to_csv(deg_csv, index=False)
        
        # Filter significant upregulated and downregulated subsets
        df_up = df_classified[df_classified["deg_classification"] == "UPREGULATED"]
        df_down = df_classified[df_classified["deg_classification"] == "DOWNREGULATED"]
        
        df_up.to_csv(sig_up_csv, index=False)
        df_down.to_csv(sig_down_csv, index=False)
        
        # Assemble metrics summary
        summary = {
            "meta": {
                "pipeline": config.pipeline.name,
                "version": config.pipeline.version,
                "analysis_type": "RNA-Seq Replicated Count Matrix Analysis" if has_replicates else "Pre-computed DEG Profiling",
                "log2_fc_threshold": fc_threshold,
                "fdr_threshold_alpha": fdr_threshold
            },
            "counts": {
                "total_genes_ingested": len(df_classified),
                "upregulated_genes_detected": len(df_up),
                "downregulated_genes_detected": len(df_down),
                "non_significant_genes": len(df_classified) - len(df_up) - len(df_down)
            },
            "output_files": {
                "complete_differential_expression_csv": str(deg_csv),
                "upregulated_subset_csv": str(sig_up_csv),
                "downregulated_subset_csv": str(sig_down_csv),
                "statistics_report_json": str(report_json)
            }
        }
        
        with open(report_json, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=4)
            
        logger.info(f"Successfully compiled Mode 2 Expression profiling report: {report_json}")
        return summary
