import json
import time
from pathlib import Path
from typing import Dict, Any, Optional
from loguru import logger
from pathoscope.utils.config_loader import AppConfig, load_config
from pathoscope.utils.logger import setup_logger
from pathoscope.core.preprocessor import process_sequences
from pathoscope.core.orf_predictor import process_orf_prediction

class PipelineCoordinator:
    """
    Manages the overall execution flow, directory structure,
    checkpoints, and metadata logging of the PathoScope AI pipeline.
    """
    def __init__(self, config_path: Path, output_dir: Path):
        self.config_path = Path(config_path)
        self.config: AppConfig = load_config(self.config_path)
        self.raw_output_dir = Path(output_dir)
        
        # Formulate versioned output directory
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        dir_name = f"run_{self.config.pipeline.version}_{timestamp}"
        self.output_dir = self.raw_output_dir / dir_name
        self.temp_dir = self.output_dir / "temp"
        
        # Track pipeline stage execution times
        self.execution_times: Dict[str, float] = {}

    def initialize_workspace(self):
        """Creates the versioned folder hierarchy."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize logging inside the output directory
        log_file = self.output_dir / "pipeline.log"
        setup_logger(log_file=log_file)
        logger.info("PathoScope AI Workspace initialized successfully.")
        logger.info(f"Target execution output directory: {self.output_dir}")

    def run_preprocessing(self, input_file: Path) -> Path:
        """Runs the sequence validation and cleaning stage."""
        input_file = Path(input_file)
        start_time = time.time()
        logger.info("=== Starting Preprocessing & QC Phase ===")
        
        output_fasta = self.output_dir / "preprocessed" / "cleaned.fasta"
        qc_json_path = self.output_dir / "preprocessed" / "qc_report.json"
        
        # Execute preprocessor
        qc_summary = process_sequences(input_file, output_fasta, self.config)
        
        # Save structured QC report
        qc_json_path.parent.mkdir(parents=True, exist_ok=True)
        with open(qc_json_path, "w", encoding="utf-8") as f:
            json.dump(qc_summary, f, indent=4)
            
        elapsed = time.time() - start_time
        self.execution_times["preprocessing"] = round(elapsed, 4)
        
        logger.info(f"QC Summary: Kept {qc_summary['counts']['total_kept']} / {qc_summary['counts']['total_processed']} sequences.")
        logger.info(f"Preprocessing completed in {elapsed:.2f} seconds.")
        return output_fasta

    def run_orf_prediction(self, cleaned_fasta: Path) -> Path:
        """Runs the 6-frame ORF scanning and translation stage."""
        cleaned_fasta = Path(cleaned_fasta)
        start_time = time.time()
        logger.info("=== Starting ORF Prediction & Translation Phase ===")
        
        orf_outdir = self.output_dir / "orfs"
        
        # Execute ORF predictor
        orf_stats = process_orf_prediction(cleaned_fasta, orf_outdir, self.config)
        
        elapsed = time.time() - start_time
        self.execution_times["orf_prediction"] = round(elapsed, 4)
        
        logger.info(f"ORF Summary: Predicted {orf_stats['counts']['total_orfs_predicted']} ORFs.")
        logger.info(f"ORF prediction completed in {elapsed:.2f} seconds.")
        return orf_outdir / "proteins.fasta"

    def run_annotation(self, proteins_fasta: Path) -> Path:
        """Runs the functional sequence similarity annotation stage."""
        proteins_fasta = Path(proteins_fasta)
        start_time = time.time()
        logger.info("=== Starting Sequence Similarity Annotation Phase ===")
        
        # Lazy load functional annotation to optimize startup latency
        from pathoscope.core.annotator import process_functional_annotation
        
        annotation_outdir = self.output_dir / "annotations"
        
        # Execute Annotator
        annotation_stats = process_functional_annotation(proteins_fasta, annotation_outdir, self.config)
        
        elapsed = time.time() - start_time
        self.execution_times["functional_annotation"] = round(elapsed, 4)
        
        logger.info(f"Annotation Summary: Annotated {annotation_stats['counts']['total_annotated']} proteins ({annotation_stats['counts']['total_hypothetical']} hypothetical).")
        logger.info(f"Functional annotation completed in {elapsed:.2f} seconds.")
        return annotation_outdir / "annotated_proteins.csv"

    def run_pathway_and_domain_mapping(self, proteins_fasta: Path, annotated_csv: Path) -> Path:
        """Runs the conserved domain search (HMMER/Pfam) and dynamic pathway mapping (KEGG/Reactome) stage."""
        proteins_fasta = Path(proteins_fasta)
        annotated_csv = Path(annotated_csv)
        start_time = time.time()
        logger.info("=== Starting Conserved Domain Search & Pathway Mapping Phase ===")
        
        # Lazy load pathway mapper to optimize startup latency
        from pathoscope.core.pathway_mapper import process_pathway_and_domain_mapping
        
        mapping_outdir = self.output_dir / "pathways"
        
        # Execute Domain & Pathway Mapper
        mapping_stats = process_pathway_and_domain_mapping(
            proteins_fasta, annotated_csv, mapping_outdir, self.config
        )
        
        elapsed = time.time() - start_time
        self.execution_times["pathway_and_domain_mapping"] = round(elapsed, 4)
        
        logger.info(f"Pathway Mapping Summary: Mapped {mapping_stats['pathways']['total_proteins_mapped']} proteins to {mapping_stats['pathways']['unique_pathways_implicated']} unique pathways. Detected {mapping_stats['domains']['total_domains_detected']} Pfam domains.")
        logger.info(f"Pathway and domain mapping completed in {elapsed:.2f} seconds.")
        return mapping_outdir / "mapped_pathways.csv"

    def run_enrichment_analysis(self, mapped_pathways_csv: Path) -> Path:
        """Runs the statistically valid pathway enrichment analysis stage."""
        mapped_pathways_csv = Path(mapped_pathways_csv)
        start_time = time.time()
        logger.info("=== Starting Hypergeometric Pathway Enrichment Analysis Phase ===")
        
        # Lazy load statistics engine to optimize startup latency
        from pathoscope.core.statistics import process_pathway_enrichment
        
        enrichment_outdir = self.output_dir / "enrichment"
        
        # Execute Enrichment Analyzer
        enrich_stats = process_pathway_enrichment(mapped_pathways_csv, enrichment_outdir, self.config)
        
        elapsed = time.time() - start_time
        self.execution_times["pathway_enrichment"] = round(elapsed, 4)
        
        logger.info(f"Pathway Enrichment Summary: Tested {enrich_stats['counts']['total_pathways_tested']} pathways. Found {enrich_stats['counts']['significant_pathways_enriched']} statistically significant enriched pathways (FDR <= {self.config.statistics.fdr_threshold}).")
        logger.info(f"Pathway enrichment completed in {elapsed:.2f} seconds.")
        return enrichment_outdir / "significant_pathways.csv"

    def run_scientific_visualizations(self) -> Dict[str, Any]:
        """Orchestrates generation of all scientific visualizations (plots and bipartite networks)."""
        start_time = time.time()
        logger.info("=== Starting Scientific Visualizations Phase ===")
        
        # Lazy load visualizer to optimize startup latency
        from pathoscope.reporting.visualizer import run_all_visualizations
        
        vis_stats = run_all_visualizations(self.output_dir, self.config)
        
        elapsed = time.time() - start_time
        self.execution_times["scientific_visualizations"] = round(elapsed, 4)
        
        logger.info(f"Scientific visualizations completed in {elapsed:.2f} seconds. Generated {len(vis_stats.get('generated_files', {}))} publication-quality figures.")
        return vis_stats

    def run_ai_synthesis(self) -> Dict[str, Any]:
        """Runs the context-grounded AI biological interpretation layer."""
        start_time = time.time()
        logger.info("=== Starting AI Biological Interpretation Phase ===")
        
        # Lazy load AI interpreter to optimize startup latency
        from pathoscope.interpretation.ai_interpreter import run_ai_biological_interpretation
        
        ai_stats = run_ai_biological_interpretation(self.output_dir, self.config)
        
        elapsed = time.time() - start_time
        self.execution_times["ai_synthesis"] = round(elapsed, 4)
        
        logger.info(f"AI biological interpretation completed in {elapsed:.2f} seconds.")
        return ai_stats

    def run_report_generation_stage(self) -> Dict[str, Any]:
        """Orchestrates generation of HTML dashboard and PDF analytical report."""
        start_time = time.time()
        logger.info("=== Starting Report Generation Phase ===")
        
        # Lazy load PDF reporter to optimize startup latency
        from pathoscope.reporting.reporter import run_report_generation
        
        report_stats = run_report_generation(self.output_dir)
        
        elapsed = time.time() - start_time
        self.execution_times["report_generation"] = round(elapsed, 4)
        
        logger.info(f"Report generation completed in {elapsed:.2f} seconds. HTML Dashboard saved: {report_stats['html_report']}")
        return report_stats

    def write_run_metadata(self, input_file: Path, status: str = "SUCCESS"):
        """Saves execution metadata to metadata.json for audit and reproducibility."""
        metadata_path = self.output_dir / "metadata.json"
        
        # Compute input file MD5 hash if exists
        input_md5 = "UNKNOWN"
        try:
            if input_file.exists():
                h = hashlib.md5()
                with open(input_file, "rb") as f:
                    for chunk in iter(lambda: f.read(4096), b""):
                        h.update(chunk)
                input_md5 = h.hexdigest()
        except Exception:
            pass

        import hashlib
        import sys
        
        metadata = {
            "pipeline": {
                "name": self.config.pipeline.name,
                "version": self.config.pipeline.version,
                "schema_version": "1.0.0"
            },
            "environment": {
                "python_version": sys.version,
                "platform": sys.platform
            },
            "run": {
                "status": status,
                "input_file": str(input_file),
                "input_md5": input_md5,
                "output_dir": str(self.output_dir),
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                "command_line": "python " + " ".join(sys.argv),
                "execution_times_seconds": self.execution_times
            },
            "config_snapshot": self.config.model_dump()
        }
        
        with open(metadata_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=4)
        logger.info(f"Run metadata written to: {metadata_path}")

    def clean_temporary_files(self):
        """Safely removes the temp directory to clean up space."""
        try:
            if self.temp_dir.exists():
                # Recursively delete files
                for child in self.temp_dir.iterdir():
                    if child.is_file():
                        child.unlink()
                self.temp_dir.rmdir()
                logger.info("Temporary workspace directory cleaned up.")
        except Exception as e:
            logger.warning(f"Failed to clean up temporary workspace: {e}")

    def run_expression_profiling(self, input_file: Path) -> Path:
        """
        Coordinates the Mode 2 Differential Gene Expression (DGE) and count normalization profiling.

        Purpose:
            To ingest RNA-Seq count matrices or gene expression sheets, apply CPM depth normalizations, 
            run variance-stabilizing log2 transformations, perform Welch's t-test calculations across replicates, 
            and classify significant upregulated and downregulated Differentially Expressed Genes (DEGs).

        Inputs:
            input_file (Path): Absolute or relative path to the count matrix or expression CSV/TSV table.

        Outputs:
            Path: Path to the generated standard 'differential_expression.csv' table containing normalized values,
                  fold changes, raw p-values, FDR adjusted q-values, and classifications.

        Biological Rationale:
            Transcriptome profiling requires normalization for sequencing depth (CPM library size scaling) and 
            log-transformation to stabilize expression variance across replicates before statistical DGE evaluation. 
            Welch's t-test is biologically robust because it assumes unequal variances between replicates, 
            preventing false positives due to varying sequencing depth or transcript dispersion.
        """
        input_file = Path(input_file)
        start_time = time.time()
        logger.info("=== Starting Mode 2: Differential Gene Expression Phase ===")
        
        # Lazy load ExpressionAnalyzer to optimize startup latency
        from pathoscope.core.expression import ExpressionAnalyzer
        
        expression_outdir = self.output_dir / "expression"
        analyzer = ExpressionAnalyzer()
        
        # Ingest cohorts if configured, otherwise fuzzy auto-detect replicates inside the analyzer
        control_cols = getattr(self.config.statistics, "control_replicates", None)
        treated_cols = getattr(self.config.statistics, "treated_replicates", None)
        
        # Execute DGE analysis
        expr_stats = analyzer.run_hybrid_expression_analysis(
            input_file, expression_outdir, self.config, control_cols, treated_cols
        )
        
        elapsed = time.time() - start_time
        self.execution_times["differential_expression"] = round(elapsed, 4)
        
        logger.info(f"DGE Summary: Ingested {expr_stats['counts']['total_genes_ingested']} genes. Detected {expr_stats['counts']['upregulated_genes_detected']} upregulated and {expr_stats['counts']['downregulated_genes_detected']} downregulated DEGs.")
        logger.info(f"Differential Gene Expression completed in {elapsed:.2f} seconds.")
        
        # Remap output path for downstream statistics enrichment
        return expression_outdir / "differential_expression.csv"

    def execute_pipeline(self, input_file: Path):
        """
        Coordinates full end-to-end execution, dynamically routing inputs.

        Purpose:
            To automatically route raw inputs (sequences or expression tables) through their appropriate 
            computational workflows: Mode 1 (Viral Sequence Analysis) or Mode 2 (Functional Genomics Expression Analysis), 
            and generate complete analytical reports.

        Inputs:
            input_file (Path): Path to the input file (FASTA, FASTQ, CSV, TSV, or TXT).

        Outputs:
            None: Compiles GFF3 files, intermediate tables, Plotly figures, and MultiQC-style HTML/PDF reports on disk.

        Biological Rationale:
            Biological sequence analysis and host transcriptomic profiling are fundamentally distinct processes. 
            Sequence analysis seeks to identify coding regions and evolutionary homologs (Mode 1), while expression 
            analysis evaluates host-system cellular hijacks and pathway shifts (Mode 2). Routing inputs dynamically 
            preserves mathematical and biological validity by applying the correct statistical workflows.
        """
        input_file = Path(input_file)
        self.initialize_workspace()
        
        status = "FAILED"
        try:
            # Detect pipeline mode based on file suffixes and headers
            suffix = input_file.suffix.lower()
            is_expression = False
            
            if suffix in [".csv", ".tsv", ".txt"]:
                # Check headers to fuzzy-match expression matrix/DEG files vs sequence tables
                try:
                    df_check = pd.read_csv(input_file, nrows=5)
                    # Check for logFC, counts, replication keywords, or gene lists without nucleotides
                    cols_lower = [c.lower() for c in df_check.columns]
                    has_sequence_col = any("seq" in c or "nuc" in c for c in cols_lower)
                    
                    if not has_sequence_col:
                        # Exclude sequence tables
                        has_fc = any("fc" in c or "fold" in c or "pval" in c for c in cols_lower)
                        has_counts = any("control" in c or "ctrl" in c or "virus" in c or "infect" in c for c in cols_lower)
                        is_expression = has_fc or has_counts or len(df_check.columns) <= 3
                except Exception:
                    pass
            
            if is_expression:
                logger.info(f"Input '{input_file.name}' routed to Mode 2: Functional Genomics branch.")
                
                # 1. Run Differential Gene Expression profiling (includes CPM & Welch's tests)
                deg_csv = self.run_expression_profiling(input_file)
                
                # 2. Run Statistics ORA & ssGSEA Enrichment
                # To maintain compatibility, we redirect the mapping table to DGE outputs
                significant_pathways_csv = self.run_enrichment_analysis(deg_csv)
                
                # 3. Run Visualizations (generates responsive Volcano, Bubble, and PCA Plots)
                vis_stats = self.run_scientific_visualizations()
                
                # 4. Run RAG-Grounded AI literature explanations
                ai_stats = self.run_ai_synthesis()
                
            else:
                logger.info(f"Input '{input_file.name}' routed to Mode 1: Viral Sequence branch.")
                
                # 1. Run Preprocessing & QC
                cleaned_fasta = self.run_preprocessing(input_file)
                
                # 2. Run ORF Prediction & Translation
                proteins_fasta = self.run_orf_prediction(cleaned_fasta)
                
                # 3. Run Swiss-Prot Functional similarity search
                annotated_csv = self.run_annotation(proteins_fasta)
                
                # 4. Run Conserved domain (HMMER) and dynamic pathway mapping
                mapped_pathways_csv = self.run_pathway_and_domain_mapping(proteins_fasta, annotated_csv)
                
                # 5. Run Pathway Enrichment
                significant_pathways_csv = self.run_enrichment_analysis(mapped_pathways_csv)
                
                # 6. Run Visualizations
                vis_stats = self.run_scientific_visualizations()
                
                # 7. Run AI biological interpretation
                ai_stats = self.run_ai_synthesis()
                
            # Write run metadata report footprint
            self.write_run_metadata(input_file, status="SUCCESS")
            
            # 8. Run HTML dashboard and PDF compilation Exporter
            report_stats = self.run_report_generation_stage()
            
            status = "SUCCESS"
            
        except Exception as e:
            logger.error(f"Pipeline execution encountered an unhandled error: {e}")
            raise e
            
        finally:
            self.write_run_metadata(input_file, status=status)
            self.clean_temporary_files()
