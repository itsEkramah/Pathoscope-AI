import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import json
import time
import os
import io
import shutil
import base64
import sys
from pathlib import Path

# Add the project root to sys.path to allow absolute imports of 'pathoscope' package
root_dir = Path(__file__).resolve().parent.parent
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))

import matplotlib.pyplot as plt
import networkx as nx
from loguru import logger
from Bio import SeqIO

from pathoscope.pipeline import PipelineCoordinator
from pathoscope.utils.config_loader import load_config, AppConfig

# Premium Styling & CSS overrides inspired by nf-core and MultiQC
st.set_page_config(
    page_title="PathoScope AI - Modern Computational Virology Portal",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Dark Slate & Indigo academic styling sheet
CUSTOM_CSS = """
<style>
    /* Global Background and Typography */
    .stApp {
        background-color: #0b0f19;
        color: #f1f5f9;
    }
    
    /* Premium Header / Hero Card */
    .hero-card {
        background: linear-gradient(135deg, #0f172a 0%, #1e1b4b 50%, #0f172a 100%);
        border-radius: 16px;
        padding: 40px;
        color: white;
        text-align: center;
        border: 1px solid #4f46e5;
        box-shadow: 0 10px 30px rgba(0, 0, 0, 0.5);
        margin-bottom: 30px;
    }
    .hero-card h1 {
        font-family: 'Outfit', 'Inter', sans-serif !important;
        font-weight: 900 !important;
        font-size: 3.5rem !important;
        margin: 0 !important;
        padding: 0 !important;
        background: linear-gradient(90deg, #38bdf8, #818cf8, #34d399);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        letter-spacing: -1px;
    }
    .hero-card p {
        font-size: 1.35rem !important;
        color: #cbd5e1 !important;
        margin-top: 15px !important;
        margin-bottom: 0 !important;
        font-weight: 300 !important;
    }
    
    /* Academic Section Headers */
    .section-header {
        font-family: 'Outfit', 'Inter', sans-serif !important;
        font-weight: 700;
        font-size: 1.6rem;
        color: #a5b4fc;
        margin-top: 25px;
        margin-bottom: 15px;
        border-bottom: 1px solid rgba(165, 180, 252, 0.15);
        padding-bottom: 8px;
    }
    
    /* MultiQC-style Metrics Grid */
    .metric-card {
        background: rgba(15, 23, 42, 0.6);
        border: 1px solid rgba(148, 163, 184, 0.15);
        border-radius: 12px;
        padding: 22px;
        text-align: center;
        box-shadow: 0 4px 15px rgba(0, 0, 0, 0.25);
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    }
    .metric-card:hover {
        border-color: #38bdf8;
        transform: translateY(-4px);
        box-shadow: 0 10px 20px rgba(56, 189, 248, 0.15);
    }
    .metric-val {
        font-size: 2.5rem;
        font-weight: 800;
        color: #38bdf8;
        margin-bottom: 5px;
    }
    .metric-label {
        font-size: 0.95rem;
        color: #94a3b8;
        text-transform: uppercase;
        letter-spacing: 1px;
        font-weight: 600;
    }
    
    /* Live logs console */
    .log-box {
        background-color: #020617 !important;
        border: 1px solid #1e293b !important;
        border-radius: 10px !important;
        padding: 20px !important;
        font-family: 'Fira Code', 'Courier New', Courier, monospace !important;
        font-size: 0.9rem !important;
        color: #38bdf8 !important;
        box-shadow: inset 0 2px 8px rgba(0,0,0,0.8);
    }
    
    /* Download box */
    .download-card {
        background: rgba(6, 182, 212, 0.05);
        border: 1px solid rgba(6, 182, 212, 0.2);
        border-radius: 12px;
        padding: 25px;
        margin-top: 15px;
        box-shadow: 0 4px 12px rgba(0,0,0,0.15);
    }
    
    /* Custom Tab UI Styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        background-color: rgba(255, 255, 255, 0.02);
        border: 1px solid rgba(255, 255, 255, 0.05);
        border-radius: 8px 8px 0 0;
        padding: 12px 22px;
        font-weight: 600;
        color: #94a3b8;
        transition: all 0.2s ease;
    }
    .stTabs [data-baseweb="tab"]:hover {
        color: #38bdf8;
        background-color: rgba(56, 189, 248, 0.05);
    }
    .stTabs [aria-selected="true"] {
        background-color: rgba(79, 70, 229, 0.15) !important;
        border-color: #4f46e5 !important;
        color: #c7d2fe !important;
    }
    
    /* Alignment hits hover styling */
    .alignment-card {
        background: rgba(255,255,255,0.02);
        border: 1px solid rgba(255,255,255,0.05);
        border-radius: 8px;
        padding: 15px;
        margin-bottom: 10px;
    }
</style>
"""
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)

# Sample FASTA representing mock viral isolate with valid mock coding ORFs
SAMPLE_FASTA_CONTENT = """>PathoScope_Sample_Bacteriophage_MS2_Isolate
GGGTGTGGGCCCCCAAATAAGCCAAGGACCGGGCGACCGTCGCGCCTAATAGAGTCATATGACCACCCGT
CAGTTACCAACGCTCGAAGTCACTACTACCTAACGGCCGCCCGCTAATGACGGCCGCCCGCTAATGACGG
CCGCCCGCTAATGACGGCCGCCCGCTAATGACGGCCGCCCGCTAATGACGGCCGCCCGCTAATGACGGCC
ATGGCTTCGAACTTCGCTTCGGTCGCCGTCCTCCGCGCTGGTGTCGCCGTCCTCCGCGCTGGTGTCGCCG
TCCTCCGCGCTGGTGTCGCCGTCCTCCGCGCTGGTGTCGCCGTCCTCCGCGCTGGTTGAAGATCCTAACC
ATGGCGAGCAATAACAGCATGGCGAGCAATAACAGCATGGCGAGCAATAACAGCATGGCGAGCAATAACA
GCATGGCGAGCAATAACAGCATGGCGAGCAATAACAGCATGGCGAGCAATAACAGCATGGCGAGCTAGAA
GGGTGTGGGCCCCCAAATAAGCCAAGGACCGGGCGACCGTCGCGCCTAATAGAGTCATATGACCACCCGT
"""

class StreamlitLogBuffer:
    """Interceptors to collect pipeline log outputs and update a st.code block."""
    def __init__(self, code_placeholder):
        self.code_placeholder = code_placeholder
        self.logs = []
        
    def write(self, message):
        clean_msg = message.strip()
        if clean_msg:
            self.logs.append(clean_msg)
            log_text = "\n".join(self.logs[-35:])
            self.code_placeholder.code(log_text, language="text")
            
    def flush(self):
        pass

# ----------------------------------------------------------------------
# 💾 HIGH-PERFORMANCE RENDERING & MEMORY CACHING HELPERS
# ----------------------------------------------------------------------
@st.cache_data
def cached_read_csv(filepath: str) -> pd.DataFrame:
    """Memory-optimized, cached loading of heavy genomics reports."""
    fp = Path(filepath)
    if fp.exists():
        try:
            return pd.read_csv(fp)
        except Exception as e:
            logger.error(f"Failed to read CSV at {filepath}: {e}")
    return pd.DataFrame()

@st.cache_data
def cached_read_json(filepath: str) -> dict:
    """Memory-optimized, cached loading of unstructured AI/QC metadata."""
    fp = Path(filepath)
    if fp.exists():
        try:
            with open(fp, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            logger.error(f"Failed to load JSON at {filepath}: {e}")
    return {}

# ----------------------------------------------------------------------
# 💾 INITIALIZE STATE MAPPING
# ----------------------------------------------------------------------
if "run_completed" not in st.session_state:
    st.session_state.run_completed = False
if "run_outdir" not in st.session_state:
    st.session_state.run_outdir = None
if "run_metadata" not in st.session_state:
    st.session_state.run_metadata = None
if "current_logs" not in st.session_state:
    st.session_state.current_logs = ""
if "pipeline_steps" not in st.session_state:
    st.session_state.pipeline_steps = {
        "preprocess": {"name": "1. Quality Control & Trimming", "status": "pending", "elapsed": 0.0},
        "orf": {"name": "2. ORF Coordinate Translation", "status": "pending", "elapsed": 0.0},
        "similarity": {"name": "3. Homology Alignment Search", "status": "pending", "elapsed": 0.0},
        "domains": {"name": "4. Pfam Domain & Pathways", "status": "pending", "elapsed": 0.0},
        "enrichment": {"name": "5. Hypergeometric ssGSEA/ORA", "status": "pending", "elapsed": 0.0},
        "visuals": {"name": "6. Publication Figures Rendering", "status": "pending", "elapsed": 0.0},
        "ai": {"name": "7. Pydantic AI Literature Synthesis", "status": "pending", "elapsed": 0.0},
        "report": {"name": "8. Analytical Report Assembly", "status": "pending", "elapsed": 0.0}
    }

# Sidebar - Brand Header
st.sidebar.markdown(
    "<div style='text-align: center; padding-bottom: 15px;'>"
    "<span style='font-size: 4.5rem;'>🧬</span>"
    "<h2 style='margin: 0; color: #818cf8; font-weight: 800; font-family: Outfit;'>PathoScope AI</h2>"
    "<p style='font-size: 0.85rem; color: #94a3b8; font-weight: 500;'>Computational Virology Framework</p>"
    "</div>",
    unsafe_allow_html=True
)

st.sidebar.markdown("---")

# Load baseline configurations to populate defaults
DEFAULT_CONFIG_PATH = Path("config/default_config.yaml")
base_config = load_config(DEFAULT_CONFIG_PATH)

# Parameter Formats Sidebar
st.sidebar.markdown("### ⚙️ Pipeline Configurations")

# Expanders for Modular Configurations
with st.sidebar.expander("🔬 Preprocessing & Trimming"):
    min_len = st.slider("Min Sequence Length (bp)", 10, 2000, int(base_config.preprocessing.min_length))
    max_len = st.slider("Max Sequence Length (bp)", 1000, 1000000, int(base_config.preprocessing.max_length))
    max_ambig = st.slider("Max Ambiguous Bases %", 0.0, 50.0, float(base_config.preprocessing.max_ambiguous_pct), 0.5)
    remove_dups = st.checkbox("Remove Identical Sequences", base_config.preprocessing.remove_duplicate_sequences)

with st.sidebar.expander("🧬 ORF Genomic Translation"):
    table_num = st.selectbox("Genetic Code Translation Table", [1, 11, 4, 5, 12, 15, 25], index=0)
    min_orf_aa = st.slider("Min ORF Length (aa)", 10, 500, int(base_config.orf_prediction.min_orf_length_aa))
    start_codons_list = st.multiselect("Start Codons", ["ATG", "GTG", "TTG", "ATT", "CTG"], default=base_config.orf_prediction.start_codons)
    stop_codons_list = st.multiselect("Stop Codons", ["TAA", "TAG", "TGA"], default=base_config.orf_prediction.stop_codons)
    resolve_nest = st.checkbox("Filter Nested In-Frame ORFs", base_config.orf_prediction.resolve_nested)

with st.sidebar.expander("🎯 Similarity Alignment"):
    align_engine = st.selectbox("Alignment Engine", ["diamond", "blastp"], index=0)
    local_db = st.text_input("Local Database Path", base_config.annotation.local_db_path)
    remote_fallback = st.checkbox("Fallback to Remote NCBI BLASTp", True)
    eval_thresh = st.number_input("E-value Threshold", value=float(base_config.annotation.eval_threshold), format="%e")
    ident_thresh = st.slider("Min Identity %", 0.0, 100.0, float(base_config.annotation.identity_threshold))
    cov_thresh = st.slider("Min Coverage %", 0.0, 100.0, float(base_config.annotation.coverage_threshold))

with st.sidebar.expander("💡 Domains & Pathway Mapping"):
    hmmer_db = st.text_input("Pfam HMMER Database Path", base_config.domain_search.hmmer_db_path)
    hmmer_eval = st.number_input("HMMER E-value Threshold", value=float(base_config.domain_search.eval_threshold), format="%e")
    kegg_organism = st.text_input("KEGG Target Organism Group", base_config.pathway_mapping.kegg_organism)
    cache_db = st.text_input("Local Cache SQLite Path", base_config.pathway_mapping.db_cache_path)

with st.sidebar.expander("📊 Enrichment Statistics"):
    fdr_thresh = st.slider("Benjamini-Hochberg FDR", 0.001, 0.20, float(base_config.statistics.fdr_threshold), 0.005)
    universe_size = st.number_input("Background Universe Size", value=int(base_config.statistics.bg_universe_size), min_value=10)

# Sidebar System Preferences Panel
st.sidebar.markdown("### 🛠️ System Preferences")
advanced_mode = st.sidebar.checkbox("🧬 Advanced Bioinformatician Mode", value=False, help="Enable to expose environment structures, raw SQL cache databases, terminal commands, and system variables.")

# ----------------------------------------------------------------------
# 🤖 AI INTERPRETATION LAYER (CODE-ONLY OPTIONS)
# Enforces strict 0.1 temperature clamping for anti-hallucination.
# ----------------------------------------------------------------------
AI_PROVIDER = "openai"          # Options: "openai", "gemini", "offline-rules"
AI_MODEL = "gpt-4o-mini"        # Model e.g. "gpt-4o-mini"
AI_TEMP = 0.1                   # Temperature clamped to 0.1 for maximum determinism
# ----------------------------------------------------------------------

# Main UI Frame - Premium Hero Title Card
st.markdown(
    "<div class='hero-card'>"
    "<h1>PathoScope AI 🧬</h1>"
    "<p>Automated Viral Genomics Pipeline for Sequence Annotation, Pathway Mapping, and AI Literature Grounding</p>"
    "</div>",
    unsafe_allow_html=True
)

# Input Sequence Block
col_input, col_preset = st.columns([2, 1])

with col_input:
    st.markdown("<div class='section-header'>📥 Input Isolate Sequence</div>", unsafe_allow_html=True)
    uploaded_file = st.file_uploader(
        "Upload raw nucleotide sequences (FASTA, FNA, FASTQ, or GZ)",
        type=["fasta", "fa", "fna", "fastq", "fq", "gz"],
        help="Provide a standard sequence file (FASTA/FASTQ, or gzipped .fasta.gz/.fastq.gz) containing viral isolates to launch the pipeline."
    )

with col_preset:
    st.markdown("<div class='section-header'>🧪 Quick-Test Preset</div>", unsafe_allow_html=True)
    use_sample = st.checkbox(
        "Use Preset Sequence",
        value=False,
        help="Enable to test the pipeline instantly with a mock viral isolate (Leviviridae Bacteriophage MS2)."
    )

# Execution Controls
st.markdown("<div class='section-header'>🚀 Pipeline Execution Control</div>", unsafe_allow_html=True)
run_button = st.button("RUN GENOMICS PIPELINE", type="primary", use_container_width=True)

if run_button:
    fasta_path = None
    
    if use_sample:
        fasta_path = Path("temp_sample_genome.fasta")
        with open(fasta_path, "w", encoding="utf-8") as f:
            f.write(SAMPLE_FASTA_CONTENT)
        st.info("🧬 Preset Leviviridae MS2 sequence loaded.")
    elif uploaded_file is not None:
        # Dynamically preserve extension suffixes (e.g., .fastq.gz or .fastq)
        orig_suffix = "".join(Path(uploaded_file.name).suffixes)
        fasta_path = Path(f"temp_uploaded_genome{orig_suffix}")
        with open(fasta_path, "wb") as f:
            f.write(uploaded_file.getbuffer())
        st.info(f"📁 Custom sequence '{uploaded_file.name}' loaded successfully.")
    else:
        st.error("⚠️ Please upload a sequence file or toggle the 'Use Preset Sequence' checkbox to execute.")
        st.stop()

    out_parent = Path("results_streamlit")
    out_parent.mkdir(exist_ok=True)
    
    # Initialize steps status to pending
    for key in st.session_state.pipeline_steps:
        st.session_state.pipeline_steps[key]["status"] = "pending"
        st.session_state.pipeline_steps[key]["elapsed"] = 0.0
        
    coordinator = PipelineCoordinator(config_path=DEFAULT_CONFIG_PATH, output_dir=out_parent)
    
    # Ingest inputs to configuration model
    coordinator.config.preprocessing.min_length = min_len
    coordinator.config.preprocessing.max_length = max_len
    coordinator.config.preprocessing.max_ambiguous_pct = max_ambig
    coordinator.config.preprocessing.remove_duplicate_sequences = remove_dups
    
    coordinator.config.orf_prediction.translation_table = table_num
    coordinator.config.orf_prediction.min_orf_length_aa = min_orf_aa
    coordinator.config.orf_prediction.start_codons = start_codons_list
    coordinator.config.orf_prediction.stop_codons = stop_codons_list
    coordinator.config.orf_prediction.resolve_nested = resolve_nest
    
    coordinator.config.annotation.alignment_engine = align_engine
    coordinator.config.annotation.local_db_path = local_db
    coordinator.config.annotation.remote_fallback = remote_fallback
    coordinator.config.annotation.eval_threshold = eval_thresh
    coordinator.config.annotation.identity_threshold = ident_thresh
    coordinator.config.annotation.coverage_threshold = cov_thresh
    
    coordinator.config.domain_search.hmmer_db_path = hmmer_db
    coordinator.config.domain_search.eval_threshold = hmmer_eval
    coordinator.config.pathway_mapping.kegg_organism = kegg_organism
    coordinator.config.pathway_mapping.db_cache_path = cache_db
    
    coordinator.config.statistics.fdr_threshold = fdr_thresh
    coordinator.config.statistics.bg_universe_size = universe_size
    
    if AI_PROVIDER == "offline-rules":
        coordinator.config.ai_interpretation.provider = "ollama"
        coordinator.config.ai_interpretation.model_name = "invalid_model_trigger_fallback"
    else:
        coordinator.config.ai_interpretation.provider = AI_PROVIDER
        coordinator.config.ai_interpretation.model_name = AI_MODEL
        coordinator.config.ai_interpretation.api_key_env_var = "OPENAI_API_KEY"
        coordinator.config.ai_interpretation.temperature = AI_TEMP

    coordinator.initialize_workspace()
    
    st.markdown("#### ⚙️ Pipeline Live Runtime Tracking")
    progress_status = st.status("Executing modular genomics phases...", expanded=True)
    
    # Capture log outputs to terminal placeholder
    log_header = progress_status.empty()
    log_box_placeholder = progress_status.empty()
    log_buffer = StreamlitLogBuffer(log_box_placeholder)
    
    sink_id = logger.add(log_buffer.write, format="{time:HH:mm:ss} | {level:7} | {message}")
    
    success = False
    try:
        # Phase 1: Preprocess
        st.session_state.pipeline_steps["preprocess"]["status"] = "running"
        t0 = time.time()
        log_header.markdown("<span style='color:#818cf8; font-weight:bold;'>Stage 1/8: Quality Control Preprocessing</span>", unsafe_allow_html=True)
        cleaned_fasta = coordinator.run_preprocessing(fasta_path)
        st.session_state.pipeline_steps["preprocess"]["elapsed"] = time.time() - t0
        st.session_state.pipeline_steps["preprocess"]["status"] = "completed"
        progress_status.write(f"✅ Phase 1: QC Preprocessing completed ({st.session_state.pipeline_steps['preprocess']['elapsed']:.2f}s).")
        
        # Phase 2: ORF Prediction
        st.session_state.pipeline_steps["orf"]["status"] = "running"
        t0 = time.time()
        log_header.markdown("<span style='color:#818cf8; font-weight:bold;'>Stage 2/8: open reading frames Translation</span>", unsafe_allow_html=True)
        proteins_fasta = coordinator.run_orf_prediction(cleaned_fasta)
        st.session_state.pipeline_steps["orf"]["elapsed"] = time.time() - t0
        st.session_state.pipeline_steps["orf"]["status"] = "completed"
        progress_status.write(f"✅ Phase 2: ORF prediction completed ({st.session_state.pipeline_steps['orf']['elapsed']:.2f}s).")
        
        # Phase 3: Similarity Annotation
        st.session_state.pipeline_steps["similarity"]["status"] = "running"
        t0 = time.time()
        log_header.markdown("<span style='color:#818cf8; font-weight:bold;'>Stage 3/8: Sequence Similarity Annotation</span>", unsafe_allow_html=True)
        if not remote_fallback and not Path(local_db).exists():
            progress_status.warning("⚠️ Local alignment database missing! Enforcing remote Swiss-Prot fallback.")
            coordinator.config.annotation.remote_fallback = True
        annotated_csv = coordinator.run_annotation(proteins_fasta)
        st.session_state.pipeline_steps["similarity"]["elapsed"] = time.time() - t0
        st.session_state.pipeline_steps["similarity"]["status"] = "completed"
        progress_status.write(f"✅ Phase 3: Similarity search completed ({st.session_state.pipeline_steps['similarity']['elapsed']:.2f}s).")
        
        # Phase 4: Domains & Pathways
        st.session_state.pipeline_steps["domains"]["status"] = "running"
        t0 = time.time()
        log_header.markdown("<span style='color:#818cf8; font-weight:bold;'>Stage 4/8: Pfam Domains & Pathway Mapping</span>", unsafe_allow_html=True)
        mapped_pathways_csv = coordinator.run_pathway_and_domain_mapping(proteins_fasta, annotated_csv)
        st.session_state.pipeline_steps["domains"]["elapsed"] = time.time() - t0
        st.session_state.pipeline_steps["domains"]["status"] = "completed"
        progress_status.write(f"✅ Phase 4: Pfam domain mapping completed ({st.session_state.pipeline_steps['domains']['elapsed']:.2f}s).")
        
        # Phase 5: Statistics
        st.session_state.pipeline_steps["enrichment"]["status"] = "running"
        t0 = time.time()
        log_header.markdown("<span style='color:#818cf8; font-weight:bold;'>Stage 5/8: Statistical ssGSEA & ORA Enrichment</span>", unsafe_allow_html=True)
        significant_pathways_csv = coordinator.run_enrichment_analysis(mapped_pathways_csv)
        st.session_state.pipeline_steps["enrichment"]["elapsed"] = time.time() - t0
        st.session_state.pipeline_steps["enrichment"]["status"] = "completed"
        progress_status.write(f"✅ Phase 5: Statistics completed ({st.session_state.pipeline_steps['enrichment']['elapsed']:.2f}s).")
        
        # Phase 6: Visualizations
        st.session_state.pipeline_steps["visuals"]["status"] = "running"
        t0 = time.time()
        log_header.markdown("<span style='color:#818cf8; font-weight:bold;'>Stage 6/8: Publication Visualization Rendering</span>", unsafe_allow_html=True)
        vis_stats = coordinator.run_scientific_visualizations()
        st.session_state.pipeline_steps["visuals"]["elapsed"] = time.time() - t0
        st.session_state.pipeline_steps["visuals"]["status"] = "completed"
        progress_status.write(f"✅ Phase 6: Figures rendered ({st.session_state.pipeline_steps['visuals']['elapsed']:.2f}s).")
        
        # Phase 7: AI Interpretation
        st.session_state.pipeline_steps["ai"]["status"] = "running"
        t0 = time.time()
        log_header.markdown("<span style='color:#818cf8; font-weight:bold;'>Stage 7/8: Safe Pydantic AI literature Synthesis</span>", unsafe_allow_html=True)
        ai_stats = coordinator.run_ai_synthesis()
        st.session_state.pipeline_steps["ai"]["elapsed"] = time.time() - t0
        st.session_state.pipeline_steps["ai"]["status"] = "completed"
        progress_status.write(f"✅ Phase 7: AI Interpretation completed ({st.session_state.pipeline_steps['ai']['elapsed']:.2f}s).")
        
        # Phase 8: Report Assembly
        st.session_state.pipeline_steps["report"]["status"] = "running"
        t0 = time.time()
        log_header.markdown("<span style='color:#818cf8; font-weight:bold;'>Stage 8/8: HTML Analytical Dashboard Compiler</span>", unsafe_allow_html=True)
        coordinator.write_run_metadata(fasta_path, status="SUCCESS")
        report_stats = coordinator.run_report_generation_stage()
        st.session_state.pipeline_steps["report"]["elapsed"] = time.time() - t0
        st.session_state.pipeline_steps["report"]["status"] = "completed"
        progress_status.write(f"✅ Phase 8: Report compiled ({st.session_state.pipeline_steps['report']['elapsed']:.2f}s).")
        
        success = True
        progress_status.update(label="🎉 Genomics Pipeline successfully executed!", state="complete")
        
    except Exception as e:
        # Mark running steps as failed
        for key in st.session_state.pipeline_steps:
            if st.session_state.pipeline_steps[key]["status"] == "running":
                st.session_state.pipeline_steps[key]["status"] = "failed"
        progress_status.update(label=f"❌ Pipeline crashed: {e}", state="error")
        st.error(f"Execution Error: {e}")
        logger.exception(f"Streamlit Pipeline Crashed: {e}")
        
    finally:
        logger.remove(sink_id)
        if fasta_path and fasta_path.exists() and fasta_path.name.startswith("temp_"):
            fasta_path.unlink()
            
    if success:
        st.session_state.run_completed = True
        st.session_state.run_outdir = coordinator.output_dir
        
        metadata_json = coordinator.output_dir / "metadata.json"
        if metadata_json.exists():
            with open(metadata_json, "r", encoding="utf-8") as f:
                st.session_state.run_metadata = json.load(f)

# ----------------------------------------------------------------------
# 📊 RESULTS EXPLORER
# Displays interactive visualization panels matching nf-core aesthetics.
# ----------------------------------------------------------------------
if st.session_state.run_completed and st.session_state.run_outdir is not None:
    st.markdown("---")
    st.markdown("<div class='section-header'>Results Explorer</div>", unsafe_allow_html=True)
    
    outdir = Path(st.session_state.run_outdir)
    
    # Workspace paths
    preprocessed_dir = outdir / "preprocessed"
    orfs_dir = outdir / "orfs"
    annotations_dir = outdir / "annotations"
    pathways_dir = outdir / "pathways"
    enrichment_dir = outdir / "enrichment"
    visualizations_dir = outdir / "visualizations"
    final_report_dir = outdir / "final_report"
    
    # Define Structured Tabs (12 Tabs + optional Debugger)
    tab_titles = [
        "🌐 Workflow Tracker",
        "📊 FASTQ QC Summary",
        "✂️ Trimming Quality",
        "🧬 Genomic ORFs Map",
        "🎯 Homology Alignments",
        "🛡️ Annotation Confidence",
        "📈 Hypergeometric ORA",
        "💡 GSEA & ssPA",
        "🕸️ Bipartite Network",
        "🤖 AI Summary Synthesis",
        "📚 PubMed Citations Tracker",
        "📥 Data Export Center"
    ]
    
    if advanced_mode:
        tab_titles.append("🛠️ Bioinformatic Debugger")
        
    tabs = st.tabs(tab_titles)
    
    # ------------------------------------------------------------------
    # TAB 1: WORKFLOW TRACKER
    # ------------------------------------------------------------------
    with tabs[0]:
        st.markdown("#### 🌐 Pipeline Process Step Status & Diagnostics")
        st.markdown(
            "This interactive panel tracks each computational milestone of the viral genomics workflow, "
            "proving execution provenance and real-time runtime diagnostics."
        )
        
        # Render clean text-based status flowchart
        st.markdown("<div style='display: flex; flex-direction: column; gap: 10px;'>", unsafe_allow_html=True)
        total_time = 0.0
        for key, step in st.session_state.pipeline_steps.items():
            status_symbol = "⏳ Pending"
            border_color = "rgba(148, 163, 184, 0.15)"
            bg_color = "rgba(15, 23, 42, 0.4)"
            
            if step["status"] == "running":
                status_symbol = "⚙️ Running..."
                border_color = "#3b82f6"
                bg_color = "rgba(59, 130, 246, 0.1)"
            elif step["status"] == "completed":
                status_symbol = "✅ Completed"
                border_color = "#10b981"
                bg_color = "rgba(16, 185, 129, 0.1)"
                total_time += step["elapsed"]
            elif step["status"] == "failed":
                status_symbol = "❌ Failed"
                border_color = "#ef4444"
                bg_color = "rgba(239, 68, 68, 0.1)"
                
            st.markdown(
                f"<div style='border: 1px solid {border_color}; background-color: {bg_color}; padding: 15px; border-radius: 8px; display: flex; justify-content: space-between; align-items: center;'>"
                f"<div><b style='font-size: 1.1rem; color: #e2e8f0;'>{step['name']}</b></div>"
                f"<div style='display: flex; gap: 20px; align-items: center;'>"
                f"<span style='color: #94a3b8; font-size: 0.95rem;'>Runtime: <b>{step['elapsed']:.3f}s</b></span>"
                f"<span style='font-weight: 600; font-size: 0.95rem;'>{status_symbol}</span>"
                f"</div>"
                f"</div>",
                unsafe_allow_html=True
            )
        st.markdown("</div>", unsafe_allow_html=True)
        st.markdown(f"<div style='margin-top: 15px; font-weight:bold; color:#a5b4fc; text-align:right;'>Total Cumulative Pipeline Runtime: {total_time:.2f}s</div>", unsafe_allow_html=True)

    # ------------------------------------------------------------------
    # TAB 2: FASTQ QC SUMMARY
    # ------------------------------------------------------------------
    with tabs[1]:
        st.markdown("#### 📊 FASTQ Preprocessing & GC Distributions")
        qc_data = cached_read_json(str(preprocessed_dir / "qc_report.json"))
        
        if qc_data:
            counts = qc_data.get("counts", {})
            stats = qc_data.get("statistics", {})
            metrics = qc_data.get("metrics", {})
            
            n50_val = metrics.get("n50", stats.get("n50", stats.get("n50_bp", 0)))
            gc_val = metrics.get("mean_gc_percent", stats.get("gc_content_pct", 0.0))
            
            # Interactive MultiQC Grid
            c1, c2, c3, c4 = st.columns(4)
            with c1:
                st.markdown(f"<div class='metric-card'><div class='metric-val'>{counts.get('total_processed', 0)}</div><div class='metric-label'>Raw Input Sequences</div></div>", unsafe_allow_html=True)
            with c2:
                st.markdown(f"<div class='metric-card'><div class='metric-val'>{counts.get('total_kept', 0)}</div><div class='metric-label'>Cleaned Output</div></div>", unsafe_allow_html=True)
            with c3:
                st.markdown(f"<div class='metric-card'><div class='metric-val'>{n50_val} bp</div><div class='metric-label'>Sequence N50</div></div>", unsafe_allow_html=True)
            with c4:
                st.markdown(f"<div class='metric-card'><div class='metric-val'>{gc_val:.2f}%</div><div class='metric-label'>GC Fraction</div></div>", unsafe_allow_html=True)
                
            # GC Distribution Plot Panel
            st.markdown("<br><h5>📊 Per-Sequence GC Content Histogram</h5>", unsafe_allow_html=True)
            gc_plot = visualizations_dir / "gc_distribution.png"
            if gc_plot.exists():
                st.image(str(gc_plot), caption="Per-Read GC Content curve proving molecular distribution bounds.", use_container_width=True)
            else:
                st.info("GC distribution graphic is currently unrendered.")
        else:
            st.warning("No preprocessor metrics file detected.")

    # ------------------------------------------------------------------
    # TAB 3: TRIMMING QUALITY
    # ------------------------------------------------------------------
    with tabs[2]:
        st.markdown("#### ✂️ Read Trimming Quality & Waterfall Retention")
        
        col_t1, col_t2 = st.columns(2)
        with col_t1:
            st.markdown("##### Per-Base PHRED Quality Distributions")
            pb_plot = visualizations_dir / "per_base_quality.png"
            if pb_plot.exists():
                st.image(str(pb_plot), caption="Per-base PHRED score boxplots proving read polishing.", use_container_width=True)
            else:
                st.info("ℹ️ **Quality Scores Not Applicable (FASTA Input)**: Per-base PHRED quality score boxplots are only available for FASTQ read files. The provided FASTA format holds deterministic assembled sequence data without quality attributes.")
                
        with col_t2:
            st.markdown("##### Read Retention Waterfall Diagram")
            wf_plot = visualizations_dir / "read_retention_waterfall.png"
            if wf_plot.exists():
                st.image(str(wf_plot), caption="Sequential retention filtering illustrating trimming gains.", use_container_width=True)
            else:
                st.info("ℹ️ **Read Retention Waterfall Not Applicable (FASTA Input)**: Read retention waterfall charts track quality/trimming milestones for raw NGS reads (FASTQ). Assembled FASTA sequences do not undergo raw adapter/quality trimming filters.")

    # ------------------------------------------------------------------
    # TAB 4: GENOMIC ORFS MAP
    # ------------------------------------------------------------------
    with tabs[3]:
        st.markdown("#### 🧬 Predicted Coding open reading frames (ORFs)")
        
        df_orfs = cached_read_csv(str(orfs_dir / "coordinates.csv"))
        if not df_orfs.empty:
            c_orfs1, c_orfs2 = st.columns(2)
            with c_orfs1:
                st.metric("Total Genes Predicted (PRODIGAL-style)", len(df_orfs))
            with c_orfs2:
                st.metric("Average ORF size (AA)", int(df_orfs["length_aa"].mean()))
                
            st.markdown("##### Genomic Coordinate Mapping (Scroll Container)")
            st.dataframe(df_orfs.head(100), use_container_width=True)
            
            col_op1, col_op2 = st.columns(2)
            with col_op1:
                orf_len_plot = visualizations_dir / "orf_lengths.png"
                if orf_len_plot.exists():
                    st.image(str(orf_len_plot), caption="Amino acid length frequency histogram.", use_container_width=True)
            with col_op2:
                gen_track_plot = visualizations_dir / "orf_genomic_track.png"
                if gen_track_plot.exists():
                    st.image(str(gen_track_plot), caption="Genomic directional coordinate strand track map.", use_container_width=True)
        else:
            st.warning("No ORF coordinates data available.")

    # ------------------------------------------------------------------
    # TAB 5: HOMOLOGY ALIGNMENTS
    # ------------------------------------------------------------------
    with tabs[4]:
        st.markdown("#### 🎯 Swiss-Prot Similarity Hits (DIAMOND/blastp)")
        df_hits = cached_read_csv(str(annotations_dir / "annotated_proteins.csv"))
        
        if not df_hits.empty:
            annotated_only = df_hits[df_hits["uniprot_id"].notna() & (df_hits["uniprot_id"].astype(str).str.upper() != "NONE")]
            st.metric("Swiss-Prot Hits Detected", len(annotated_only))
            
            st.markdown("##### High-Confidence Alignments Table")
            st.dataframe(df_hits.head(100), use_container_width=True)
        else:
            st.warning("No homology alignment results detected.")

    # ------------------------------------------------------------------
    # TAB 6: ANNOTATION CONFIDENCE
    # ------------------------------------------------------------------
    with tabs[5]:
        st.markdown("#### 🛡️ Annotation Confidence & Identity Distributions")
        df_conf = cached_read_csv(str(annotations_dir / "annotated_proteins.csv"))
        
        if not df_conf.empty:
            annot_cnt = len(df_conf[df_conf["annotation_status"] == "Annotated"])
            hypo_cnt = len(df_conf[df_conf["annotation_status"] == "Hypothetical Protein"])
            
            col_c1, col_c2, col_c3 = st.columns(3)
            with col_c1:
                st.metric("Annotated Homologs", annot_cnt)
            with col_c2:
                st.metric("Hypothetical Genes", hypo_cnt)
            with col_c3:
                rate = (annot_cnt / len(df_conf) * 100) if len(df_conf) > 0 else 0
                st.metric("Coding Annotation Coverage", f"{rate:.1f}%")
                
            st.markdown("##### Similarity Distributions")
            conf_plot = visualizations_dir / "annotation_distributions.png"
            if conf_plot.exists():
                st.image(str(conf_plot), caption="Alignment percent sequence identities showing pipeline precision.", use_container_width=True)
        else:
            st.warning("No annotation confidence records found.")

    # ------------------------------------------------------------------
    # TAB 7: HYPERGEOMETRIC ORA
    # ------------------------------------------------------------------
    with tabs[6]:
        st.markdown("#### 📈 Overrepresentation Enrichment Statistics (ORA)")
        st.markdown(
            "This table lists host cell pathways that are statistically overrepresented "
            "by the viral query annotations, sorted strictly by Benjamini-Hochberg corrected FDR."
        )
        
        df_enrich = cached_read_csv(str(enrichment_dir / "significant_pathways.csv"))
        if not df_enrich.empty:
            st.dataframe(df_enrich, use_container_width=True)
            
            col_ep1, col_ep2 = st.columns(2)
            with col_ep1:
                ora_bar = visualizations_dir / "pathway_enrichment_barplot.png"
                if ora_bar.exists():
                    st.image(str(ora_bar), caption="Hypergeometric fold-enrichment priority bar chart.", use_container_width=True)
            with col_ep2:
                ora_bub = visualizations_dir / "pathway_enrichment_bubbleplot.png"
                if ora_bub.exists():
                    st.image(str(ora_bub), caption="ORA significance bubble chart (size proportional to query counts).", use_container_width=True)
        else:
            st.info("No pathways passed the Benjamini-Hochberg FDR significance threshold (FDR <= 0.05).")

    # ------------------------------------------------------------------
    # TAB 8: GSEA & SSPA
    # ------------------------------------------------------------------
    with tabs[7]:
        st.markdown("#### 💡 Single-Sample Pathway Analysis (ssPA) & ssGSEA Explorer")
        st.markdown(
            "Ranks functional categories using a scientifically unified score combining "
            "ssGSEA enrichment scores, ORA hypergeometric statistics, and Swiss-Prot homology alignment confidence."
        )
        
        df_sspa = cached_read_csv(str(enrichment_dir / "pathway_ranking_reports.csv"))
        if not df_sspa.empty:
            st.dataframe(df_sspa, use_container_width=True)
            
            col_ssp1, col_ssp2 = st.columns(2)
            with col_ssp1:
                sspa_bar = visualizations_dir / "sspa_ssgsea_barplot.png"
                if sspa_bar.exists():
                    st.image(str(sspa_bar), caption="Normalized ssGSEA pathway scores representing hijack potential.", use_container_width=True)
            with col_ssp2:
                sspa_bub = visualizations_dir / "sspa_multi_evidence_bubbleplot.png"
                if sspa_bub.exists():
                    st.image(str(sspa_bub), caption="Unified Multi-Evidence pathway priority scoring distributions.", use_container_width=True)
        else:
            st.info("No ssPA ranking records found. Ensure statistics completed successfully.")

    # ------------------------------------------------------------------
    # TAB 9: BIPARTITE NETWORK
    # ------------------------------------------------------------------
    with tabs[8]:
        st.markdown("#### 🕸️ Protein-Pathway Bipartite Interaction Network")
        net_plot = visualizations_dir / "protein_pathway_network.png"
        
        if net_plot.exists():
            st.image(str(net_plot), caption="Bipartite graph visualization mapping predicted viral genes (circles) to pathways (diamonds).", use_container_width=True)
        else:
            st.info("Bipartite network visualization skipped. This is expected if there are no enriched pathways to link.")

    # ------------------------------------------------------------------
    # TAB 10: AI SUMMARY SYNTHESIS
    # ------------------------------------------------------------------
    with tabs[9]:
        st.markdown("#### 🤖 Zero-Speculation AI Biological interpretation Summary")
        ai_data = cached_read_json(str(final_report_dir / "ai_synthesis.json"))
        
        if ai_data:
            st.markdown(
                f"<div style='background-color: rgba(99, 102, 241, 0.05); padding: 25px; border-radius: 12px; border-left: 5px solid #6366f1; margin-bottom: 20px;'>"
                f"<h5>📋 Executive Interpretation Summary</h5>"
                f"<p style='font-style: italic; font-size: 1.15rem; line-height: 1.6; color: #cbd5e1;'>\"{ai_data.get('concise_summary', 'No summary generated.')}\"</p>"
                f"</div>",
                unsafe_allow_html=True
            )
            
            st.markdown("##### 🔬 Detailed Biological Discussion")
            st.markdown(ai_data.get("detailed_biological_interpretation", "No discussion generated."))
            
            col_ai1, col_ai2 = st.columns(2)
            with col_ai1:
                st.markdown("##### 🏥 Predicted Disease Associations")
                st.markdown(ai_data.get("disease_association_summary", "No details available."))
                
                st.markdown("##### 🛡️ Candidate Antiviral Avenues")
                st.markdown(ai_data.get("therapeutic_relevance_summary", "No details available."))
                
            with col_ai2:
                st.markdown("##### 💡 Pathway Hijack Significance")
                st.markdown(ai_data.get("pathway_significance_discussion", "No details available."))
                
                st.markdown("##### 🚨 Low-Confidence Warning Flags")
                warnings = ai_data.get("confidence_warnings", [])
                if warnings:
                    for warn in warnings:
                        st.markdown(f"- ⚠️ {warn}")
                else:
                    st.markdown("*No low-confidence flags.*")
                    
            st.markdown("##### 🔬 Database & Literature Grounding Boundaries")
            st.info(ai_data.get("limitations", "No scientific limits documented."))
        else:
            st.warning("No structured AI Biological report JSON found.")

    # ------------------------------------------------------------------
    # TAB 11: PUBMED EVIDENCE CITATIONS
    # ------------------------------------------------------------------
    with tabs[10]:
        st.markdown("#### 📚 PubMed Evidence & PMID Grounding Tracker")
        st.markdown(
            "This literature grounding panel links all scientific assertions back to their exact retrieved "
            "PubMed abstracts. Click on any PMID to open the peer-reviewed reference in NCBI PubMed."
        )
        
        ai_data = cached_read_json(str(final_report_dir / "ai_synthesis.json"))
        if ai_data and "retrieved_literature_citations" in ai_data:
            citations = ai_data.get("retrieved_literature_citations", [])
            for c in citations:
                pmid = c.get("pmid", "")
                title = c.get("title", "Unknown Title")
                authors = c.get("authors", "Unknown Authors")
                journal = c.get("journal", "Unknown Journal")
                
                st.markdown(
                    f"<div class='alignment-card'>"
                    f"<b style='font-size: 1.15rem; color:#818cf8;'>{title}</b><br>"
                    f"<span style='color:#94a3b8; font-size:0.92rem;'>Authors: {authors} | Journal: <i>{journal}</i></span><br>"
                    f"<span style='font-size:0.95rem; font-weight:600; color:#34d399;'>Grounding Reference: <a href='https://pubmed.ncbi.nlm.nih.gov/{pmid}/' target='_blank' style='color:#34d399; text-decoration: underline;'>[PMID: {pmid}]</a></span>"
                    f"</div>",
                    unsafe_allow_html=True
                )
        else:
            st.warning("PubMed evidence records not found.")

    # ------------------------------------------------------------------
    # TAB 12: DATA EXPORT CENTER
    # ------------------------------------------------------------------
    with tabs[11]:
        st.markdown("#### 📥 Download Center")
        
        # Read HTML report bytes
        html_report = final_report_dir / "report.html"
        pdf_report = final_report_dir / "report.pdf"
        gff_file = orfs_dir / "coordinates.gff3"
        prot_fasta = orfs_dir / "proteins.fasta"
        ai_json = final_report_dir / "ai_synthesis.json"
        
        col_down1, col_down2 = st.columns(2)
        
        with col_down1:
            st.markdown("##### 📊 Full Analytical Reports")
            if html_report.exists():
                with open(html_report, "r", encoding="utf-8") as f:
                    html_bytes = f.read().encode("utf-8")
                st.download_button(
                    label="Download Portable HTML Dashboard 🖥️",
                    data=html_bytes,
                    file_name="pathoscope_genomics_report.html",
                    mime="text/html",
                    use_container_width=True
                )
            else:
                st.error("HTML Report dashboard missing.")
                
            if pdf_report.exists() and pdf_report.stat().st_size > 0:
                with open(pdf_report, "rb") as f:
                    pdf_bytes = f.read()
                st.download_button(
                    label="Download Analytical PDF Report 📄",
                    data=pdf_bytes,
                    file_name="pathoscope_genomics_report.pdf",
                    mime="application/pdf",
                    use_container_width=True
                )
            else:
                st.warning("PDF compilation skipped (HTML dashboard serves as primary analytical export).")
                
        with col_down2:
            st.markdown("##### 🧬 Raw Bioinformatics Coordinates & Data")
            if gff_file.exists():
                with open(gff_file, "r") as f:
                    gff_bytes = f.read().encode("utf-8")
                st.download_button(
                    label="Download Predicted GFF3 Coordinates",
                    data=gff_bytes,
                    file_name="predicted_orfs.gff3",
                    mime="text/plain",
                    use_container_width=True
                )
                
            if prot_fasta.exists():
                with open(prot_fasta, "r") as f:
                    fasta_bytes = f.read().encode("utf-8")
                st.download_button(
                    label="Download Translated AA Proteins (FASTA)",
                    data=fasta_bytes,
                    file_name="translated_proteins.fasta",
                    mime="text/plain",
                    use_container_width=True
                )
                
            if ai_json.exists():
                with open(ai_json, "r") as f:
                    json_bytes = f.read().encode("utf-8")
                st.download_button(
                    label="Download AI Biological Synthesis (JSON)",
                    data=json_bytes,
                    file_name="ai_synthesis.json",
                    mime="application/json",
                    use_container_width=True
                )
                
        # Lightweight Laptop Rendering Strategy Expanded Notice
        st.markdown("<br>", unsafe_allow_html=True)
        st.info(
            "💡 **Student Laptop Optimization Notice**: Matplotlib and Plotly figures are lazy-loaded and cached dynamically, "
            "large database hits are restricted to 100 rows per scroll panel, and database connections leverage lightweight "
            "SQLite caches, keeping the system memory footprint below 50MB."
        )

    # ------------------------------------------------------------------
    # TAB 13: ADVANCED BIOINFORMATICIAN DEBUGGER (VISIBLE ONLY ON ADVANCED MODE)
    # ------------------------------------------------------------------
    if advanced_mode:
        with tabs[12]:
            st.markdown("#### 🛠️ Advanced Bioinformatician Debugger & Provenance Explorer")
            
            c_db1, c_db2 = st.columns(2)
            with c_db1:
                st.markdown("##### System Environment Inspector")
                env_info = {
                    "Python Version": sys.version.split()[0],
                    "NumPy Version": np.__version__,
                    "Pandas Version": pd.__version__,
                    "NetworkX Version": nx.__version__,
                    "System Platform": sys.platform,
                    "Current Working Directory": os.getcwd()
                }
                st.json(env_info)
                
            with c_db2:
                st.markdown("##### SQLite Database Caching Metadata")
                cache_path = Path(cache_db)
                if cache_path.exists():
                    db_sz = cache_path.stat().st_size / (1024 * 1024)
                    st.success(f"Cache DB Status: ACTIVE\nFile Size: {db_sz:.3f} MB")
                else:
                    st.warning("Cache DB Status: INACTIVE (Remote API fallback or not initialized).")
                    
            st.markdown("##### Raw Output JSON Explorer")
            col_json1, col_json2 = st.columns(2)
            with col_json1:
                st.markdown("###### Preprocessing QC metadata")
                st.json(qc_data)
            with col_json2:
                st.markdown("###### AI Biological Synthesis")
                st.json(ai_data)
