# PathoScope AI v2.0 — Complete Foundation Package
## National Centre for Bioinformatics (NCB), Quaid-i-Azam University
### Functional Genomics Course — Sir Ghulam Abbas
**Team:** Muhammad Ekramah | Taha Javed | Asad Imam | Mishal Tariq

---

# ═══════════════════════════════════════════════════
# DOCUMENT 1 — PRODUCT REQUIREMENTS DOCUMENT (PRD)
# ═══════════════════════════════════════════════════

---

## 1. EXECUTIVE SUMMARY

PathoScope AI v2.0 is an automated, dual-workflow bioinformatics platform designed for academic use in functional genomics education and research. The system accepts multiple biological input types including FASTA sequences, FASTQ reads, gene lists (TXT/CSV), and expression matrices (CSV/TSV), and processes them through scientifically validated computational pipelines to produce biological insights, visualizations, and AI-assisted interpretations.

The platform addresses a critical educational gap: students and early-career researchers lack access to a unified tool that integrates quality control, sequence analysis, gene normalization, differential expression, pathway mapping, and biological interpretation in a single automated workflow.

PathoScope AI v2.0 replaces the previous architecture by introducing a scientifically corrected dual-workflow design that separates viral genomics analysis from functional genomics analysis, resolving the biological inconsistency of applying GSEA and ORA enrichment to raw viral FASTA files.

**Platform Type:** Command-line Python pipeline with web-based dashboard interface
**Primary Language:** Python 3.10+
**Target Environment:** Student laptops, university workstations, Linux/macOS/Windows
**Deployment Mode:** Local installation (offline-first, online-enhanced)

---

## 2. PRODUCT VISION

PathoScope AI aims to become the reference student-facing bioinformatics platform for functional genomics education at QAU, demonstrating that scientifically rigorous, publication-quality analysis can be achieved in an academic setting through automated, reproducible pipelines.

**Vision Statement:**
> "Give any biology student a raw genomic file or gene list and return a complete, scientifically defensible biological analysis report in minutes — with zero manual tool switching."

---

## 3. PROBLEM STATEMENT

### Problem 3.1 — Workflow Fragmentation
Current functional genomics analysis requires students to manually chain FASTQC, Trimmomatic, BLAST, KEGG REST API, Bioconductor R packages, and custom scripts. This fragmentation prevents reproducibility and consumes time that should be spent on biological interpretation.

### Problem 3.2 — Scientific Inconsistency in Previous Architecture
The previous PathoScope AI architecture applied GSEA and hypergeometric ORA enrichment directly to viral FASTA files. This is biologically invalid because:
- GSEA requires ranked gene lists derived from expression data (fold changes, t-statistics).
- ORA requires a defined background gene universe, which small viral genomes cannot provide.
- Viral genomes contain 4–100 proteins; enrichment statistics become meaningless at this scale.

### Problem 3.3 — Missing Functional Genomics Support
Professor requirements explicitly include gene list processing, expression matrices, Gene ID normalization (HGNC, Ensembl, Entrez), differential expression analysis (logFC, p-value, FDR), and classification of upregulated/downregulated genes. The previous architecture lacked this entire workflow branch.

### Problem 3.4 — Missing Gene ID Normalization
Biological datasets use heterogeneous identifier systems. A pipeline that cannot normalize between HGNC symbols, Ensembl IDs, and Entrez IDs will fail on real-world data.

### Problem 3.5 — Reproducibility Deficit
No configuration system, no environment specification, no run tracking, and no local database caching existed in the previous design, making demonstration and viva defense unreliable.

---

## 4. TARGET USERS

### Primary Users
- Bioinformatics MSc/BSc students at NCB, QAU
- Functional genomics course participants
- Computational biology researchers performing preliminary analysis

### Secondary Users
- Research supervisors reviewing student work
- Virology researchers needing rapid annotation
- Course instructors demonstrating analysis workflows

### User Personas

**Persona A — Beginner Student (Amna)**
- No command-line experience
- Understands basic molecular biology
- Needs a guided web interface to upload files and receive results
- Cannot interpret raw BLAST output
- Needs AI-plain-language explanations

**Persona B — Advanced Student (Bilal)**
- Comfortable with Python
- Can modify config files
- Wants to adjust thresholds and rerun analysis
- Interested in pathway-level interpretation
- Will present results in viva

**Persona C — Research Supervisor (Sir Ghulam Abbas)**
- Evaluates scientific validity
- Checks reproducibility
- Expects correct statistical methodology
- Will ask about biological justification of design decisions

---

## 5. FUNCTIONAL REQUIREMENTS

### FR-01 — Input Handling
- System SHALL accept FASTA files (.fa, .fasta, .fna)
- System SHALL accept FASTQ files (.fq, .fastq, .fastq.gz)
- System SHALL accept gene lists (.txt, .csv) containing HGNC symbols, Ensembl IDs, or Entrez IDs
- System SHALL accept expression matrices (.csv, .tsv) with gene identifiers and numerical expression values
- System SHALL auto-detect file format from extension and content
- System SHALL reject malformed files with descriptive error messages

### FR-02 — Workflow Routing
- System SHALL detect input type and route to Workflow A (Viral Genomics) or Workflow B (Functional Genomics)
- System SHALL allow manual workflow selection via config.yaml or CLI flag
- System SHALL log which workflow was selected and why

### FR-03 — Workflow A: Viral Genomics
- System SHALL validate FASTA/FASTQ sequences for format correctness
- System SHALL remove duplicate sequences
- System SHALL filter sequences with ambiguity ratio above configurable threshold (default: N > 10%)
- System SHALL perform six-frame ORF prediction
- System SHALL translate predicted ORFs to protein sequences
- System SHALL perform functional annotation via local BLAST or API fallback
- System SHALL map proteins to Pfam domains
- System SHALL map proteins to KEGG Orthology entries
- System SHALL map proteins to Reactome pathways (where coverage allows)
- System SHALL generate functional categorization tables (not statistical enrichment on <50 proteins)
- System SHALL generate visualizations: ORF map, domain frequency chart, KEGG category bar chart
- System SHALL send validated results to Claude API for biological interpretation
- System SHALL generate PDF/HTML report

### FR-04 — Workflow B: Functional Genomics
- System SHALL normalize gene identifiers to HGNC symbols as primary format
- System SHALL convert Ensembl IDs to HGNC via local mapping table + Ensembl REST API
- System SHALL convert Entrez IDs to HGNC via NCBI E-utilities
- System SHALL handle mixed identifier lists (auto-detect per row)
- System SHALL apply statistical filtering: logFC threshold, p-value cutoff, FDR correction
- System SHALL classify genes as: Upregulated (logFC > 1, padj < 0.05), Downregulated (logFC < -1, padj < 0.05), Not Significant
- System SHALL perform ORA enrichment against KEGG pathways using valid background gene universe
- System SHALL perform ORA enrichment against GO Biological Process terms
- System SHALL generate volcano plots, enrichment bar charts, GO bubble plots
- System SHALL generate DEG classification tables
- System SHALL send validated results to Claude API for biological interpretation
- System SHALL generate PDF/HTML report

### FR-05 — Statistical Analysis
- System SHALL apply Benjamini-Hochberg FDR correction to all multiple hypothesis tests
- System SHALL implement Fisher's exact test for ORA
- System SHALL use configurable thresholds (config.yaml) for all cutoffs
- System SHALL NEVER fabricate statistical values

### FR-06 — AI Interpretation
- System SHALL use Claude API (claude-sonnet-4-20250514) for biological interpretation
- AI SHALL only interpret validated computational outputs
- AI SHALL NOT generate statistics, p-values, or pathway assignments
- System SHALL include a disclaimer on every AI interpretation output

### FR-07 — Reporting
- System SHALL generate HTML report for every run
- System SHALL generate CSV tables for annotation, DEG, and pathway results
- System SHALL generate PDF report (optional, wkhtmltopdf/WeasyPrint)
- System SHALL generate JSON metadata file for every run
- System SHALL save all outputs to timestamped run directory

### FR-08 — Offline/Online Mode
- System SHALL function in offline mode using locally cached databases
- System SHALL use online APIs when internet is available and cache responses in SQLite
- System SHALL log whether each data source was online or offline
- System SHALL never fail silently; if API fails, system uses local fallback and logs the event

### FR-09 — Reproducibility
- System SHALL read all parameters from config.yaml
- System SHALL write run metadata (input file, parameters, tool versions, timestamp) to metadata.json
- System SHALL support re-running any previous analysis by providing its run ID
- System SHALL log all steps to run.log

---

## 6. NON-FUNCTIONAL REQUIREMENTS

- NFR-01: Complete pipeline must finish within 15 minutes on a 4-core student laptop for standard-sized inputs
- NFR-02: System must work on Python 3.10+ on Linux, macOS, and Windows (WSL)
- NFR-03: All third-party dependencies must be specified in requirements.txt and environment.yml
- NFR-04: System must operate without internet using local fallback databases
- NFR-05: All database versions must be recorded in metadata.json
- NFR-06: Code must be modular (one function per biological task, max 50 lines per function)
- NFR-07: All modules must include unit tests

---

## 7. MVP SCOPE (Version 1.0)

**In Scope:**
- Workflow A: FASTA input, QC, ORF prediction, protein translation, Pfam annotation, KEGG mapping, functional categorization, visualization, AI interpretation, HTML report
- Workflow B: Gene list input (HGNC/Ensembl/Entrez), Gene ID normalization, statistical filtering, DEG classification, KEGG ORA, GO ORA, volcano plot, enrichment chart, AI interpretation, HTML report
- Offline mode with SQLite caching
- config.yaml system
- Logging and metadata tracking

**Out of Scope for v1.0:**
- GSEA (requires large expression datasets not typical in course setting)
- Multi-sample RNA-seq normalization (DESeq2/edgeR methods — requires R integration)
- Custom database uploads
- User authentication system
- Cloud deployment
- Real-time collaborative analysis

---

## 8. SUCCESS METRICS

| Metric | Target |
|---|---|
| FASTA file processed end-to-end | < 10 minutes for SARS-CoV-2 genome |
| Gene list normalized (500 genes) | < 2 minutes |
| Offline mode functional | Yes — all core features work without internet |
| Report generated per run | Yes — HTML + CSV always generated |
| Tests passing | > 90% unit test pass rate |
| Reproducibility | Same input produces identical outputs |
| Viva demonstration | Live demo works on student laptop |

---

## 9. SCIENTIFIC OBJECTIVES

1. Demonstrate biologically valid ORF prediction using six-frame translation with configurable minimum ORF length
2. Perform functional annotation using Pfam domain mapping (scientifically appropriate for viral proteins)
3. Apply KEGG Orthology categorization (not enrichment) for small viral proteomes
4. Apply statistically valid ORA enrichment only to functional genomics datasets with appropriate background universes
5. Correctly implement Benjamini-Hochberg FDR correction for multiple hypothesis testing
6. Correctly classify DEGs using validated logFC and adjusted p-value thresholds

---

## 10. ACADEMIC OBJECTIVES

1. Demonstrate mastery of Python-based bioinformatics pipeline development
2. Show understanding of functional genomics data types and their appropriate analyses
3. Demonstrate ability to integrate biological databases (Pfam, KEGG, GO, Ensembl)
4. Produce reproducible, documented scientific software
5. Prepare for viva defense with scientifically defensible design decisions

---

## 11. RISKS AND LIMITATIONS

| Risk | Likelihood | Mitigation |
|---|---|---|
| KEGG API rate limiting | High | Local SQLite cache, rate limiter |
| Ensembl API downtime | Medium | Local mapping table fallback |
| BLAST local installation | Medium | Pre-built local database included |
| Small viral proteome enrichment weakness | High | Replace enrichment with functional categorization |
| Student hardware performance | Medium | Efficient Python, no unnecessary loops |

---

## 12. FEATURES EXCLUDED FROM VERSION 1

- GSEA analysis (biologically inappropriate for available input types without large expression datasets)
- Smith-Waterman implementation from scratch (educational value does not justify engineering cost; use BioPython pairwise2 instead)
- D3.js visualizations (engineering complexity not justified; use Plotly/Matplotlib instead)
- Multi-omics integration
- Cloud API deployment
- Real-time streaming results

---
---

# ═══════════════════════════════════════════════════
# DOCUMENT 2 — TECHNICAL REQUIREMENTS DOCUMENT (TRD)
# ═══════════════════════════════════════════════════

---

## 1. SYSTEM ARCHITECTURE OVERVIEW

PathoScope AI v2.0 follows a modular pipeline architecture organized into three tiers:

```
┌──────────────────────────────────────────────────────────────────┐
│                         TIER 1 — INPUT LAYER                      │
│    InputManager → FileDetector → Validator → WorkflowRouter       │
└──────────────────────────┬───────────────────────────────────────┘
                           │
           ┌───────────────┴───────────────┐
           ▼                               ▼
┌──────────────────────┐       ┌──────────────────────────────┐
│  TIER 2A             │       │  TIER 2B                     │
│  WORKFLOW A          │       │  WORKFLOW B                  │
│  Viral Genomics      │       │  Functional Genomics         │
│                      │       │                              │
│  QC Engine           │       │  GeneIDNormalizer            │
│  ORFPredictor        │       │  StatisticalFilter           │
│  ProteinTranslator   │       │  DEGClassifier               │
│  PfamMapper          │       │  ORAEnrichment               │
│  KEGGMapper          │       │  GOEnrichment                │
│  FunctionalCategorizer│       │  PathwayMapper               │
└──────────────────────┘       └──────────────────────────────┘
           │                               │
           └───────────────┬───────────────┘
                           ▼
┌──────────────────────────────────────────────────────────────────┐
│                      TIER 3 — OUTPUT LAYER                        │
│  VisualizationEngine → AIInterpreter → ReportGenerator           │
│  MetadataWriter → Logger → DatabaseCache                          │
└──────────────────────────────────────────────────────────────────┘
```

---

## 2. COMPLETE FOLDER STRUCTURE

```
pathoscope_ai/
│
├── main.py                          # Entry point — CLI
├── config.yaml                      # All configurable parameters
├── requirements.txt                 # Python dependencies
├── environment.yml                  # Conda environment
├── README.md                        # Installation and usage guide
├── .env.example                     # API key template (never commit .env)
│
├── core/                            # Core pipeline modules
│   ├── __init__.py
│   ├── input_manager.py             # File detection and routing
│   ├── file_detector.py             # Format auto-detection
│   ├── validator.py                 # Input validation
│   ├── workflow_router.py           # Routes to Workflow A or B
│   └── config_loader.py             # Loads and validates config.yaml
│
├── workflow_a/                      # Viral Genomics Pipeline
│   ├── __init__.py
│   ├── qc_engine.py                 # Quality control for FASTA/FASTQ
│   ├── sequence_cleaner.py          # Normalization, duplicate removal
│   ├── orf_predictor.py             # Six-frame ORF prediction
│   ├── protein_translator.py        # ORF → protein sequences
│   ├── annotation_engine.py         # BLAST-based annotation
│   ├── pfam_mapper.py               # Pfam domain mapping
│   ├── kegg_mapper_a.py             # KEGG Orthology mapping for proteins
│   └── functional_categorizer.py   # Non-statistical categorization
│
├── workflow_b/                      # Functional Genomics Pipeline
│   ├── __init__.py
│   ├── gene_id_normalizer.py        # HGNC/Ensembl/Entrez normalization
│   ├── expression_loader.py         # CSV/TSV expression matrix loading
│   ├── statistical_filter.py        # logFC, p-value, FDR filtering
│   ├── deg_classifier.py            # DEG classification logic
│   ├── ora_engine.py                # ORA Fisher's exact test
│   ├── go_enrichment.py             # GO term enrichment
│   └── kegg_mapper_b.py             # KEGG pathway enrichment for genes
│
├── databases/                       # Local database files
│   ├── pfam_local/                  # Pfam HMM profiles (subset)
│   ├── kegg_cache.sqlite            # KEGG API response cache
│   ├── ensembl_hgnc_map.tsv         # Ensembl → HGNC mapping table
│   ├── entrez_hgnc_map.tsv          # Entrez → HGNC mapping table
│   ├── go_annotations.sqlite        # GO term annotations cache
│   └── gene_universe.tsv            # Background gene universe for ORA
│
├── output/                          # All run outputs (auto-created)
│   └── run_YYYYMMDD_HHMMSS/
│       ├── metadata.json
│       ├── run.log
│       ├── qc_report.html
│       ├── annotation_results.csv
│       ├── deg_results.csv
│       ├── pathway_results.csv
│       ├── visualizations/
│       │   ├── orf_map.png
│       │   ├── volcano_plot.png
│       │   ├── enrichment_chart.png
│       │   └── domain_frequency.png
│       ├── ai_interpretation.txt
│       └── final_report.html
│
├── visualization/                   # Visualization engine
│   ├── __init__.py
│   ├── orf_visualizer.py
│   ├── volcano_plotter.py
│   ├── enrichment_plotter.py
│   ├── domain_plotter.py
│   └── report_builder.py
│
├── ai/                              # AI interpretation module
│   ├── __init__.py
│   ├── claude_interpreter.py        # Claude API integration
│   └── prompt_templates.py          # Structured prompts for Claude
│
├── utils/                           # Shared utilities
│   ├── __init__.py
│   ├── logger.py                    # Logging system
│   ├── metadata_writer.py           # metadata.json generation
│   ├── db_cache.py                  # SQLite caching layer
│   ├── api_client.py                # HTTP client with retry logic
│   └── exceptions.py                # Custom exception classes
│
├── tests/                           # Test suite
│   ├── __init__.py
│   ├── test_data/                   # Test datasets
│   │   ├── fasta/
│   │   │   ├── ms2_bacteriophage.fasta
│   │   │   ├── phix174.fasta
│   │   │   ├── sars_cov2.fasta
│   │   │   ├── influenza_a.fasta
│   │   │   ├── invalid_chars.fasta
│   │   │   ├── high_ambiguity.fasta
│   │   │   └── duplicates.fasta
│   │   ├── fastq/
│   │   │   ├── valid_reads.fastq
│   │   │   └── corrupted.fastq
│   │   └── gene_lists/
│   │       ├── hgnc_symbols.txt
│   │       ├── ensembl_ids.txt
│   │       ├── entrez_ids.txt
│   │       ├── mixed_ids.csv
│   │       └── expression_matrix.csv
│   ├── test_qc_engine.py
│   ├── test_orf_predictor.py
│   ├── test_gene_normalizer.py
│   ├── test_statistical_filter.py
│   ├── test_deg_classifier.py
│   ├── test_ora_engine.py
│   └── test_report_generator.py
│
└── docs/                            # Documentation
    ├── PRD.md
    ├── TRD.md
    ├── viva_preparation.md
    └── scientific_justifications.md
```

---

## 3. TECHNOLOGY STACK

### Core Language
- Python 3.10+ (mandatory)

### Bioinformatics Libraries
| Library | Version | Purpose |
|---|---|---|
| BioPython | >=1.81 | FASTA/FASTQ parsing, sequence manipulation, ORF detection |
| pyhmmer | >=0.10 | Pfam HMM domain search (replaces HMMER binary dependency) |
| requests | >=2.31 | API calls (KEGG, Ensembl, NCBI) |
| pandas | >=2.0 | Data handling, expression matrices, results tables |
| numpy | >=1.24 | Numerical operations |
| scipy | >=1.11 | Fisher's exact test, statistical functions |
| statsmodels | >=0.14 | Benjamini-Hochberg FDR correction |

### Visualization Libraries
| Library | Version | Purpose |
|---|---|---|
| matplotlib | >=3.7 | ORF maps, domain frequency charts |
| plotly | >=5.15 | Interactive volcano plots, enrichment charts |
| seaborn | >=0.12 | Heatmaps, additional charts |

### Database and Caching
| Library | Version | Purpose |
|---|---|---|
| sqlite3 | stdlib | Local caching database |
| PyYAML | >=6.0 | config.yaml parsing |
| python-dotenv | >=1.0 | .env file loading for API keys |

### Reporting
| Library | Version | Purpose |
|---|---|---|
| jinja2 | >=3.1 | HTML report templating |
| WeasyPrint | >=60.0 | PDF generation from HTML |

### AI Integration
| Library | Purpose |
|---|---|
| anthropic | Official Claude API Python client |

### Testing
| Library | Purpose |
|---|---|
| pytest | Test runner |
| pytest-cov | Coverage reporting |

---

## 4. MODULE SPECIFICATIONS

### Module 4.1 — InputManager (core/input_manager.py)

**Responsibility:** Accept file path from CLI or web interface, detect format, validate existence, route to appropriate handler.

**Functions:**
```
load_input(file_path: str) -> InputObject
detect_format(file_path: str) -> str  # returns: "fasta", "fastq", "csv", "txt"
validate_file_exists(file_path: str) -> bool
route_to_workflow(input_obj: InputObject) -> str  # returns: "A" or "B"
```

**Input:** File path string
**Output:** InputObject with format, content, and workflow designation

### Module 4.2 — QCEngine (workflow_a/qc_engine.py)

**Responsibility:** Validate and filter FASTA/FASTQ sequences.

**FASTA QC Functions:**
```
validate_fasta_format(sequences: list) -> QCReport
check_invalid_characters(seq: str) -> bool
calculate_ambiguity_ratio(seq: str) -> float
filter_by_length(sequences: list, min_len: int) -> list
remove_duplicates(sequences: list) -> list
generate_qc_summary(sequences: list) -> dict
```

**FASTQ QC Functions:**
```
validate_fastq_format(records: list) -> QCReport
calculate_mean_quality(quality_scores: list) -> float
filter_by_quality(records: list, min_q: int) -> list
trim_low_quality_ends(record, min_q: int) -> record
```

**Thresholds (configurable via config.yaml):**
- min_sequence_length: 50 bp
- max_ambiguity_ratio: 0.10
- min_fastq_quality: 20

### Module 4.3 — ORFPredictor (workflow_a/orf_predictor.py)

**Responsibility:** Identify all Open Reading Frames across six reading frames.

**Algorithm:**
1. Translate sequence in 6 frames (3 forward, 3 reverse complement)
2. Find all ATG start codons
3. Find next in-frame stop codon (TAA, TAG, TGA)
4. Extract ORF if length >= min_orf_length (default: 100 bp / 33 aa)
5. Score ORFs by length
6. Return ranked ORF list

**Functions:**
```
predict_orfs(sequence: str, min_length: int) -> list[ORF]
translate_six_frames(sequence: str) -> dict
find_start_codons(frame_seq: str) -> list[int]
find_stop_codons(frame_seq: str, start: int) -> int
extract_orf(sequence: str, start: int, stop: int, frame: int) -> ORF
rank_orfs(orfs: list) -> list[ORF]
```

**Output:** List of ORF objects with: sequence, start, stop, frame, length, strand

### Module 4.4 — GeneIDNormalizer (workflow_b/gene_id_normalizer.py)

**Responsibility:** Convert any gene identifier format to HGNC symbol.

**Supported Input Formats:**
- HGNC gene symbols (e.g., TP53, BRCA1)
- Ensembl gene IDs (e.g., ENSG00000141510)
- Entrez gene IDs (e.g., 7157)
- Mixed (auto-detected per row)

**Conversion Workflow:**
```
normalize_gene_list(genes: list) -> NormalizationResult
detect_id_type(gene_id: str) -> str  # "hgnc", "ensembl", "entrez", "unknown"
ensembl_to_hgnc(ensembl_id: str) -> str  # local table first, API fallback
entrez_to_hgnc(entrez_id: str) -> str    # local table first, API fallback
validate_hgnc_symbol(symbol: str) -> bool
handle_unmapped(gene_id: str) -> str     # logs and returns None
```

**Data Sources (priority order):**
1. Local TSV mapping tables (databases/ensembl_hgnc_map.tsv, entrez_hgnc_map.tsv)
2. Ensembl REST API (https://rest.ensembl.org/xrefs/id/)
3. NCBI E-utilities (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/)
4. MyGene.info API (https://mygene.info/v3/gene/)

**Output:** NormalizationResult with: original_id, normalized_symbol, conversion_method, success_flag

### Module 4.5 — StatisticalFilter (workflow_b/statistical_filter.py)

**Responsibility:** Apply statistical thresholds to expression data.

**Functions:**
```
load_expression_matrix(file_path: str) -> pd.DataFrame
validate_matrix_columns(df: pd.DataFrame) -> bool
apply_logfc_threshold(df: pd.DataFrame, threshold: float) -> pd.DataFrame
apply_pvalue_filter(df: pd.DataFrame, cutoff: float) -> pd.DataFrame
apply_fdr_correction(df: pd.DataFrame) -> pd.DataFrame  # Benjamini-Hochberg
get_significant_genes(df: pd.DataFrame) -> pd.DataFrame
```

**Required Columns in Expression Matrix:**
- gene_id (any format, normalized by GeneIDNormalizer)
- log2FoldChange (or logFC)
- pvalue
- padj (or computed by pipeline via BH correction)

**Thresholds (config.yaml):**
```yaml
statistical:
  logfc_up_threshold: 1.0        # logFC > 1.0 → Upregulated
  logfc_down_threshold: -1.0     # logFC < -1.0 → Downregulated
  pvalue_cutoff: 0.05
  fdr_cutoff: 0.05
```

### Module 4.6 — DEGClassifier (workflow_b/deg_classifier.py)

**Responsibility:** Classify genes into upregulated, downregulated, not significant.

**Functions:**
```
classify_genes(df: pd.DataFrame) -> pd.DataFrame
get_upregulated(df: pd.DataFrame) -> pd.DataFrame
get_downregulated(df: pd.DataFrame) -> pd.DataFrame
get_not_significant(df: pd.DataFrame) -> pd.DataFrame
generate_deg_summary(df: pd.DataFrame) -> dict
```

**Classification Logic:**
```
Upregulated:   logFC > 1.0  AND padj < 0.05
Downregulated: logFC < -1.0 AND padj < 0.05
Not Significant: all others
```

### Module 4.7 — ORAEngine (workflow_b/ora_engine.py)

**Responsibility:** Perform Over-Representation Analysis using Fisher's exact test.

**Algorithm:**
1. Input: list of DEGs, background gene universe
2. For each pathway/term: construct 2x2 contingency table
3. Apply Fisher's exact test (scipy.stats.fisher_exact)
4. Collect all p-values
5. Apply Benjamini-Hochberg FDR correction (statsmodels)
6. Filter by FDR < 0.05
7. Return significant pathways ranked by FDR

**Background Universe:** databases/gene_universe.tsv (Human Gene Universe, ~20,000 genes)

**Functions:**
```
run_ora(deg_list: list, gene_sets: dict, universe: list) -> ORAResult
build_contingency_table(deg_list, pathway_genes, universe) -> np.array
fisher_exact_test(table: np.array) -> float
apply_bh_correction(pvalues: list) -> list
filter_significant(results: pd.DataFrame) -> pd.DataFrame
```

---

## 5. API INTEGRATIONS

### KEGG REST API
- Base URL: https://rest.kegg.jp
- Endpoints: /list/pathway, /get/{pathway_id}, /link/pathway/{gene_id}
- Rate limit: 3 requests/second
- Caching: SQLite (kegg_cache.sqlite), TTL: 7 days

### Ensembl REST API
- Base URL: https://rest.ensembl.org
- Endpoints: /xrefs/id/{ensembl_id}
- Rate limit: 15 requests/second
- Caching: SQLite, TTL: 30 days (gene IDs are stable)

### NCBI E-utilities
- Base URL: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
- Endpoints: efetch.fcgi, esearch.fcgi
- Rate limit: 3 requests/second (10 with API key)
- Caching: SQLite, TTL: 30 days

### Claude API (Anthropic)
- Model: claude-sonnet-4-20250514
- Max tokens: 2000 per interpretation request
- Input: Structured JSON of validated biological results
- Output: Plain language biological interpretation paragraph
- No caching (interpretations may vary by context)

---

## 6. DATABASE DESIGN (SQLite)

### Table: kegg_cache
```sql
CREATE TABLE kegg_cache (
    cache_key TEXT PRIMARY KEY,
    response_data TEXT NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    expires_at TIMESTAMP NOT NULL
);
```

### Table: gene_id_cache
```sql
CREATE TABLE gene_id_cache (
    original_id TEXT PRIMARY KEY,
    id_type TEXT NOT NULL,
    hgnc_symbol TEXT,
    conversion_method TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

### Table: go_annotations
```sql
CREATE TABLE go_annotations (
    go_id TEXT NOT NULL,
    go_term TEXT NOT NULL,
    go_category TEXT NOT NULL,
    gene_symbol TEXT NOT NULL,
    evidence_code TEXT,
    PRIMARY KEY (go_id, gene_symbol)
);
```

---

## 7. CONFIGURATION SYSTEM (config.yaml)

```yaml
# PathoScope AI v2.0 Configuration
project:
  name: "PathoScope AI"
  version: "2.0.0"
  output_dir: "./output"
  log_level: "INFO"

workflow:
  mode: "auto"  # auto | workflow_a | workflow_b

quality_control:
  min_sequence_length: 50
  max_ambiguity_ratio: 0.10
  min_fastq_quality: 20
  remove_duplicates: true

orf_prediction:
  min_orf_length: 100
  genetic_code: 1  # Standard genetic code
  both_strands: true

statistical:
  logfc_up_threshold: 1.0
  logfc_down_threshold: -1.0
  pvalue_cutoff: 0.05
  fdr_cutoff: 0.05
  fdr_method: "benjamini-hochberg"

annotation:
  blast_evalue: 0.001
  blast_identity: 30.0
  blast_coverage: 50.0

api:
  kegg_base_url: "https://rest.kegg.jp"
  ensembl_base_url: "https://rest.ensembl.org"
  ncbi_base_url: "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  request_timeout: 30
  retry_attempts: 3
  offline_mode: false

cache:
  sqlite_path: "./databases/cache.sqlite"
  kegg_ttl_days: 7
  gene_id_ttl_days: 30

ai:
  model: "claude-sonnet-4-20250514"
  max_tokens: 2000
  enabled: true
```

---

## 8. ERROR HANDLING FRAMEWORK

All exceptions are handled through custom exception classes in utils/exceptions.py:

```
PathoScopeError (base)
├── InputValidationError
├── FileFormatError
├── SequenceQualityError
├── GeneNormalizationError
├── APIConnectionError
├── DatabaseError
├── StatisticalError
└── ReportGenerationError
```

Every module catches specific exceptions, logs with full traceback, and either:
- Falls back to local data source (API errors)
- Skips invalid records and logs count (sequence errors)
- Raises to main pipeline with informative message (critical errors)

---

## 9. LOGGING STRATEGY

Every run creates a timestamped run directory with run.log:

```
[2026-05-31 10:23:45] [INFO]    PathoScope AI v2.0 started
[2026-05-31 10:23:45] [INFO]    Input file: sars_cov2.fasta
[2026-05-31 10:23:45] [INFO]    Detected format: FASTA
[2026-05-31 10:23:45] [INFO]    Routed to Workflow A (Viral Genomics)
[2026-05-31 10:23:46] [INFO]    QC: 1 sequence loaded, 0 removed
[2026-05-31 10:23:46] [INFO]    ORF Prediction: 12 ORFs found
[2026-05-31 10:23:50] [WARNING] KEGG API unavailable, using local cache
[2026-05-31 10:23:51] [INFO]    KEGG Mapping: 8/12 proteins mapped
[2026-05-31 10:24:10] [INFO]    Report generated: output/run_20260531_102345/
[2026-05-31 10:24:10] [INFO]    Total runtime: 25.3 seconds
```

---

## 10. REPRODUCIBILITY FRAMEWORK

Every run generates metadata.json:

```json
{
  "run_id": "run_20260531_102345",
  "timestamp": "2026-05-31T10:23:45Z",
  "input_file": "sars_cov2.fasta",
  "input_md5": "a3f2e1...",
  "workflow": "A",
  "config_snapshot": { ... },
  "software_versions": {
    "pathoscope_ai": "2.0.0",
    "python": "3.10.14",
    "biopython": "1.81",
    "pandas": "2.0.3",
    "scipy": "1.11.0"
  },
  "database_versions": {
    "pfam": "36.0",
    "kegg_cache_date": "2026-05-28",
    "ensembl_release": "110"
  },
  "parameters_used": { ... },
  "outputs_generated": [...]
}
```

---

*End of Document 2 — TRD*
