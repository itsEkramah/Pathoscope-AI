# PathoScope AI v2.0
# DOCUMENT 5 — BACKEND SCHEMA DOCUMENT
# DOCUMENT 6 — IMPLEMENTATION PLAN

---

# ═══════════════════════════════════════════════════
# DOCUMENT 5 — BACKEND SCHEMA DOCUMENT
# ═══════════════════════════════════════════════════

---

## 1. DATABASE ARCHITECTURE OVERVIEW

PathoScope AI uses a hybrid storage strategy:
- **SQLite** for all local caching, run tracking, and offline databases
- **File system** for large data outputs (FASTA sequences, reports, visualizations)
- **JSON** for run metadata and configuration snapshots
- **TSV/CSV** for static biological mapping tables (Ensembl→HGNC, Entrez→HGNC)

All SQLite databases live in the `databases/` directory.
All run outputs live in `output/run_YYYYMMDD_HHMMSS/` directories.

---

## 2. SQLITE DATABASE SCHEMAS

### Database File: `databases/pathoscope.sqlite`
Primary operational database for run tracking and caching.

---

### TABLE: runs

Tracks every pipeline run executed.

```sql
CREATE TABLE runs (
    run_id          TEXT PRIMARY KEY,
    -- Format: run_YYYYMMDD_HHMMSS
    -- Example: run_20260531_102345

    input_file      TEXT NOT NULL,
    -- Original filename provided by user

    input_md5       TEXT NOT NULL,
    -- MD5 hash of input file for reproducibility verification

    input_format    TEXT NOT NULL,
    -- Values: 'fasta', 'fastq', 'csv', 'tsv', 'txt'

    workflow        TEXT NOT NULL,
    -- Values: 'A', 'B'

    status          TEXT NOT NULL DEFAULT 'running',
    -- Values: 'running', 'complete', 'failed', 'cancelled'

    started_at      TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    completed_at    TIMESTAMP,
    runtime_seconds REAL,

    output_dir      TEXT NOT NULL,
    -- Absolute path to run output directory

    config_snapshot TEXT NOT NULL,
    -- JSON string of complete config.yaml used for this run

    error_message   TEXT,
    -- NULL if successful, error text if failed

    python_version  TEXT NOT NULL,
    pathoscope_version TEXT NOT NULL DEFAULT '2.0.0'
);

CREATE INDEX idx_runs_started_at ON runs(started_at DESC);
CREATE INDEX idx_runs_status ON runs(status);
```

---

### TABLE: qc_results

Stores quality control metrics for each run.

```sql
CREATE TABLE qc_results (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,

    -- Input counts
    sequences_input     INTEGER NOT NULL,
    sequences_passed    INTEGER NOT NULL,
    sequences_removed   INTEGER NOT NULL,

    -- Removal reasons
    removed_ambiguity   INTEGER NOT NULL DEFAULT 0,
    removed_duplicates  INTEGER NOT NULL DEFAULT 0,
    removed_length      INTEGER NOT NULL DEFAULT 0,
    removed_invalid     INTEGER NOT NULL DEFAULT 0,

    -- Sequence statistics
    avg_length          REAL,
    min_length          INTEGER,
    max_length          INTEGER,
    total_bases         INTEGER,
    gc_content_percent  REAL,
    avg_ambiguity_ratio REAL,

    -- FASTQ-specific (NULL for FASTA runs)
    avg_quality_score   REAL,
    reads_below_quality INTEGER,

    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

---

### TABLE: orf_results

Stores all predicted ORFs for Workflow A runs.

```sql
CREATE TABLE orf_results (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,

    orf_id          TEXT NOT NULL,
    -- Format: orf_001, orf_002...

    source_sequence TEXT NOT NULL,
    -- Sequence ID from input FASTA

    start_position  INTEGER NOT NULL,
    stop_position   INTEGER NOT NULL,
    length_bp       INTEGER NOT NULL,
    length_aa       INTEGER NOT NULL,
    reading_frame   INTEGER NOT NULL,
    -- Values: 1, 2, 3, -1, -2, -3

    strand          TEXT NOT NULL,
    -- Values: 'forward', 'reverse'

    dna_sequence    TEXT NOT NULL,
    protein_sequence TEXT NOT NULL,

    start_codon     TEXT NOT NULL DEFAULT 'ATG',
    stop_codon      TEXT,

    UNIQUE (run_id, orf_id)
);

CREATE INDEX idx_orf_run_id ON orf_results(run_id);
CREATE INDEX idx_orf_length ON orf_results(length_bp DESC);
```

---

### TABLE: annotation_results

Stores protein functional annotation results.

```sql
CREATE TABLE annotation_results (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,
    orf_id          TEXT NOT NULL,
    -- References orf_results.orf_id

    annotation_method TEXT NOT NULL,
    -- Values: 'blast_local', 'blast_api', 'uniprot_api', 'none'

    hit_accession   TEXT,
    hit_description TEXT,
    hit_organism    TEXT,
    hit_length      INTEGER,

    identity_percent REAL,
    coverage_percent REAL,
    evalue          REAL,
    bitscore        REAL,

    -- KEGG mapping (from annotation hit)
    kegg_orthology  TEXT,
    -- KO number, e.g., K02164

    kegg_pathways   TEXT,
    -- JSON array of KEGG pathway IDs

    -- Pfam mapping
    pfam_domains    TEXT,
    -- JSON array: [{"domain_id": "PF00552", "domain_name": "...", "evalue": ...}]

    annotation_status TEXT NOT NULL DEFAULT 'annotated',
    -- Values: 'annotated', 'partial', 'unannotated'

    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    UNIQUE (run_id, orf_id)
);

CREATE INDEX idx_annotation_run_id ON annotation_results(run_id);
```

---

### TABLE: gene_normalization_results

Stores gene ID normalization results for Workflow B.

```sql
CREATE TABLE gene_normalization_results (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,

    original_id     TEXT NOT NULL,
    detected_type   TEXT NOT NULL,
    -- Values: 'hgnc', 'ensembl', 'entrez', 'unknown'

    hgnc_symbol     TEXT,
    -- NULL if normalization failed

    conversion_method TEXT,
    -- Values: 'local_table', 'ensembl_api', 'ncbi_api', 'mygene_api'

    conversion_success BOOLEAN NOT NULL DEFAULT FALSE,
    error_reason    TEXT,
    -- NULL if successful

    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_gene_norm_run_id ON gene_normalization_results(run_id);
CREATE INDEX idx_gene_norm_success ON gene_normalization_results(conversion_success);
```

---

### TABLE: deg_results

Stores differential expression analysis results for Workflow B.

```sql
CREATE TABLE deg_results (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,

    gene_symbol     TEXT NOT NULL,
    -- HGNC symbol (normalized)

    original_id     TEXT NOT NULL,
    -- Original ID from input file

    log2_fold_change REAL NOT NULL,
    pvalue          REAL NOT NULL,
    padj            REAL NOT NULL,
    -- Benjamini-Hochberg corrected

    classification  TEXT NOT NULL,
    -- Values: 'upregulated', 'downregulated', 'not_significant'

    -- Thresholds applied (from config, for reproducibility)
    logfc_threshold_applied REAL NOT NULL,
    pvalue_threshold_applied REAL NOT NULL,
    fdr_threshold_applied   REAL NOT NULL,

    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    UNIQUE (run_id, gene_symbol)
);

CREATE INDEX idx_deg_run_id ON deg_results(run_id);
CREATE INDEX idx_deg_classification ON deg_results(classification);
CREATE INDEX idx_deg_padj ON deg_results(padj);
```

---

### TABLE: pathway_results

Stores ORA enrichment and pathway mapping results.

```sql
CREATE TABLE pathway_results (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,

    pathway_id      TEXT NOT NULL,
    -- Example: hsa04110 (KEGG), GO:0006915 (GO)

    pathway_name    TEXT NOT NULL,
    pathway_database TEXT NOT NULL,
    -- Values: 'KEGG', 'GO_BP', 'GO_MF', 'GO_CC', 'Reactome'

    analysis_type   TEXT NOT NULL,
    -- Values: 'ORA', 'functional_categorization'

    -- ORA statistics
    gene_count_in_pathway  INTEGER NOT NULL,
    -- Number of DEGs in this pathway

    gene_count_in_universe INTEGER NOT NULL,
    -- Total genes in pathway (background)

    total_degs      INTEGER NOT NULL,
    -- Total DEGs tested

    universe_size   INTEGER NOT NULL,
    -- Background universe size

    pvalue          REAL NOT NULL,
    padj            REAL NOT NULL,
    -- BH-corrected

    gene_ratio      REAL NOT NULL,
    -- gene_count_in_pathway / total_degs

    background_ratio REAL NOT NULL,
    -- gene_count_in_universe / universe_size

    gene_list       TEXT NOT NULL,
    -- JSON array of gene symbols in this pathway

    significant     BOOLEAN NOT NULL DEFAULT FALSE,
    -- TRUE if padj < 0.05

    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    UNIQUE (run_id, pathway_id, pathway_database)
);

CREATE INDEX idx_pathway_run_id ON pathway_results(run_id);
CREATE INDEX idx_pathway_padj ON pathway_results(padj);
CREATE INDEX idx_pathway_significant ON pathway_results(significant);
```

---

### TABLE: visualization_metadata

Tracks all generated visualization files.

```sql
CREATE TABLE visualization_metadata (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,

    chart_type      TEXT NOT NULL,
    -- Values: 'orf_map', 'volcano_plot', 'enrichment_bar',
    --         'go_bubble', 'domain_frequency', 'kegg_bar'

    file_path       TEXT NOT NULL,
    -- Relative path within output directory

    file_format     TEXT NOT NULL,
    -- Values: 'png', 'svg', 'html'

    chart_title     TEXT,
    data_points     INTEGER,
    -- Number of features plotted

    generation_library TEXT NOT NULL,
    -- Values: 'matplotlib', 'plotly', 'seaborn'

    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

---

### TABLE: ai_interpretations

Stores AI-generated biological interpretations.

```sql
CREATE TABLE ai_interpretations (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,

    model_used      TEXT NOT NULL DEFAULT 'claude-sonnet-4-20250514',
    interpretation_type TEXT NOT NULL,
    -- Values: 'workflow_a_summary', 'workflow_b_summary',
    --         'pathway_interpretation', 'deg_interpretation'

    input_prompt    TEXT NOT NULL,
    -- Structured prompt sent to Claude

    output_text     TEXT NOT NULL,
    -- Claude's response

    tokens_used     INTEGER,
    generation_time_seconds REAL,

    disclaimer_text TEXT NOT NULL DEFAULT
        'This interpretation was generated by an AI language model (Claude).
         All statistical values and pathway assignments were derived from
         validated computational analysis. AI provides contextual
         biological explanation only and should not be treated as
         expert clinical or research opinion.',

    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

---

### TABLE: kegg_cache

Caches KEGG API responses to support offline mode.

```sql
CREATE TABLE kegg_cache (
    cache_key       TEXT PRIMARY KEY,
    -- Format: "kegg_{endpoint}_{parameter}"
    -- Example: "kegg_pathway_hsa04110"

    endpoint        TEXT NOT NULL,
    parameter       TEXT NOT NULL,
    response_data   TEXT NOT NULL,
    -- JSON string of API response

    http_status     INTEGER NOT NULL,
    response_size   INTEGER NOT NULL,

    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    expires_at      TIMESTAMP NOT NULL,

    access_count    INTEGER NOT NULL DEFAULT 1,
    last_accessed   TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_kegg_expires ON kegg_cache(expires_at);
```

---

### TABLE: gene_id_cache

Caches gene ID conversion results.

```sql
CREATE TABLE gene_id_cache (
    original_id     TEXT PRIMARY KEY,
    id_type         TEXT NOT NULL,
    hgnc_symbol     TEXT,
    conversion_method TEXT,
    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    expires_at      TIMESTAMP NOT NULL
    -- Long TTL: 30 days (gene IDs are stable)
);
```

---

### TABLE: go_term_cache

Caches GO term annotation data.

```sql
CREATE TABLE go_term_cache (
    go_id           TEXT NOT NULL,
    go_term         TEXT NOT NULL,
    go_category     TEXT NOT NULL,
    -- Values: 'biological_process', 'molecular_function', 'cellular_component'

    gene_symbol     TEXT NOT NULL,
    evidence_code   TEXT,
    taxon_id        INTEGER NOT NULL DEFAULT 9606,
    -- 9606 = Homo sapiens

    PRIMARY KEY (go_id, gene_symbol)
);

CREATE INDEX idx_go_symbol ON go_term_cache(gene_symbol);
CREATE INDEX idx_go_category ON go_term_cache(go_category);
```

---

## 3. FILE SYSTEM SCHEMA

### Run Output Directory Structure:
```
output/
└── run_20260531_102345/
    │
    ├── metadata.json              # Complete run metadata
    ├── run.log                    # Full execution log
    ├── config_snapshot.yaml       # Config used for this run
    │
    ├── qc/
    │   ├── qc_report.html
    │   ├── qc_metrics.csv
    │   └── qc_summary.json
    │
    ├── sequences/                 # Workflow A only
    │   ├── cleaned_sequences.fasta
    │   ├── predicted_orfs.fasta
    │   └── translated_proteins.fasta
    │
    ├── annotation/                # Workflow A only
    │   ├── annotation_results.csv
    │   ├── pfam_results.csv
    │   └── kegg_mapping.csv
    │
    ├── deg/                       # Workflow B only
    │   ├── normalized_genes.csv
    │   ├── deg_results.csv
    │   ├── upregulated.csv
    │   └── downregulated.csv
    │
    ├── enrichment/                # Workflow B only
    │   ├── kegg_ora_results.csv
    │   └── go_ora_results.csv
    │
    ├── visualizations/
    │   ├── orf_map.png
    │   ├── orf_map.svg
    │   ├── volcano_plot.png
    │   ├── volcano_plot.html      # Interactive Plotly
    │   ├── enrichment_chart.png
    │   ├── enrichment_chart.html
    │   ├── domain_frequency.png
    │   └── go_bubble_chart.html
    │
    ├── ai/
    │   └── ai_interpretation.txt
    │
    └── reports/
        ├── final_report.html
        └── final_report.pdf       # If WeasyPrint available
```

---

## 4. STATIC DATABASE FILES

### databases/ensembl_hgnc_map.tsv
```
# Format: Ensembl_ID \t HGNC_Symbol \t Entrez_ID \t Gene_Name
# Source: Ensembl Release 110, Homo sapiens
# Date: 2026-01-15
ENSG00000141510    TP53    7157    Tumor Protein P53
ENSG00000012048    BRCA1   672     BRCA1 DNA Repair Associated
...
```

### databases/entrez_hgnc_map.tsv
```
# Format: Entrez_ID \t HGNC_Symbol \t Ensembl_ID \t Gene_Name
7157    TP53    ENSG00000141510    Tumor Protein P53
672     BRCA1   ENSG00000012048    BRCA1 DNA Repair Associated
...
```

### databases/gene_universe.tsv
```
# Human protein-coding gene universe for ORA background
# Source: HGNC Complete Dataset, 2026-01-15
# Total genes: 19,194
TP53
BRCA1
BRCA2
...
```

---

## 5. ENTITY RELATIONSHIP DIAGRAM (Text Format)

```
runs (1) ─────────────── (1) qc_results
  │
  ├──────────────────── (many) orf_results
  │                              │
  │                              └── (1) annotation_results
  │
  ├──────────────────── (many) gene_normalization_results
  │
  ├──────────────────── (many) deg_results
  │
  ├──────────────────── (many) pathway_results
  │
  ├──────────────────── (many) visualization_metadata
  │
  └──────────────────── (many) ai_interpretations

kegg_cache ─── (standalone, shared across runs)
gene_id_cache ─── (standalone, shared across runs)
go_term_cache ─── (standalone, shared across runs)
```

---

## 6. DATA RETENTION POLICIES

| Table | Retention | Notes |
|---|---|---|
| runs | Forever | Core audit trail |
| qc_results | Forever | Linked to runs |
| orf_results | 90 days | Can be regenerated from input |
| annotation_results | 90 days | Can be regenerated |
| gene_normalization_results | 90 days | Cheap to regenerate |
| deg_results | Forever | Scientific outputs |
| pathway_results | Forever | Scientific outputs |
| ai_interpretations | 30 days | Non-deterministic |
| kegg_cache | TTL: 7 days | Auto-purged on startup |
| gene_id_cache | TTL: 30 days | Auto-purged on startup |

---
---

# ═══════════════════════════════════════════════════
# DOCUMENT 6 — IMPLEMENTATION PLAN
# ═══════════════════════════════════════════════════

---

## IMPLEMENTATION PHILOSOPHY

Build order follows **biological dependency chain** — each phase depends on the previous being validated before the next begins. Never skip validation steps. Every phase ends with a working, testable deliverable.

**Golden Rule:** If it doesn't work on MS2 Bacteriophage or a 10-gene list first, don't proceed to SARS-CoV-2 or 500-gene lists.

---

## PHASE 0 — ENVIRONMENT AND ARCHITECTURE SETUP
**Duration:** 1 day
**Goal:** Working project skeleton that runs without errors

### 0.1 Objectives
- Set up Python virtual environment
- Install all dependencies
- Create complete folder structure
- Create config.yaml with all defaults
- Set up logging system
- Set up SQLite database with all schemas
- Verify all imports work

### 0.2 Files to Create
```
pathoscope_ai/
├── main.py                     (CLI entry point)
├── config.yaml                 (all parameters)
├── requirements.txt            (all dependencies)
├── environment.yml             (conda spec)
├── .env.example                (API key template)
├── core/__init__.py
├── core/config_loader.py
├── utils/__init__.py
├── utils/logger.py
├── utils/exceptions.py
├── utils/db_cache.py
├── utils/metadata_writer.py
└── databases/schema.sql        (all CREATE TABLE statements)
```

### 0.3 Key Functions (Phase 0)
```python
# utils/logger.py
setup_logger(run_id: str, output_dir: str) -> logging.Logger
log_step(logger, step_name: str, status: str, details: dict)

# core/config_loader.py
load_config(config_path: str) -> dict
validate_config(config: dict) -> bool
get_threshold(config: dict, key: str) -> any

# utils/db_cache.py
initialize_database(db_path: str) -> sqlite3.Connection
create_all_tables(conn: sqlite3.Connection)
purge_expired_cache(conn: sqlite3.Connection)
```

### 0.4 Validation Criteria (Phase 0)
- [ ] `python main.py --help` runs without error
- [ ] config.yaml loads successfully
- [ ] SQLite database created with all tables
- [ ] Logger creates run directory and log file
- [ ] All imports succeed: biopython, pandas, scipy, statsmodels, plotly, matplotlib, requests, anthropic

---

## PHASE 1 — INPUT HANDLING AND WORKFLOW ROUTING
**Duration:** 2 days
**Goal:** Any supported file is correctly identified and validated

### 1.1 Objectives
- Detect file format from extension + content inspection
- Validate file existence and non-empty
- Route to Workflow A or B correctly
- Reject malformed inputs with clear errors

### 1.2 Files to Create
```
core/input_manager.py
core/file_detector.py
core/validator.py
core/workflow_router.py
```

### 1.3 Key Functions (Phase 1)
```python
# core/file_detector.py
detect_format(file_path: str) -> str
    # Returns: 'fasta', 'fastq', 'csv', 'tsv', 'txt_genes', 'txt_expression'
    # Method: check extension, then read first 5 lines to confirm format
    
is_fasta(file_path: str) -> bool
    # Check: starts with '>', no quality lines

is_fastq(file_path: str) -> bool
    # Check: lines 1,5,9... start with '@', line 3 starts with '+'

is_expression_matrix(file_path: str) -> bool
    # Check: CSV/TSV with headers containing 'logFC' or 'log2FoldChange' and 'pvalue'

is_gene_list(file_path: str) -> bool
    # Check: TXT or CSV, first column contains gene identifiers

# core/workflow_router.py
determine_workflow(detected_format: str) -> str
    # 'fasta' or 'fastq' → 'A'
    # 'csv', 'tsv', 'txt_genes', 'txt_expression' → 'B'
```

### 1.4 Test Data Required
- tests/test_data/fasta/ms2_bacteriophage.fasta (valid)
- tests/test_data/fasta/invalid_chars.fasta (contains X, Z, J)
- tests/test_data/fasta/high_ambiguity.fasta (>10% N)
- tests/test_data/fastq/valid_reads.fastq
- tests/test_data/gene_lists/hgnc_symbols.txt
- tests/test_data/gene_lists/expression_matrix.csv

### 1.5 Unit Tests (Phase 1)
```python
# tests/test_input_manager.py
test_detect_fasta_format()
test_detect_fastq_format()
test_detect_csv_expression()
test_detect_txt_gene_list()
test_reject_unsupported_format()
test_route_fasta_to_workflow_a()
test_route_gene_list_to_workflow_b()
test_route_expression_to_workflow_b()
test_reject_empty_file()
test_reject_nonexistent_file()
```

### 1.6 Validation Criteria (Phase 1)
- [ ] MS2 FASTA file → detected as 'fasta' → routed to Workflow A
- [ ] gene_symbols.txt → detected as 'txt_genes' → routed to Workflow B
- [ ] expression_matrix.csv → detected as 'csv_expression' → routed to Workflow B
- [ ] invalid_chars.fasta → detected as 'fasta' → ValidationError raised with list of invalid chars
- [ ] All 10 unit tests pass

---

## PHASE 2 — QUALITY CONTROL ENGINE (Workflow A)
**Duration:** 2 days
**Goal:** FASTA/FASTQ sequences are validated and cleaned

### 2.1 Objectives
- Validate FASTA format (BioPython SeqIO)
- Remove duplicate sequences
- Filter sequences with invalid characters
- Filter sequences with excessive ambiguity (N > threshold)
- Filter sequences below minimum length
- Generate QC report
- Handle FASTQ quality scores

### 2.2 Files to Create
```
workflow_a/qc_engine.py
workflow_a/sequence_cleaner.py
```

### 2.3 Key Functions (Phase 2)
```python
# workflow_a/qc_engine.py
run_fasta_qc(file_path: str, config: dict) -> QCResult
    # Returns: QCResult(passed_sequences, removed, qc_metrics)

validate_fasta_format(sequences: list) -> list[str]
    # Returns list of validation errors (empty if valid)

check_invalid_characters(seq: str) -> list[str]
    # Returns list of invalid chars found
    # Valid chars: ACGTRYSWKMBDHVN (IUPAC nucleotide codes)
    # Invalid: anything else

calculate_ambiguity_ratio(seq: str) -> float
    # Returns N count / total length

filter_by_length(sequences, min_len: int) -> tuple[list, int]
    # Returns (passed_sequences, removed_count)

remove_duplicates(sequences) -> tuple[list, int]
    # Uses sequence MD5 for deduplication
    # Returns (unique_sequences, removed_count)

# workflow_a/qc_engine.py — FASTQ specific
run_fastq_qc(file_path: str, config: dict) -> QCResult
calculate_mean_phred_quality(quality_string: str) -> float
filter_by_quality(records, min_q: int) -> tuple[list, int]
```

### 2.4 Algorithm: FASTA QC Pipeline
```
Step 1: Parse all sequences with BioPython SeqIO.parse()
Step 2: For each sequence:
        a. Convert to uppercase
        b. Check for invalid characters → REJECT if found
        c. Calculate ambiguity ratio → REJECT if N > max_ambiguity_ratio
        d. Check length → REJECT if < min_sequence_length
        e. Calculate MD5 of sequence → REJECT if duplicate
Step 3: Generate QC report with all metrics
Step 4: Write cleaned sequences to output FASTA
Step 5: Log: X sequences in, Y passed, Z removed (reasons)
```

### 2.5 Test Datasets
- MS2 Bacteriophage (NC_001417): 1 sequence, 3,569 bp — should pass QC cleanly
- PhiX174 (NC_001422): 1 sequence, 5,386 bp — should pass QC cleanly
- invalid_chars.fasta: ATGCXYZ... — should fail character check
- high_ambiguity.fasta: ATGNNNNNNNNN... (>10% N) — should fail ambiguity check
- duplicates.fasta: Same sequence repeated 3x — should be collapsed to 1

### 2.6 Validation Criteria (Phase 2)
- [ ] MS2 FASTA passes QC with 0 sequences removed
- [ ] invalid_chars.fasta: all sequences rejected, error message lists invalid chars
- [ ] high_ambiguity.fasta: sequences above threshold removed, count logged
- [ ] duplicates.fasta: 3 sequences in → 1 sequence out
- [ ] QC metrics CSV generated with correct values

---

## PHASE 3 — ORF PREDICTION AND PROTEIN TRANSLATION (Workflow A)
**Duration:** 3 days
**Goal:** All valid ORFs predicted and translated to protein sequences

### 3.1 Objectives
- Implement six-frame ORF prediction
- Filter ORFs by minimum length
- Translate to protein sequences
- Write output FASTA files
- Visualize ORF positions on genome

### 3.2 Files to Create
```
workflow_a/orf_predictor.py
workflow_a/protein_translator.py
```

### 3.3 Algorithm: Six-Frame ORF Prediction
```
Input: DNA sequence (uppercase, validated)

Step 1: Generate 6 reading frames:
        Frame +1: sequence[0::1]
        Frame +2: sequence[1::1]
        Frame +3: sequence[2::1]
        Frame -1: reverse_complement(sequence)[0::1]
        Frame -2: reverse_complement(sequence)[1::1]
        Frame -3: reverse_complement(sequence)[2::1]

Step 2: For each frame:
        a. Scan for ATG (start codon)
        b. From ATG, scan for next in-frame stop (TAA, TAG, TGA)
        c. Extract sequence between start and stop
        d. Calculate ORF length
        e. If length >= min_orf_length: record ORF

Step 3: Sort all ORFs by length (descending)
Step 4: Remove ORFs fully contained within longer ORFs (optional, configurable)
Step 5: Assign ORF IDs: orf_001, orf_002...
Step 6: Write orf_predictions.fasta (DNA sequences)
Step 7: Write predicted_proteins.fasta (AA sequences)
```

### 3.4 Key Functions (Phase 3)
```python
# workflow_a/orf_predictor.py
predict_orfs(sequence: str, config: dict) -> list[ORF]
get_reverse_complement(sequence: str) -> str
translate_six_frames(sequence: str) -> dict[str, str]
find_orfs_in_frame(frame_seq: str, frame_id: str, offset: int) -> list[ORF]
filter_orfs_by_length(orfs: list, min_len: int) -> list[ORF]
remove_nested_orfs(orfs: list) -> list[ORF]
assign_orf_ids(orfs: list) -> list[ORF]

# workflow_a/protein_translator.py
translate_orf(dna_sequence: str, genetic_code: int) -> str
translate_all_orfs(orfs: list, genetic_code: int) -> list[ORF]
write_protein_fasta(orfs: list, output_path: str)
```

### 3.5 BioPython Usage
```python
from Bio.Seq import Seq
from Bio import SeqIO

# Use BioPython's genetic code tables
seq = Seq(dna_string)
protein = seq.translate(table=genetic_code, to_stop=True)
reverse_comp = seq.reverse_complement()
```

### 3.6 Validation Against Known Genomes
| Genome | Known ORF Count | Expected Range |
|---|---|---|
| MS2 Bacteriophage (3,569 bp) | 4 canonical ORFs | 4–10 (including small ORFs) |
| PhiX174 (5,386 bp) | 11 canonical ORFs | 10–20 |
| SARS-CoV-2 (29,903 bp) | ~12 canonical ORFs | 12–30 |

### 3.7 Validation Criteria (Phase 3)
- [ ] MS2: minimum 4 ORFs predicted in forward frames
- [ ] SARS-CoV-2: Spike protein ORF detected (>3,000 bp)
- [ ] Protein sequences start with M (Methionine) for all ORFs
- [ ] DNA ORF FASTA and protein FASTA written correctly
- [ ] ORF count logged with frame distribution

---

## PHASE 4 — GENE ID NORMALIZATION (Workflow B)
**Duration:** 3 days
**Goal:** Any gene identifier format converted to HGNC symbol

### 4.1 Objectives
- Detect identifier type (HGNC/Ensembl/Entrez/unknown)
- Convert to HGNC using priority: local table → API
- Handle unmapped genes gracefully
- Generate normalization report

### 4.2 Files to Create
```
workflow_b/gene_id_normalizer.py
databases/ensembl_hgnc_map.tsv    (download from Ensembl BioMart)
databases/entrez_hgnc_map.tsv     (download from NCBI)
```

### 4.3 Download Local Mapping Tables
```bash
# Ensembl → HGNC mapping (run once, save to databases/)
# Source: Ensembl BioMart or FTP
wget https://ftp.ensembl.org/pub/current_tsv/homo_sapiens/
# Or use MyGene.info bulk download

# Entrez → HGNC mapping
wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
```

### 4.4 Algorithm: Gene ID Detection
```python
ENSEMBL_PATTERN = r'^ENSG\d{11}$'
ENTREZ_PATTERN  = r'^\d+$'
HGNC_PATTERN    = r'^[A-Z][A-Z0-9\-]+$'  # simplified

def detect_id_type(gene_id: str) -> str:
    if re.match(ENSEMBL_PATTERN, gene_id):
        return 'ensembl'
    elif re.match(ENTREZ_PATTERN, gene_id):
        return 'entrez'
    elif re.match(HGNC_PATTERN, gene_id):
        return 'hgnc'
    else:
        return 'unknown'
```

### 4.5 Conversion Priority
```
1. Check SQLite gene_id_cache (fastest)
2. Check local TSV mapping table (fast, offline)
3. Query Ensembl REST API (online, cache result)
4. Query NCBI E-utilities (online, cache result)
5. Query MyGene.info (online, last resort, cache result)
6. Log as 'unmapped' (never silently drop genes)
```

### 4.6 Validation Criteria (Phase 4)
- [ ] ENSG00000141510 → TP53 (via local table)
- [ ] 7157 → TP53 (via local table)
- [ ] TP53 → TP53 (identity, recognized as HGNC)
- [ ] Mixed list: 100 genes normalized, unmapped count logged
- [ ] Offline mode: all conversions use local table, no API calls

---

## PHASE 5 — STATISTICAL FILTERING AND DEG CLASSIFICATION (Workflow B)
**Duration:** 2 days
**Goal:** Expression data filtered and genes correctly classified

### 5.1 Objectives
- Load expression matrix (CSV/TSV)
- Validate required columns
- Apply BH FDR correction if padj not present
- Apply logFC and padj thresholds
- Classify as upregulated/downregulated/not significant
- Generate DEG summary statistics

### 5.2 Algorithm: Statistical Filtering
```python
import pandas as pd
from statsmodels.stats.multitest import multipletests

def run_statistical_analysis(df: pd.DataFrame, config: dict) -> pd.DataFrame:
    
    # Step 1: Validate columns
    required = ['gene_id', 'log2FoldChange', 'pvalue']
    
    # Step 2: Apply BH correction if padj not present
    if 'padj' not in df.columns:
        reject, padj, _, _ = multipletests(df['pvalue'], method='fdr_bh')
        df['padj'] = padj
    
    # Step 3: Classify
    logfc_up = config['statistical']['logfc_up_threshold']      # default 1.0
    logfc_dn = config['statistical']['logfc_down_threshold']    # default -1.0
    fdr_cut  = config['statistical']['fdr_cutoff']             # default 0.05
    
    def classify(row):
        if row['log2FoldChange'] > logfc_up and row['padj'] < fdr_cut:
            return 'upregulated'
        elif row['log2FoldChange'] < logfc_dn and row['padj'] < fdr_cut:
            return 'downregulated'
        else:
            return 'not_significant'
    
    df['classification'] = df.apply(classify, axis=1)
    return df
```

### 5.3 Validation Criteria (Phase 5)
- [ ] TP53 with logFC=2.5, padj=0.001 → classified as upregulated
- [ ] BRCA1 with logFC=-2.1, padj=0.003 → classified as downregulated
- [ ] GAPDH with logFC=0.3, padj=0.4 → classified as not_significant
- [ ] FDR correction applied when padj column absent
- [ ] DEG summary counts match manual calculation

---

## PHASE 6 — FUNCTIONAL ANNOTATION AND DOMAIN MAPPING (Workflow A)
**Duration:** 3 days
**Goal:** Predicted proteins annotated with Pfam domains and KEGG Orthology

### 6.1 Pfam Domain Mapping

**Tool:** pyhmmer (Python-native HMMER, no binary dependency)

```python
import pyhmmer

def search_pfam(protein_sequences: list, pfam_database_path: str) -> list[PfamResult]:
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(pfam_database_path) as hmm_file:
        hmms = list(hmm_file)
    
    sequences = [
        pyhmmer.easel.TextSequence(name=s.id.encode(), sequence=str(s.seq))
        for s in protein_sequences
    ]
    
    results = []
    for hits in pyhmmer.hmmsearch(hmms, sequences):
        for hit in hits:
            if hit.evalue < evalue_threshold:
                results.append(PfamResult(...))
    return results
```

**Note:** Download Pfam-A.hmm database (subset for viral proteins to reduce size).
URL: https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

### 6.2 KEGG Orthology Mapping

For small viral proteomes: map proteins to KEGG Orthology (KO) numbers and functional categories.
Do NOT perform pathway enrichment — sample size too small for statistical validity.

```python
def map_to_kegg_orthology(protein_annotation: str, conn: sqlite3.Connection) -> KEGGResult:
    # Step 1: Check cache
    cached = get_from_cache(conn, protein_annotation)
    if cached:
        return cached
    
    # Step 2: Query KEGG API
    url = f"https://rest.kegg.jp/find/genes/{protein_annotation}"
    response = requests.get(url, timeout=30)
    
    # Step 3: Cache result
    cache_result(conn, protein_annotation, response.json())
    return parse_kegg_response(response)
```

**Output:** Functional categorization table (not enrichment statistics):
- KEGG category: Metabolism, Information Processing, Cellular Processes, etc.
- Counts per category
- Bar chart of category distribution

### 6.3 Validation Criteria (Phase 6)
- [ ] SARS-CoV-2 Spike protein: Pfam domain PF09408 (Spike protein receptor-binding domain) detected
- [ ] MS2 coat protein: Pfam domain PF01819 detected
- [ ] KEGG categorization: at least 1 category mapped for annotated proteins
- [ ] Unannotated proteins logged as 'unannotated', not dropped

---

## PHASE 7 — ORA ENRICHMENT ANALYSIS (Workflow B)
**Duration:** 3 days
**Goal:** Statistically valid pathway enrichment using Fisher's exact test

### 7.1 Scientific Justification for ORA

ORA is appropriate for Workflow B because:
- Input: list of DEGs (logFC > 1 or < -1, padj < 0.05)
- Background: defined human gene universe (~19,000 protein-coding genes)
- Test: Fisher's exact test per pathway
- Correction: Benjamini-Hochberg FDR

This is biologically valid unlike GSEA on viral FASTA (which requires ranked expression data).

### 7.2 Algorithm: ORA
```python
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

def run_ora(deg_list: list, gene_sets: dict, universe: list) -> pd.DataFrame:
    results = []
    universe_size = len(universe)
    total_degs = len(deg_list)
    
    for pathway_name, pathway_genes in gene_sets.items:
        # Genes in pathway
        pathway_in_universe = [g for g in pathway_genes if g in universe]
        pathway_in_degs     = [g for g in pathway_genes if g in deg_list]
        
        # 2x2 contingency table
        # | in_pathway | not_in_pathway |
        # | in_degs    |   a            |   b    |  total_degs
        # | not_in_degs|   c            |   d    |  universe-total_degs
        a = len(pathway_in_degs)
        b = total_degs - a
        c = len(pathway_in_universe) - a
        d = universe_size - total_degs - c
        
        _, pvalue = fisher_exact([[a, b], [c, d]], alternative='greater')
        
        results.append({
            'pathway': pathway_name,
            'genes_in_pathway': a,
            'pvalue': pvalue,
            'gene_list': pathway_in_degs
        })
    
    # BH correction
    df = pd.DataFrame(results)
    _, padj, _, _ = multipletests(df['pvalue'], method='fdr_bh')
    df['padj'] = padj
    df['significant'] = df['padj'] < 0.05
    
    return df.sort_values('padj')
```

### 7.3 Gene Sets Sources
- KEGG pathways: download gene sets from KEGG REST API (cached locally)
- GO Biological Process: download from Gene Ontology consortium (goa_human.gaf)
- Minimum genes in pathway for testing: 15 (configurable)
- Maximum genes in pathway: 500 (configurable, avoids over-broad terms)

### 7.4 Validation Criteria (Phase 7)
- [ ] Test with known dataset: TP53 pathway gene list → Cell cycle pathway enriched
- [ ] FDR correction applied correctly (verify with known BH formula)
- [ ] Pathways with < 15 genes excluded from testing
- [ ] Background universe applied correctly (not just pathway genes)
- [ ] Results table: pathway, gene_count, pvalue, padj, gene_list

---

## PHASE 8 — VISUALIZATION ENGINE
**Duration:** 3 days
**Goal:** Publication-quality charts for all analysis outputs

### 8.1 Chart 1: ORF Position Map (Matplotlib)
```
Purpose: Show positions of all predicted ORFs along genome
Type: Horizontal bar chart (like NCBI genome viewer)
X-axis: Genome position (bp)
Y-axis: Reading frame (+1, +2, +3, -1, -2, -3)
Colors: Forward frames = blues, Reverse frames = reds
Annotations: ORF ID, length label for major ORFs
Output: orf_map.png, orf_map.svg
```

### 8.2 Chart 2: Volcano Plot (Plotly)
```
Purpose: Show all genes plotted by fold change and significance
X-axis: log2FoldChange
Y-axis: -log10(padj)
Point colors:
  Red = Upregulated (logFC > 1, padj < 0.05)
  Blue = Downregulated (logFC < -1, padj < 0.05)
  Gray = Not significant
Threshold lines:
  Vertical: x = -1 and x = 1
  Horizontal: y = -log10(0.05)
Hover: gene name, logFC, padj
Output: volcano_plot.png, volcano_plot.html (interactive)
```

### 8.3 Chart 3: KEGG Enrichment Bar Chart (Plotly)
```
Purpose: Show top enriched KEGG pathways
Type: Horizontal bar chart
X-axis: -log10(padj)
Y-axis: Pathway name (sorted by padj ascending)
Color: bar color encodes gene count (colorscale)
Show: top 20 significant pathways
Output: kegg_enrichment.png, kegg_enrichment.html
```

### 8.4 Chart 4: GO Bubble Chart (Plotly)
```
Purpose: Show GO term enrichment
X-axis: GeneRatio (genes_in_term / total_degs)
Y-axis: GO term name
Bubble size: gene count
Bubble color: padj value (colorscale from red to blue)
Output: go_bubble.html
```

### 8.5 Chart 5: Domain Frequency Chart (Matplotlib/Plotly)
```
Purpose: Show most common Pfam domains in viral proteome
Type: Vertical bar chart
X-axis: Pfam domain name
Y-axis: Count of proteins with domain
Output: domain_frequency.png
```

### 8.6 Validation Criteria (Phase 8)
- [ ] All charts generated without errors on test data
- [ ] Volcano plot: upregulated genes appear in upper right quadrant
- [ ] Enrichment chart: sorted correctly by padj
- [ ] All PNG outputs are minimum 300 DPI (publication quality)
- [ ] Interactive HTML charts render correctly in browser

---

## PHASE 9 — AI INTERPRETATION (Claude API Integration)
**Duration:** 1 day
**Goal:** Validated results sent to Claude for plain-language biological interpretation

### 9.1 Integration Design
```python
# ai/claude_interpreter.py
import anthropic

def generate_interpretation(analysis_results: dict, workflow: str) -> str:
    client = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY from environment
    
    prompt = build_interpretation_prompt(analysis_results, workflow)
    
    message = client.messages.create(
        model="claude-sonnet-4-20250514",
        max_tokens=2000,
        messages=[{"role": "user", "content": prompt}]
    )
    
    return message.content[0].text
```

### 9.2 Prompt Structure
```
WORKFLOW A PROMPT:
"You are a bioinformatics research assistant providing plain-language biological 
interpretation of computationally validated results. 

Analysis results:
- Genome: [name, length, GC content]
- ORFs predicted: [count, size range]
- Annotated proteins: [list with function descriptions]
- Pfam domains detected: [list of domains]
- KEGG functional categories: [category distribution]

Provide a 3-5 paragraph biological interpretation of these results, explaining:
1. What the genome structure suggests about this virus
2. What the identified proteins indicate about viral function
3. What the detected domains suggest about molecular mechanisms
4. Key biological insights from this analysis

IMPORTANT: Do not generate any statistics, p-values, or pathway assignments.
Only explain the biological significance of the provided results.
Do not make claims beyond what the data supports."
```

### 9.3 Validation Criteria (Phase 9)
- [ ] Claude API call succeeds with test data
- [ ] Interpretation references actual genes/domains from results (not hallucinated)
- [ ] AI does not generate fake statistics
- [ ] Disclaimer appears in all AI output sections
- [ ] Graceful failure if API unavailable (log warning, skip AI section, continue report)

---

## PHASE 10 — REPORT GENERATION
**Duration:** 2 days
**Goal:** Complete HTML and CSV reports generated automatically

### 10.1 HTML Report Structure (Jinja2 Template)
```html
<!DOCTYPE html>
<html>
<head>PathoScope AI Report — {{ run_id }}</head>
<body>
  <section id="header">
    Logo, run ID, timestamp, input file, workflow
  </section>
  
  <section id="qc-summary">
    QC metrics table, sequence statistics
  </section>
  
  <!-- Workflow A sections -->
  <section id="orf-results">
    ORF table, ORF map image
  </section>
  
  <section id="annotation">
    Protein annotation table, domain frequency chart
  </section>
  
  <!-- Workflow B sections -->
  <section id="deg-results">
    DEG classification table, volcano plot
  </section>
  
  <section id="enrichment">
    KEGG ORA table, GO ORA table, enrichment charts
  </section>
  
  <!-- Both workflows -->
  <section id="ai-interpretation">
    AI interpretation text with disclaimer
  </section>
  
  <section id="methods">
    Automated methods text (thresholds, tools, versions)
  </section>
  
  <section id="metadata">
    Run parameters, software versions, database versions
  </section>
</body>
</html>
```

### 10.2 CSV Outputs
- annotation_results.csv: orf_id, protein_sequence, best_hit, identity, evalue, pfam_domains, kegg_ko
- deg_results.csv: gene_symbol, logFC, pvalue, padj, classification
- pathway_results.csv: pathway_id, pathway_name, gene_count, pvalue, padj, significant, gene_list
- qc_metrics.csv: sequences_in, passed, removed_ambiguity, removed_duplicates, avg_length, gc_content

### 10.3 Validation Criteria (Phase 10)
- [ ] HTML report renders correctly in Chrome and Firefox
- [ ] All CSV files have correct headers and data
- [ ] Charts embedded in HTML as base64 (no external file dependencies)
- [ ] Methods section accurately describes parameters used
- [ ] Report file size < 10MB for standard analyses

---

## PHASE 11 — COMPLETE INTEGRATION TESTING
**Duration:** 3 days
**Goal:** Full end-to-end pipeline tested on all validation datasets

### 11.1 Test Dataset Matrix

| Test File | Workflow | Expected Result |
|---|---|---|
| ms2_bacteriophage.fasta | A | 4+ ORFs, coat protein annotated |
| phix174.fasta | A | 11+ ORFs, multiple domains |
| sars_cov2.fasta | A | 12+ ORFs, spike protein detected |
| influenza_a.fasta | A | 8+ ORFs (segmented genome) |
| invalid_chars.fasta | A | ValidationError, informative message |
| high_ambiguity.fasta | A | Sequences filtered, report shows removal |
| hgnc_symbols.txt | B | All symbols normalized, ORA run |
| ensembl_ids.txt | B | All converted to HGNC, ORA run |
| entrez_ids.txt | B | All converted to HGNC, ORA run |
| mixed_ids.csv | B | All types normalized, classified |
| expression_matrix.csv | B | DEG classified, volcano plot, ORA |
| corrupted.fastq | A | FileFormatError, informative message |

### 11.2 Offline Mode Test
```bash
# Disconnect internet, run:
python main.py --input sars_cov2.fasta --offline
# Verify: pipeline completes, log shows "using local cache"
```

### 11.3 Reproducibility Test
```bash
# Run same analysis twice:
python main.py --input expression_matrix.csv
python main.py --input expression_matrix.csv
# Compare: outputs must be byte-identical (except timestamps)
```

### 11.4 Performance Benchmarks
| Input | Target Runtime |
|---|---|
| MS2 FASTA (3.5 kb) | < 2 minutes |
| SARS-CoV-2 FASTA (30 kb) | < 10 minutes |
| Gene list (500 genes) | < 3 minutes |
| Expression matrix (2000 genes) | < 5 minutes |

---

## PHASE 12 — VIVA PREPARATION AND DOCUMENTATION
**Duration:** 2 days
**Goal:** Project ready for academic defense

### 12.1 Viva Demonstration Script
```
Demo 1 — Workflow A (5 minutes):
  python main.py --input tests/test_data/fasta/sars_cov2.fasta
  Show: live progress, open HTML report, show ORF map, domain chart
  
Demo 2 — Workflow B (5 minutes):
  python main.py --input tests/test_data/gene_lists/expression_matrix.csv
  Show: live progress, open HTML report, show volcano plot, enrichment chart

Demo 3 — Offline Mode (2 minutes):
  python main.py --input tests/test_data/fasta/ms2_bacteriophage.fasta --offline
  Show: runs without internet, uses cache

Demo 4 — Error Handling (2 minutes):
  python main.py --input tests/test_data/fasta/invalid_chars.fasta
  Show: informative error message, no crash
```

### 12.2 Scientific Justification Points (for Viva Questions)

**Q: Why don't you do GSEA?**
A: GSEA requires pre-ranked gene lists with continuous expression values (fold changes from RNA-seq experiments). Our Workflow A input is viral genomic sequence, not differential expression data. Applying GSEA to a list of annotated viral proteins would violate its core statistical assumptions. We use functional categorization for viral proteomes and reserve ORA (Fisher's exact test) for Workflow B where we have proper DEG lists with defined background universe.

**Q: Why is enrichment inappropriate for small viral genomes?**
A: Standard ORA requires sufficient gene counts for Fisher's exact test to be meaningful (typically > 500 background genes). MS2 bacteriophage has 4 proteins, SARS-CoV-2 has ~12. A "significant" pathway with 2/12 proteins has no statistical power and cannot be compared to pathway enrichment studies with thousands of genes. We therefore report functional categorization for Workflow A, not enrichment p-values.

**Q: How do you handle missing gene IDs?**
A: Our GeneIDNormalizer attempts conversion via: (1) local TSV mapping tables, (2) Ensembl REST API, (3) NCBI E-utilities, (4) MyGene.info. Genes that cannot be mapped are logged as 'unmapped' with reasons. They are excluded from enrichment analysis but remain in the DEG table with their original IDs, and the report shows the unmapped count and rate.

**Q: How do you ensure reproducibility?**
A: Every run produces a metadata.json recording the input file MD5, all parameters from config.yaml, software versions, and database versions. Any run can be re-executed by providing the same input and metadata.json. Our SQLite caching means offline runs return identical results regardless of API availability at time of rerun.

---

*End of Document 6 — Implementation Plan*
