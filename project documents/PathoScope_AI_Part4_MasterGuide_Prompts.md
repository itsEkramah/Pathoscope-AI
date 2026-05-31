# PathoScope AI v2.0
# MASTER BUILD GUIDE + COMPLETE PHASE-BY-PHASE PROMPTS
# Build From Zero — Detailed Instructions

---

# ═══════════════════════════════════════════════════
# MASTER BUILD GUIDE
# ═══════════════════════════════════════════════════

## HOW TO USE THIS GUIDE

This guide gives you one prompt per phase. Copy each prompt EXACTLY into Claude.
Before copying any prompt, read the SETUP RULES below.

## SETUP RULES (NEVER SKIP THESE)

Rule 1: Start every new chat with the SYSTEM CONTEXT PROMPT below.
Rule 2: Before each coding phase, paste the corresponding Phase Prompt.
Rule 3: After each phase, run the validation tests before moving to the next phase.
Rule 4: If something breaks, use the DEBUG PROMPT (section at bottom) — do NOT skip phases.
Rule 5: Never ask Claude to write all phases at once. One phase per conversation section.

---

## MASTER SYSTEM CONTEXT PROMPT
(Paste this FIRST in every new conversation about PathoScope AI)

```
I am a bioinformatics MSc student at the National Centre for Bioinformatics (NCB), 
Quaid-i-Azam University, Pakistan. I am building PathoScope AI v2.0 for my 
Functional Genomics course under Professor Sir Ghulam Abbas.

PROJECT IDENTITY:
Name: PathoScope AI v2.0
Type: Automated dual-workflow bioinformatics pipeline
Language: Python 3.10+
Purpose: Academic university project + viva demonstration

DUAL WORKFLOW DESIGN:
- Workflow A: Viral Genomics (Input: FASTA/FASTQ → QC → ORF → Protein → Pfam/KEGG → Report)
- Workflow B: Functional Genomics (Input: Gene List/Expression Matrix → Normalize → DEG → ORA → Report)

KEY SCIENTIFIC CONSTRAINTS (never violate these):
1. Do NOT apply GSEA or enrichment analysis to viral FASTA files
2. Do NOT apply ORA to small viral proteomes (< 50 proteins)
3. ORA is ONLY valid for Workflow B with proper DEG lists and background universe
4. AI (Claude API) interprets VALIDATED outputs only — never generates statistics
5. All statistical thresholds are configurable via config.yaml
6. System must work in offline mode using SQLite cache
7. Every run generates: HTML report + CSV results + metadata.json + run.log

TECH STACK:
Python 3.10+ | BioPython | pyhmmer | pandas | scipy | statsmodels
matplotlib | plotly | seaborn | sqlite3 | requests | PyYAML
anthropic | jinja2 | pytest

FOLDER STRUCTURE:
pathoscope_ai/
├── main.py
├── config.yaml
├── requirements.txt
├── core/ (input_manager, file_detector, validator, workflow_router, config_loader)
├── workflow_a/ (qc_engine, sequence_cleaner, orf_predictor, protein_translator, 
│               annotation_engine, pfam_mapper, kegg_mapper_a, functional_categorizer)
├── workflow_b/ (gene_id_normalizer, expression_loader, statistical_filter,
│               deg_classifier, ora_engine, go_enrichment, kegg_mapper_b)
├── visualization/ (orf_visualizer, volcano_plotter, enrichment_plotter, report_builder)
├── ai/ (claude_interpreter, prompt_templates)
├── utils/ (logger, metadata_writer, db_cache, api_client, exceptions)
├── databases/ (pfam_local/, kegg_cache.sqlite, ensembl_hgnc_map.tsv, 
│              entrez_hgnc_map.tsv, go_annotations.sqlite, gene_universe.tsv)
├── tests/ (test_data/, unit tests per module)
└── output/ (run_YYYYMMDD_HHMMSS/ per run)

Act as a senior bioinformatics software engineer.
Write clean, modular Python code with:
- docstrings for every function
- type hints
- proper error handling (try/except with custom exceptions)
- logging at every major step
- no hardcoded values (everything from config.yaml)
Do NOT generate fake biological data or fake statistics.
```

---

# ═══════════════════════════════════════════════════
# PHASE-BY-PHASE PROMPTS
# ═══════════════════════════════════════════════════

---

## PHASE 0 PROMPT — ENVIRONMENT AND SKELETON SETUP

```
[PHASE 0 — ENVIRONMENT AND ARCHITECTURE SETUP]

Task: Create the complete project skeleton for PathoScope AI v2.0.
Do NOT implement any biology yet. Create only the structure and infrastructure.

CREATE THE FOLLOWING FILES IN ORDER:

FILE 1: requirements.txt
Include exact versions for:
biopython>=1.81
pyhmmer>=0.10.0
pandas>=2.0.0
numpy>=1.24.0
scipy>=1.11.0
statsmodels>=0.14.0
matplotlib>=3.7.0
plotly>=5.15.0
seaborn>=0.12.0
requests>=2.31.0
PyYAML>=6.0
python-dotenv>=1.0.0
anthropic>=0.25.0
jinja2>=3.1.0
WeasyPrint>=60.0
pytest>=7.4.0
pytest-cov>=4.1.0
tqdm>=4.65.0

FILE 2: environment.yml
Conda environment named: pathoscope_ai
Python 3.10
Include all packages from requirements.txt

FILE 3: config.yaml
Include ALL of the following sections with these exact keys:
project:
  name: "PathoScope AI"
  version: "2.0.0"
  output_dir: "./output"
  log_level: "INFO"

workflow:
  mode: "auto"

quality_control:
  min_sequence_length: 50
  max_ambiguity_ratio: 0.10
  min_fastq_quality: 20
  remove_duplicates: true

orf_prediction:
  min_orf_length: 100
  genetic_code: 1
  both_strands: true
  remove_nested: false

statistical:
  logfc_up_threshold: 1.0
  logfc_down_threshold: -1.0
  pvalue_cutoff: 0.05
  fdr_cutoff: 0.05
  fdr_method: "benjamini-hochberg"
  min_pathway_genes: 15
  max_pathway_genes: 500

annotation:
  blast_evalue: 0.001
  blast_identity: 30.0
  blast_coverage: 50.0
  pfam_evalue: 0.01

api:
  kegg_base_url: "https://rest.kegg.jp"
  ensembl_base_url: "https://rest.ensembl.org"
  ncbi_base_url: "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  mygene_base_url: "https://mygene.info/v3"
  request_timeout: 30
  retry_attempts: 3
  retry_delay: 2
  offline_mode: false

cache:
  sqlite_path: "./databases/pathoscope.sqlite"
  kegg_ttl_days: 7
  gene_id_ttl_days: 30

ai:
  model: "claude-sonnet-4-20250514"
  max_tokens: 2000
  enabled: true

FILE 4: .env.example
ANTHROPIC_API_KEY=your_anthropic_api_key_here
NCBI_API_KEY=your_ncbi_api_key_here_optional

FILE 5: utils/exceptions.py
Create custom exception hierarchy:
- PathoScopeError (base)
- InputValidationError
- FileFormatError
- SequenceQualityError
- GeneNormalizationError
- APIConnectionError
- DatabaseError
- StatisticalError
- ReportGenerationError

Each exception should:
- Accept message and optional details dict
- Include __str__ with formatted message

FILE 6: utils/logger.py
Create a logging module with:
- setup_logger(run_id: str, output_dir: str) -> logging.Logger
  Creates logger that writes to both console (INFO level) and
  file (DEBUG level) at output_dir/run.log
  Format: [YYYY-MM-DD HH:MM:SS] [LEVEL] [module_name] message
  
- log_step(logger, step_name: str, status: str, details: dict = None)
  Standardized step logging: [STEP: QC] [STATUS: COMPLETE] sequences_in=100 passed=98

FILE 7: utils/db_cache.py
Create SQLite database manager with:
- initialize_database(db_path: str) -> sqlite3.Connection
  Creates all tables if they don't exist
  Returns connection

- create_all_tables(conn: sqlite3.Connection)
  Creates: runs, qc_results, orf_results, annotation_results,
  gene_normalization_results, deg_results, pathway_results,
  visualization_metadata, ai_interpretations, kegg_cache,
  gene_id_cache, go_term_cache
  (Use schema exactly as specified in the Backend Schema document)

- purge_expired_cache(conn: sqlite3.Connection)
  Deletes rows where expires_at < CURRENT_TIMESTAMP

- get_cached(conn, table: str, key: str) -> dict or None
- set_cached(conn, table: str, key: str, data: dict, ttl_days: int)

FILE 8: core/config_loader.py
Create config loader with:
- load_config(config_path: str) -> dict
  Loads config.yaml, returns as dict, raises FileNotFoundError if missing

- validate_config(config: dict) -> bool
  Checks all required keys present, types correct, values in valid ranges
  Raises PathoScopeError if invalid

- get_threshold(config: dict, *keys) -> any
  Safely access nested config keys: get_threshold(config, 'statistical', 'fdr_cutoff')

FILE 9: main.py
CLI entry point using argparse:
Arguments:
  --input     (required): Path to input file
  --workflow  (optional): 'A', 'B', or 'auto' (default: 'auto')
  --config    (optional): Path to config file (default: './config.yaml')
  --offline   (flag): Force offline mode
  --strict    (flag): Reject any invalid sequences (default: False = skip and log)
  --version   (flag): Print version and exit

Main flow:
  1. Parse args
  2. Load config
  3. Set up logger with new run_id
  4. Display startup banner
  5. Initialize database
  6. Route to appropriate workflow function (stub for now)
  7. Display completion message

Display startup banner:
  ╔══════════════════════════════════════╗
  ║        PathoScope AI v2.0            ║
  ║  Functional Genomics Analysis Suite  ║
  ║  NCB, Quaid-i-Azam University        ║
  ╚══════════════════════════════════════╝

Also create ALL __init__.py files for:
core/, workflow_a/, workflow_b/, visualization/, ai/, utils/, tests/

After generating all files, show me:
1. The complete folder structure
2. How to run: python main.py --help
3. Expected output of --help
4. Installation command: pip install -r requirements.txt
```

---

## PHASE 1 PROMPT — INPUT HANDLING AND WORKFLOW ROUTING

```
[PHASE 1 — INPUT HANDLING AND WORKFLOW ROUTING]

I have completed Phase 0. All skeleton files exist.
Now implement the input detection and routing system.

CREATE FILE: core/file_detector.py

Implement these functions with complete logic (NOT stubs):

1. detect_format(file_path: str) -> str
   Algorithm:
   a. Check file extension
   b. Read first 10 lines of file
   c. Determine actual format from content
   d. Return one of: 'fasta', 'fastq', 'csv_expression', 'tsv_expression',
      'txt_genes', 'csv_genes', 'unknown'
   
   Detection rules:
   - FASTA: lines starting with '>' followed by sequence lines
   - FASTQ: lines 1,5,9... start with '@', line 3 starts with '+'
   - Expression matrix: CSV/TSV with headers containing 
     ('log2FoldChange' or 'logFC') AND 'pvalue'
   - Gene list (TXT): each line is a single gene identifier
     (ENSG pattern, digit string, or gene symbol pattern)
   - Gene list (CSV): CSV with single column or first column only containing gene IDs

2. is_fasta(file_path: str) -> bool
3. is_fastq(file_path: str) -> bool  
4. is_expression_matrix(file_path: str) -> bool
5. is_gene_list_txt(file_path: str) -> bool
6. validate_file_accessible(file_path: str) -> None
   Raises InputValidationError if: not exists, empty, unreadable, > 500MB

CREATE FILE: core/workflow_router.py

1. determine_workflow(detected_format: str, override: str = None) -> str
   Returns 'A' or 'B'
   'fasta' or 'fastq' → 'A'
   everything else → 'B'
   If override is 'A' or 'B', return that directly
   If override is 'auto' or None, use format detection

CREATE FILE: core/input_manager.py

Create InputObject dataclass:
@dataclass
class InputObject:
    file_path: str
    detected_format: str
    workflow: str
    file_size_bytes: int
    file_md5: str

Implement:
1. load_input(file_path: str, workflow_override: str = 'auto') -> InputObject
   Complete flow: validate → detect → route → return InputObject

2. compute_md5(file_path: str) -> str
   Returns MD5 hash of file for reproducibility tracking

UPDATE: main.py
Connect input loading to main pipeline
After loading InputObject, print:
  Input file: [filename]
  Format detected: [format]
  File size: [X MB]
  Routing to: Workflow [A or B]

CREATE TEST FILE: tests/test_input_manager.py

Write pytest tests:
- test_detect_fasta() — use ms2 test data
- test_detect_fastq() — use valid_reads test data
- test_detect_expression_csv() — create small inline test CSV
- test_detect_gene_list_txt() — create small inline test TXT
- test_reject_empty_file() — create tmp empty file
- test_reject_nonexistent_file() — use fake path
- test_route_fasta_to_workflow_a()
- test_route_gene_list_to_workflow_b()
- test_route_expression_to_workflow_b()
- test_compute_md5_consistent() — same file same hash

CREATE: tests/test_data/gene_lists/hgnc_symbols.txt
Content: 20 real human gene symbols, one per line:
TP53, BRCA1, BRCA2, MYC, EGFR, PTEN, RB1, APC, VHL, MLH1,
MSH2, KRAS, BRAF, PIK3CA, CDKN2A, ATM, CDH1, SMAD4, STK11, PTEN

CREATE: tests/test_data/gene_lists/ensembl_ids.txt
Content: 10 Ensembl gene IDs, one per line:
ENSG00000141510 (TP53), ENSG00000012048 (BRCA1), etc.

CREATE: tests/test_data/gene_lists/expression_matrix.csv
Content: 50-row CSV with columns:
gene_id, log2FoldChange, pvalue
Include: 10 upregulated, 10 downregulated, 30 not significant

Show me the complete code for all files.
After generating, show me: pytest tests/test_input_manager.py -v
and what the expected output should be.
```

---

## PHASE 2 PROMPT — QUALITY CONTROL ENGINE

```
[PHASE 2 — QUALITY CONTROL ENGINE (Workflow A)]

Phase 0 and 1 are complete and all tests pass.
Now implement the QC engine for FASTA and FASTQ input.

CREATE FILE: workflow_a/qc_engine.py

The QCResult dataclass should contain:
@dataclass
class QCResult:
    sequences_input: int
    sequences_passed: int
    sequences_removed: int
    removed_ambiguity: int
    removed_duplicates: int  
    removed_length: int
    removed_invalid_chars: int
    avg_length: float
    min_length: int
    max_length: int
    total_bases: int
    gc_content_percent: float
    avg_ambiguity_ratio: float
    passed_sequences: list  # BioPython SeqRecord list
    warnings: list[str]

Implement ALL of these functions with COMPLETE logic:

1. run_fasta_qc(file_path: str, config: dict, logger) -> QCResult
   Complete FASTA quality control pipeline:
   Step 1: Parse sequences with BioPython SeqIO.parse(file, 'fasta')
   Step 2: Uppercase all sequences
   Step 3: check_invalid_characters() → remove, log invalid chars found
   Step 4: calculate_ambiguity_ratio() → filter if > config threshold
   Step 5: filter by length → remove if < config min_sequence_length
   Step 6: remove_duplicates() → collapse exact duplicate sequences
   Step 7: Calculate summary statistics on passed sequences
   Step 8: Return QCResult

2. check_invalid_characters(seq_str: str) -> list[str]
   Valid IUPAC nucleotide codes: ACGTRYSWKMBDHVN (uppercase)
   Returns list of unique invalid characters found
   Empty list means sequence is valid

3. calculate_ambiguity_ratio(seq_str: str) -> float
   ambiguity = count of N characters / total length
   Return ratio as float (0.0 to 1.0)

4. filter_by_ambiguity(sequences: list, max_ratio: float) -> tuple[list, int]
   Returns (passing_sequences, removed_count)

5. filter_by_length(sequences: list, min_len: int) -> tuple[list, int]
   Returns (passing_sequences, removed_count)

6. remove_duplicates(sequences: list) -> tuple[list, int]
   Use MD5 of sequence string for deduplication
   Keep first occurrence of each sequence
   Returns (unique_sequences, removed_count)

7. calculate_gc_content(sequences: list) -> float
   Returns average GC content across all sequences

8. run_fastq_qc(file_path: str, config: dict, logger) -> QCResult
   FASTQ-specific QC:
   Use BioPython SeqIO.parse(file, 'fastq')
   Additional metrics: avg quality score, reads below quality threshold
   Filter reads with mean Phred quality < config min_fastq_quality

9. calculate_mean_phred_quality(record) -> float
   Uses BioPython letter_annotations['phred_quality']

CREATE FILE: workflow_a/sequence_cleaner.py

1. clean_sequence(seq_str: str) -> str
   Uppercase, strip whitespace, remove non-sequence characters

2. write_cleaned_fasta(sequences: list, output_path: str)
   Write BioPython SeqRecord list to FASTA file

3. generate_qc_csv(qc_result: QCResult, output_path: str)
   Write QC metrics to CSV with columns:
   metric, value, threshold, status (PASS/FAIL)

CREATE TEST FILE: tests/test_qc_engine.py

Tests:
- test_valid_ms2_fasta_passes_qc()
- test_invalid_chars_rejected()
  Input: sequence with 'X' and 'J' characters
  Expected: removed_invalid_chars > 0
- test_high_ambiguity_filtered()
  Input: sequence with >10% N content
  Expected: removed_ambiguity > 0
- test_duplicates_collapsed()
  Input: 3 identical sequences
  Expected: 1 sequence returned
- test_length_filter()
  Input: sequences shorter than min_sequence_length
  Expected: removed_length > 0
- test_gc_content_calculated()
- test_fastq_quality_filter()

ALSO CREATE: tests/test_data/fasta/

Create these small test files inline in the test code using tmpfiles:

invalid_chars test: 
>seq_001
ATGCXYZATGC

high_ambiguity test (>10% N):
>seq_001
ATGCNNNNNNATGCATGCNNNNN

duplicates test (3 identical):
>seq_001
ATGCATGCATGCATGCATGC
>seq_002
ATGCATGCATGCATGCATGC
>seq_003
ATGCATGCATGCATGCATGC

After implementing, run me through:
What happens step by step when I run:
python main.py --input tests/test_data/fasta/ms2_bacteriophage.fasta
(up to the point where we have QC results)
```

---

## PHASE 3 PROMPT — ORF PREDICTION AND PROTEIN TRANSLATION

```
[PHASE 3 — ORF PREDICTION AND PROTEIN TRANSLATION]

Phases 0, 1, 2 complete. Now implement ORF prediction.

CREATE FILE: workflow_a/orf_predictor.py

Create the ORF dataclass:
@dataclass  
class ORF:
    orf_id: str          # orf_001, orf_002...
    source_sequence_id: str
    start_position: int  # 0-based
    stop_position: int   # 0-based
    length_bp: int
    reading_frame: int   # 1, 2, 3, -1, -2, -3
    strand: str          # 'forward', 'reverse'
    dna_sequence: str
    protein_sequence: str  # filled by translator
    start_codon: str     # usually 'ATG'
    stop_codon: str

Implement ALL functions:

1. predict_orfs(sequence_records: list, config: dict, logger) -> list[ORF]
   Main ORF prediction pipeline
   For each sequence record:
     generate_all_frames() → find_orfs_in_all_frames() → filter → rank → assign IDs

2. get_reverse_complement(sequence: str) -> str
   Manual implementation using complement dict:
   complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', ...}
   Then reverse the complemented string

3. generate_reading_frames(sequence: str) -> dict[str, str]
   Returns dict: {'frame_+1': seq[0:], 'frame_+2': seq[1:], 'frame_+3': seq[2:],
                  'frame_-1': rc[0:], 'frame_-2': rc[1:], 'frame_-3': rc[2:]}

4. find_orfs_in_frame(frame_sequence: str, frame_id: str, 
                       min_length: int, original_seq: str) -> list[ORF]
   Algorithm:
   START_CODONS = ['ATG']
   STOP_CODONS = ['TAA', 'TAG', 'TGA']
   
   codon_table = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
   
   i = 0
   while i < len(codon_table):
     if codon_table[i] in START_CODONS:
       start_idx = i * 3
       for j in range(i+1, len(codon_table)):
         if codon_table[j] in STOP_CODONS:
           stop_idx = (j * 3) + 3
           orf_len = stop_idx - start_idx
           if orf_len >= min_length:
             record ORF
           i = j + 1  # start scanning after stop codon
           break
     else:
       i += 1

5. filter_orfs_by_length(orfs: list, min_len: int) -> list[ORF]

6. assign_orf_ids(orfs: list) -> list[ORF]
   Sort by length descending, assign orf_001, orf_002... etc.

7. map_frame_position_to_genome(frame_position: int, frame_id: str, seq_len: int) -> tuple
   Convert position in frame back to genomic coordinates

CREATE FILE: workflow_a/protein_translator.py

1. translate_orf(dna_sequence: str, genetic_code: int = 1) -> str
   Use BioPython:
   from Bio.Seq import Seq
   protein = str(Seq(dna_sequence).translate(table=genetic_code, to_stop=True))
   Strip trailing stop codon if present

2. translate_all_orfs(orfs: list[ORF], genetic_code: int) -> list[ORF]
   Add protein_sequence to each ORF
   Skip translation if protein starts with non-M (flag as unusual)

3. write_orf_dna_fasta(orfs: list[ORF], output_path: str)
   Write DNA sequences of all ORFs to FASTA
   Header format: >orf_001 frame=+1 start=265 stop=21555 length_bp=21291

4. write_protein_fasta(orfs: list[ORF], output_path: str)
   Write protein sequences to FASTA
   Header format: >orf_001_protein frame=+1 length_aa=7097

CREATE TEST FILE: tests/test_orf_predictor.py

Tests:
- test_reverse_complement_correct()
  Input: 'ATGC' → Expected: 'GCAT'
  Input: 'AAAAAA' → Expected: 'TTTTTT'

- test_find_orf_simple()
  Input: 'ATGAAATAA' (ATG + AAA + TAA stop)
  Expected: 1 ORF found, length 9 bp, protein 'MK'

- test_six_frames_generated()
  Input: any 30bp sequence
  Expected: dict with 6 frame keys

- test_minimum_orf_length_filter()
  min_length = 30 bp
  Input: sequence with ORFs of various lengths
  Expected: only ORFs >= 30 bp returned

- test_orf_prediction_ms2()
  Input: MS2 bacteriophage known sequence fragment
  Expected: at least 1 ORF found

- test_protein_starts_with_methionine()
  All translated proteins must start with 'M'

- test_protein_translation_standard_code()
  ATG → M, TAA → stop (not in output), AAA → K

IMPORTANT IMPLEMENTATION NOTE:
Do NOT use Smith-Waterman implementation from scratch.
Use BioPython's built-in tools only.
Reason: Engineering cost > educational value for this course scope.
```

---

## PHASE 4 PROMPT — GENE ID NORMALIZATION

```
[PHASE 4 — GENE ID NORMALIZATION (Workflow B)]

Now implement Gene ID Normalization. This is one of the most important modules 
because the professor explicitly requires HGNC/Ensembl/Entrez support.

CREATE FILE: workflow_b/gene_id_normalizer.py

Create dataclasses:
@dataclass
class GeneIDResult:
    original_id: str
    detected_type: str      # 'hgnc', 'ensembl', 'entrez', 'unknown'
    hgnc_symbol: str        # None if normalization failed
    conversion_method: str  # 'identity', 'local_table', 'ensembl_api', 
                            # 'ncbi_api', 'mygene_api', 'failed'
    success: bool
    error_reason: str       # None if successful

@dataclass
class NormalizationReport:
    total_genes: int
    successfully_normalized: int
    failed_normalization: int
    hgnc_input: int
    ensembl_input: int
    entrez_input: int
    unknown_input: int
    results: list[GeneIDResult]

Implement ALL functions:

1. detect_id_type(gene_id: str) -> str
   Use regex patterns:
   ENSEMBL_PATTERN = r'^ENSG\d{11}(\.\d+)?$'  # allows version suffix
   ENTREZ_PATTERN  = r'^\d{1,8}$'
   HGNC_PATTERN    = r'^[A-Z][A-Z0-9\-]{1,}$'  # simplified
   
   Check in order: ensembl → entrez → hgnc → unknown
   Return type string

2. load_local_ensembl_map(tsv_path: str) -> dict[str, str]
   Load databases/ensembl_hgnc_map.tsv
   Return dict: {ensembl_id: hgnc_symbol}
   Handle missing file gracefully (return empty dict, log warning)

3. load_local_entrez_map(tsv_path: str) -> dict[str, str]
   Load databases/entrez_hgnc_map.tsv
   Return dict: {entrez_id: hgnc_symbol}

4. normalize_gene_list(gene_ids: list[str], config: dict, 
                        db_conn, logger) -> NormalizationReport
   Main normalization function:
   
   For each gene_id:
     Step 1: Strip whitespace, uppercase
     Step 2: detect_id_type()
     
     If type == 'hgnc':
       result = GeneIDResult(..., hgnc_symbol=gene_id, method='identity', success=True)
     
     If type == 'ensembl':
       Step 3: Check gene_id_cache in SQLite
       Step 4: Check local ensembl map
       Step 5: If not found and online: query_ensembl_api()
       Step 6: Cache result
     
     If type == 'entrez':
       Step 3: Check gene_id_cache in SQLite
       Step 4: Check local entrez map
       Step 5: If not found and online: query_ncbi_api()
       Step 6: Cache result
     
     If type == 'unknown':
       Log as unmapped, continue

5. query_ensembl_api(ensembl_id: str, config: dict) -> str or None
   URL: {ensembl_base_url}/xrefs/id/{ensembl_id}?content-type=application/json
   Parse response for HGNC symbol
   Return symbol or None
   Handle: timeout, 404, 429 (rate limit: sleep and retry)

6. query_ncbi_api(entrez_id: str, config: dict) -> str or None
   URL: {ncbi_base_url}efetch.fcgi?db=gene&id={entrez_id}&rettype=gene_table&retmode=text
   Parse response to extract gene symbol
   Return symbol or None

7. query_mygene_api(gene_id: str, id_type: str, config: dict) -> str or None
   URL: {mygene_base_url}/gene/{gene_id}?fields=symbol
   Last resort fallback
   Return symbol or None

8. generate_normalization_csv(report: NormalizationReport, output_path: str)
   Write CSV with columns:
   original_id, detected_type, hgnc_symbol, conversion_method, success

CREATE SAMPLE LOCAL MAPPING TABLE: databases/ensembl_hgnc_map.tsv
Create a sample file with at least 50 real Ensembl→HGNC mappings.
Include these essential ones:
ENSG00000141510  TP53
ENSG00000012048  BRCA1
ENSG00000139618  BRCA2
ENSG00000136997  MYC
ENSG00000146648  EGFR
ENSG00000171862  PTEN
ENSG00000139687  RB1
ENSG00000134982  APC
ENSG00000134086  VHL
(add 40 more from known cancer genes)

CREATE SAMPLE LOCAL MAPPING TABLE: databases/entrez_hgnc_map.tsv
Same genes with Entrez IDs:
7157  TP53
672   BRCA1
675   BRCA2
4609  MYC
1956  EGFR
5728  PTEN
5925  RB1
324   APC
7428  VHL
(add 40 more)

CREATE TEST FILE: tests/test_gene_normalizer.py

Tests:
- test_detect_ensembl_id()
  Input: 'ENSG00000141510' → Expected: 'ensembl'

- test_detect_entrez_id()
  Input: '7157' → Expected: 'entrez'

- test_detect_hgnc_symbol()
  Input: 'TP53' → Expected: 'hgnc'

- test_detect_unknown()
  Input: 'XXXXXXX123' → Expected: 'unknown'

- test_ensembl_to_hgnc_local_table()
  Input: 'ENSG00000141510' (TP53's Ensembl ID)
  Expected: 'TP53', method='local_table'

- test_entrez_to_hgnc_local_table()
  Input: '7157' (TP53's Entrez ID)
  Expected: 'TP53', method='local_table'

- test_hgnc_identity()
  Input: 'TP53'
  Expected: 'TP53', method='identity'

- test_normalization_report_counts()
  Input: mixed list of 5 HGNC, 5 Ensembl, 5 Entrez IDs
  Expected: 15 total, 15 successful (all in local table)

- test_offline_mode_no_api_calls()
  With offline_mode=True, api calls should not be made
  (Mock requests to verify no API calls)
```

---

## PHASE 5 PROMPT — STATISTICAL FILTERING AND DEG CLASSIFICATION

```
[PHASE 5 — STATISTICAL FILTERING AND DEG CLASSIFICATION]

Implement the core statistical analysis for Workflow B.
This is scientifically critical — implement carefully.

CREATE FILE: workflow_b/expression_loader.py

1. load_expression_matrix(file_path: str, logger) -> pd.DataFrame
   Load CSV or TSV file
   Auto-detect delimiter from extension
   
2. validate_expression_matrix(df: pd.DataFrame) -> list[str]
   Returns list of validation errors (empty = valid)
   Required checks:
   - Must have gene_id column (or first column treated as gene_id)
   - Must have log2FoldChange or logFC column
   - Must have pvalue column
   - No negative p-values
   - P-values must be 0 to 1
   - logFC must be numeric
   - At least 10 rows required
   
3. standardize_column_names(df: pd.DataFrame) -> pd.DataFrame
   Rename alternate column names to standard names:
   'logFC' → 'log2FoldChange'
   'log2fc' → 'log2FoldChange'  
   'FDR' → 'padj'
   'adj.P.Val' → 'padj'
   'p.value' → 'pvalue'
   'P.Value' → 'pvalue'

4. load_gene_list(file_path: str, logger) -> pd.DataFrame
   Load TXT or CSV gene list
   Return DataFrame with single column 'gene_id'

CREATE FILE: workflow_b/statistical_filter.py

1. apply_fdr_correction(df: pd.DataFrame, method: str = 'fdr_bh') -> pd.DataFrame
   ONLY apply if 'padj' column is NOT already in DataFrame
   Use: from statsmodels.stats.multitest import multipletests
   reject, padj_values, _, _ = multipletests(df['pvalue'].values, method=method)
   df['padj'] = padj_values
   df['bh_corrected'] = True  # flag that we calculated this
   Log: "BH FDR correction applied to X genes"
   
   IMPORTANT: Never apply FDR correction twice
   IMPORTANT: Log clearly if padj already existed in input

2. apply_logfc_filter(df: pd.DataFrame, up_threshold: float, 
                       down_threshold: float) -> tuple[pd.DataFrame, pd.DataFrame]
   Returns (up_candidates, down_candidates)
   up_candidates: log2FoldChange > up_threshold
   down_candidates: log2FoldChange < down_threshold

3. apply_significance_filter(df: pd.DataFrame, fdr_cutoff: float) -> pd.DataFrame
   Filter: padj < fdr_cutoff
   Log: how many genes pass

CREATE FILE: workflow_b/deg_classifier.py

1. classify_degs(df: pd.DataFrame, config: dict, logger) -> pd.DataFrame
   Add 'classification' column with values:
   'upregulated'     if log2FoldChange >  logfc_up AND padj < fdr_cutoff
   'downregulated'   if log2FoldChange < logfc_down AND padj < fdr_cutoff
   'not_significant' otherwise
   
   IMPORTANT: use the thresholds from config, NOT hardcoded values
   
   Log summary:
   "DEG Classification complete:"
   "  Upregulated: N genes"
   "  Downregulated: N genes"
   "  Not significant: N genes"
   "  Total: N genes"

2. get_deg_lists(df: pd.DataFrame) -> dict
   Returns: {
     'upregulated': [list of gene symbols],
     'downregulated': [list of gene symbols],
     'all_degs': [combined list],
     'not_significant': [list]
   }

3. generate_deg_csv(df: pd.DataFrame, output_path: str)
   Write CSV with columns:
   gene_symbol, log2FoldChange, pvalue, padj, classification
   Sorted by: padj ascending

4. generate_deg_summary(df: pd.DataFrame) -> dict
   Returns dict of counts and statistics for report

CREATE TEST FILE: tests/test_statistical_filter.py

Tests:
- test_fdr_correction_applied()
  Input: DataFrame without padj column
  Expected: padj column added, all values 0-1

- test_fdr_not_applied_twice()
  Input: DataFrame with existing padj column
  Expected: padj unchanged, bh_corrected flag = False

- test_upregulated_classification()
  Input: gene with logFC=2.5, padj=0.001
  Expected: 'upregulated'

- test_downregulated_classification()
  Input: gene with logFC=-2.0, padj=0.003
  Expected: 'downregulated'

- test_not_significant_high_pvalue()
  Input: gene with logFC=2.5, padj=0.2 (above threshold)
  Expected: 'not_significant'

- test_not_significant_low_logfc()
  Input: gene with logFC=0.3, padj=0.001 (below logFC threshold)
  Expected: 'not_significant'

- test_deg_counts_sum_to_total()
  up + down + not_significant must equal total input genes

- test_configurable_thresholds()
  Use logfc_threshold=2.0, verify classification changes vs default 1.0

SCIENTIFIC VALIDATION:
After implementation, test with this known case:
Input: TP53 with logFC=2.3, pvalue=0.001
Expected after BH correction with a list of 100 genes:
  padj will be approximately 0.001 * 100 / rank_position
  Classification should be 'upregulated' if padj < 0.05
```

---

## PHASE 6 PROMPT — ORA ENRICHMENT

```
[PHASE 6 — ORA ENRICHMENT (Workflow B)]

Implement statistically valid pathway enrichment.
This is ONLY for Workflow B, NEVER for viral FASTA inputs.

CREATE FILE: workflow_b/ora_engine.py

1. load_kegg_gene_sets(config: dict, db_conn, logger) -> dict[str, list[str]]
   Load KEGG pathway gene sets
   Priority: SQLite cache → KEGG API
   
   KEGG API strategy:
   Step 1: GET https://rest.kegg.jp/list/pathway/hsa (human pathways list)
   Step 2: For each pathway ID, GET https://rest.kegg.jp/link/hsa/{pathway_id}
   Step 3: Parse response to extract gene symbols
   Step 4: Cache in SQLite with 7-day TTL
   
   Return: {'hsa04110': ['CDK1', 'CDK2', 'TP53', ...], ...}

2. load_gene_universe(universe_path: str) -> list[str]
   Load databases/gene_universe.tsv
   Return list of all human protein-coding gene symbols (~19,000 genes)
   This is the background for Fisher's exact test

3. run_ora(deg_list: list[str], gene_sets: dict[str, list[str]], 
           universe: list[str], config: dict, logger) -> pd.DataFrame
   
   COMPLETE IMPLEMENTATION:
   
   from scipy.stats import fisher_exact
   from statsmodels.stats.multitest import multipletests
   
   universe_set = set(universe)
   deg_set = set(deg_list)
   total_degs = len(deg_set)
   universe_size = len(universe_set)
   
   results = []
   
   for pathway_name, pathway_genes in gene_sets.items():
       pathway_in_universe = [g for g in pathway_genes if g in universe_set]
       n_pathway_in_universe = len(pathway_in_universe)
       
       # Filter pathways too small or too large
       if n_pathway_in_universe < config['statistical']['min_pathway_genes']:
           continue
       if n_pathway_in_universe > config['statistical']['max_pathway_genes']:
           continue
       
       pathway_in_degs = [g for g in pathway_genes if g in deg_set]
       n_pathway_in_degs = len(pathway_in_degs)
       
       # 2x2 contingency table:
       # Rows: DEG / Not-DEG
       # Cols: In pathway / Not in pathway
       a = n_pathway_in_degs                              # DEG in pathway
       b = total_degs - a                                 # DEG not in pathway
       c = n_pathway_in_universe - a                      # Not-DEG in pathway
       d = universe_size - total_degs - c                 # Not-DEG not in pathway
       
       _, pvalue = fisher_exact([[a, b], [c, d]], alternative='greater')
       
       gene_ratio = a / total_degs if total_degs > 0 else 0
       bg_ratio = n_pathway_in_universe / universe_size
       
       results.append({
           'pathway_name': pathway_name,
           'genes_in_pathway': a,
           'universe_genes_in_pathway': n_pathway_in_universe,
           'total_degs': total_degs,
           'universe_size': universe_size,
           'pvalue': pvalue,
           'gene_ratio': gene_ratio,
           'background_ratio': bg_ratio,
           'gene_list': '|'.join(pathway_in_degs)
       })
   
   df = pd.DataFrame(results)
   
   if df.empty:
       logger.warning("No pathways passed size filters for ORA")
       return df
   
   # BH correction
   _, padj, _, _ = multipletests(df['pvalue'].values, method='fdr_bh')
   df['padj'] = padj
   df['significant'] = df['padj'] < config['statistical']['fdr_cutoff']
   
   # Sort by padj
   df = df.sort_values('padj').reset_index(drop=True)
   
   logger.info(f"ORA complete: {df['significant'].sum()} significant pathways (FDR < {config['statistical']['fdr_cutoff']})")
   
   return df

4. run_go_enrichment(deg_list: list[str], go_db_conn, config: dict, logger) -> pd.DataFrame
   Same algorithm as run_ora() but using GO terms from go_term_cache
   Run for: biological_process, molecular_function (separately)

5. generate_enrichment_csv(results: pd.DataFrame, output_path: str)
   Write CSV with all columns
   Include only significant results in "significant_pathways.csv"
   Include all results in "all_pathway_results.csv"

CREATE: databases/gene_universe.tsv
This is a critical file.
Create it with at least 200 real human gene symbols.
Include all major cancer genes, housekeeping genes, pathway genes.
This will be expanded by downloading full list from HGNC.

Include script comment explaining how to download the full list:
# To download complete gene universe:
# wget https://www.genenames.org/cgi-bin/download/custom?...
# Or use: from MyGene.info query all human protein-coding genes

CREATE TEST FILE: tests/test_ora_engine.py

Tests:
- test_contingency_table_construction()
  Known input → verify a, b, c, d values manually

- test_fisher_exact_test()
  Use example with known result:
  [[10, 90], [15, 885]] → expected pvalue ~ 0.0003

- test_bh_correction_applied()
  Multiple p-values → verify FDR corrected values follow BH procedure

- test_pathway_size_filter()
  Pathway with 5 genes (below min_pathway_genes=15) → excluded from results

- test_ora_complete_pipeline()
  Input: 50 known cancer genes
  Expected: at least 1 significant pathway (Cell cycle or DNA repair)

- test_empty_deg_list()
  Input: empty list → return empty DataFrame, no crash

SCIENTIFIC VALIDATION CASES:
Prepare and document this test case:
Input DEGs: TP53, CDKN1A, MDM2, CCND1, RB1 (all p53 pathway genes)
Expected: hsa04115 (p53 signaling pathway) should be enriched
Use this in your viva demonstration.
```

---

## PHASE 7 PROMPT — VISUALIZATION ENGINE

```
[PHASE 7 — VISUALIZATION ENGINE]

Implement all required charts.
Use matplotlib for static publication-quality PNGs.
Use Plotly for interactive HTML charts.
Match the dark color scheme from the UI/UX document.

CREATE FILE: visualization/orf_visualizer.py

1. plot_orf_map(orfs: list[ORF], genome_length: int, output_dir: str) -> str
   Type: Horizontal bar chart showing ORF positions on genome
   
   matplotlib implementation:
   - Figure size: (14, 6)
   - Background: #0F1117
   - X-axis: genome position in bp (0 to genome_length)
   - Y-axis: reading frame labels (+1, +2, +3, -1, -2, -3)
   - Each ORF: horizontal bar from start to stop position
   - Forward strand colors: shades of #4E9FFF
   - Reverse strand colors: shades of #FF5252
   - Annotate longest ORFs with orf_id label
   - X-axis label: "Genome Position (bp)"
   - Title: "ORF Position Map — [genome_name]"
   - DPI: 300 (publication quality)
   - Save as: output_dir/orf_map.png AND orf_map.svg
   
   Return path to PNG file

2. plot_domain_frequency(annotation_results: pd.DataFrame, output_dir: str) -> str
   Type: Vertical bar chart, top 15 Pfam domains
   
   Extract domain counts from annotation results
   Sort by frequency descending
   Bar color: #7C63FF gradient
   X-axis: Domain name (rotate 45°)
   Y-axis: Count of proteins with domain
   DPI: 300
   Return path to PNG

CREATE FILE: visualization/volcano_plotter.py

1. plot_volcano(deg_df: pd.DataFrame, config: dict, output_dir: str) -> tuple[str, str]
   Type: Scatter plot
   
   Plotly implementation:
   
   # Assign colors based on classification
   colors = deg_df['classification'].map({
       'upregulated': '#2DD4A0',
       'downregulated': '#FF5252',
       'not_significant': '#4A4F72'
   })
   
   fig = px.scatter(
       deg_df,
       x='log2FoldChange',
       y='-log10(padj)',  # calculate: -np.log10(df['padj'])
       color='classification',
       hover_data=['gene_symbol', 'log2FoldChange', 'padj'],
       color_discrete_map={...}
   )
   
   # Add threshold lines
   fig.add_vline(x=logfc_up, line_dash='dash', line_color='gray')
   fig.add_vline(x=logfc_down, line_dash='dash', line_color='gray')
   fig.add_hline(y=-np.log10(fdr_cutoff), line_dash='dash', line_color='gray')
   
   # Dark theme
   fig.update_layout(
       plot_bgcolor='#0F1117',
       paper_bgcolor='#1A1D2E',
       font_color='#E8EAF6',
       title='Volcano Plot — Differential Expression'
   )
   
   # Label top 10 most significant genes by name on the plot
   
   # Save both PNG (via kaleido) and HTML
   fig.write_image(f"{output_dir}/volcano_plot.png", width=900, height=600, scale=2)
   fig.write_html(f"{output_dir}/volcano_plot.html")
   
   Return (png_path, html_path)

CREATE FILE: visualization/enrichment_plotter.py

1. plot_kegg_enrichment(enrichment_df: pd.DataFrame, output_dir: str) -> tuple[str, str]
   Type: Horizontal bar chart
   
   Filter to significant results only (padj < 0.05)
   Take top 20 by padj ascending
   
   X-axis: -log10(padj)
   Y-axis: Pathway name (sorted by -log10(padj))
   Bar color: encodes gene_count using colorscale
   
   Plotly implementation with dark theme
   Save PNG and HTML
   Return (png_path, html_path)

2. plot_go_bubble(go_df: pd.DataFrame, output_dir: str) -> tuple[str, str]
   Type: Bubble chart
   
   X-axis: GeneRatio (gene_count / total_degs)
   Y-axis: GO term name (top 20)
   Bubble size: gene_count
   Bubble color: padj (colorscale: low padj = red, high padj = blue)
   
   Plotly with dark theme
   Save PNG and HTML
   Return (png_path, html_path)

CREATE FILE: visualization/report_builder.py

1. generate_html_report(run_data: dict, output_dir: str) -> str
   Use Jinja2 to render complete HTML report
   Include: all tables, charts embedded as base64 images, AI interpretation
   
   Chart embedding:
   import base64
   with open(chart_path, 'rb') as f:
       chart_b64 = base64.b64encode(f.read()).decode()
   html_img_tag = f'<img src="data:image/png;base64,{chart_b64}">'
   
   Report sections (conditionally shown based on workflow):
   - Header (always)
   - Run Summary (always)
   - QC Report (always)
   - ORF Results (Workflow A only)
   - Annotation Results (Workflow A only)
   - DEG Results (Workflow B only)
   - Pathway Enrichment (Workflow B only)
   - Visualizations (always, different charts per workflow)
   - AI Interpretation (if ai.enabled = True)
   - Methods (always — auto-generated from parameters)
   - Metadata (always)

2. generate_methods_text(config: dict, run_metadata: dict) -> str
   Auto-generate scientific methods paragraph:
   "Quality control was performed using PathoScope AI v2.0 with the following
   thresholds: minimum sequence length [X] bp, maximum ambiguity ratio [X]%.
   ORF prediction was performed using six-frame translation with minimum ORF
   length [X] bp and standard genetic code (NCBI table 1). Statistical
   filtering used a log2 fold change threshold of [X] and FDR-corrected
   p-value threshold of [X] (Benjamini-Hochberg correction). Over-representation
   analysis was performed using Fisher's exact test against a background universe
   of [X] human protein-coding genes."

Validation:
- All charts must render without errors on test data
- HTML report must be valid HTML5
- All PNGs must be at 300 DPI minimum
- Embedded charts must display in browser
```

---

## PHASE 8 PROMPT — AI INTERPRETATION AND FINAL REPORT

```
[PHASE 8 — AI INTERPRETATION AND FINAL REPORT GENERATION]

Implement Claude API integration and complete report generation.

CREATE FILE: ai/prompt_templates.py

1. build_workflow_a_prompt(analysis_data: dict) -> str
   Build structured prompt for viral genomics results
   Include:
   - Genome name, length, GC content
   - ORF count, size range, strand distribution
   - Top 10 annotated proteins with descriptions
   - Detected Pfam domains list
   - KEGG functional categories
   
   Prompt instructions:
   "You are a bioinformatics research assistant providing plain-language
   biological interpretation for a university functional genomics project.
   
   The following results were produced by automated computational analysis
   of a viral genome. All values are computationally derived.
   
   [Insert analysis data as structured text]
   
   Provide a 3-4 paragraph biological interpretation:
   1. Describe what the genome structure indicates about this virus
   2. Explain what the identified proteins and domains suggest about viral biology
   3. Discuss functional insights from the KEGG annotations
   4. Summarize the key biological findings
   
   CRITICAL RULES:
   - Do not generate any statistics, p-values, or numerical thresholds
   - Do not invent protein functions not supported by the provided data
   - Write in plain scientific English suitable for undergraduate-level understanding
   - Stay strictly within the scope of the provided results"

2. build_workflow_b_prompt(analysis_data: dict) -> str
   Build structured prompt for functional genomics results
   Include:
   - Total genes analyzed
   - DEG counts (up/down/not significant)
   - Top 10 upregulated genes with fold changes
   - Top 10 downregulated genes with fold changes
   - Top 5 enriched pathways with FDR values
   
   Prompt instructions similar to above, emphasizing:
   - Biological significance of DEG patterns
   - Pathway interpretation
   - Potential biological mechanisms
   - What the results suggest about the biological condition studied

CREATE FILE: ai/claude_interpreter.py

1. generate_interpretation(analysis_data: dict, workflow: str, 
                            config: dict, logger) -> str
   
   from anthropic import Anthropic
   import os
   
   if not config['ai']['enabled']:
       return "AI interpretation disabled. Enable in config.yaml."
   
   api_key = os.getenv('ANTHROPIC_API_KEY')
   if not api_key:
       logger.warning("ANTHROPIC_API_KEY not set. Skipping AI interpretation.")
       return "AI interpretation unavailable: API key not configured."
   
   client = Anthropic()
   
   if workflow == 'A':
       prompt = build_workflow_a_prompt(analysis_data)
   else:
       prompt = build_workflow_b_prompt(analysis_data)
   
   try:
       message = client.messages.create(
           model=config['ai']['model'],
           max_tokens=config['ai']['max_tokens'],
           messages=[{"role": "user", "content": prompt}]
       )
       interpretation = message.content[0].text
       logger.info(f"AI interpretation generated ({len(interpretation)} chars)")
       return interpretation
   
   except Exception as e:
       logger.warning(f"AI interpretation failed: {e}. Continuing without AI.")
       return f"AI interpretation unavailable: {str(e)}"

CREATE Jinja2 HTML Template: visualization/templates/report.html

Create a professional HTML template with:
- Dark theme matching UI/UX document colors
- Responsive CSS (no external dependencies, all inline)
- Navigation: anchor links to each section
- All charts embedded as base64
- All tables: sortable with simple JavaScript
- Footer: PathoScope AI v2.0, Run ID, Timestamp
- Print-friendly CSS media query

UPDATE: main.py

Connect all modules into complete pipeline:

def run_workflow_a(input_obj, config, db_conn, logger, run_dir):
    qc_result = run_fasta_qc(input_obj.file_path, config, logger)
    orfs = predict_orfs(qc_result.passed_sequences, config, logger)
    orfs = translate_all_orfs(orfs, config['orf_prediction']['genetic_code'])
    annotation_results = run_pfam_mapping(orfs, config, logger)
    kegg_results = run_kegg_mapping(annotation_results, config, db_conn, logger)
    charts = generate_workflow_a_charts(orfs, annotation_results, kegg_results, run_dir)
    analysis_data = compile_workflow_a_data(qc_result, orfs, annotation_results, kegg_results)
    ai_text = generate_interpretation(analysis_data, 'A', config, logger)
    report_path = generate_html_report({...}, run_dir)
    return report_path

def run_workflow_b(input_obj, config, db_conn, logger, run_dir):
    if input_obj.detected_format in ['csv_expression', 'tsv_expression']:
        df = load_expression_matrix(input_obj.file_path, logger)
        df = apply_fdr_correction(df, config['statistical']['fdr_method'])
        df = classify_degs(df, config, logger)
        norm_result = normalize_gene_list(df['gene_id'].tolist(), config, db_conn, logger)
    else:
        gene_list = load_gene_list(input_obj.file_path, logger)
        norm_result = normalize_gene_list(gene_list['gene_id'].tolist(), config, db_conn, logger)
    
    deg_lists = get_deg_lists(df) if expression_input else {'all_degs': norm_result.hgnc_symbols}
    kegg_ora = run_ora(deg_lists['all_degs'], load_kegg_gene_sets(config, db_conn, logger), 
                       load_gene_universe(), config, logger)
    go_ora = run_go_enrichment(deg_lists['all_degs'], db_conn, config, logger)
    charts = generate_workflow_b_charts(df, kegg_ora, go_ora, run_dir)
    analysis_data = compile_workflow_b_data(norm_result, df, kegg_ora, go_ora)
    ai_text = generate_interpretation(analysis_data, 'B', config, logger)
    report_path = generate_html_report({...}, run_dir)
    return report_path

Final integration test:
python main.py --input tests/test_data/fasta/ms2_bacteriophage.fasta
Expected: complete HTML report in output/

python main.py --input tests/test_data/gene_lists/expression_matrix.csv
Expected: complete HTML report with volcano plot, enrichment chart, DEG table
```

---

## DEBUG PROMPT (Use When Something Breaks)

```
[DEBUG REQUEST — PathoScope AI v2.0]

I am debugging PathoScope AI v2.0. The following error occurred:

[PASTE EXACT ERROR MESSAGE AND FULL TRACEBACK HERE]

The error occurred when running:
[PASTE THE COMMAND YOU RAN]

The file involved is: [FILENAME]

Known project context:
- Python 3.10+, BioPython, pandas, scipy, statsmodels, pyhmmer
- All parameters come from config.yaml
- Custom exceptions in utils/exceptions.py
- Logger writes to run.log

Please:
1. Explain exactly what caused this error
2. Show the corrected code for the specific function that failed
3. Explain what you changed and why
4. Show how to verify the fix worked

Do NOT rewrite modules I didn't ask about.
Only fix the specific function that is broken.
```

---

## VIVA PREPARATION PROMPT

```
[VIVA PREPARATION — PathoScope AI v2.0]

I am preparing for my viva defense of PathoScope AI v2.0 for my
Functional Genomics course at NCB, Quaid-i-Azam University.

Act as Professor Sir Ghulam Abbas and ask me 15 difficult viva questions about:
1. Why I chose this dual-workflow architecture
2. Why I do NOT use GSEA on viral FASTA files
3. Why ORA is valid for Workflow B but not Workflow A
4. How Benjamini-Hochberg FDR correction works mathematically
5. How Fisher's exact test works for ORA
6. Why I use functional categorization instead of enrichment for small viral genomes
7. How Gene ID normalization works (HGNC/Ensembl/Entrez)
8. What happens if the KEGG API is unavailable (offline mode design)
9. How I ensure reproducibility across different runs
10. What scientific weaknesses remain in the current design

For each question, after I answer, tell me:
- What I got right
- What I missed
- What the ideal answer should include

Be critical and rigorous. This is a serious academic evaluation.
```

---

*End of Master Build Guide and Prompts*
