# PathoScope AI v2.0 — DOCUMENT 2: TRD — REVISED v2.1

---

# ═══════════════════════════════════════════════════════════════
# DOCUMENT 2 — TECHNICAL REQUIREMENTS DOCUMENT (TRD) — REVISED v2.1
# ═══════════════════════════════════════════════════════════════

---

## 1. UPDATED SYSTEM ARCHITECTURE OVERVIEW

```
┌──────────────────────────────────────────────────────────────────────────┐
│                           TIER 1 — INPUT LAYER                            │
│   InputManager → FileDetector → Validator → WorkflowRouter (Safety Gate) │
└──────────────────────────────┬───────────────────────────────────────────┘
                               │
               ┌───────────────┴───────────────┐
               ▼                               ▼
┌──────────────────────────┐       ┌──────────────────────────────────┐
│  TIER 2A — WORKFLOW A    │       │  TIER 2B — WORKFLOW B            │
│  Viral Genomics          │       │  Functional Genomics             │
│                          │       │                                  │
│  QCEngine                │       │  GeneIDNormalizer                │
│  ORFPredictor            │       │  ExpressionLoader                │
│  ProteinTranslator       │       │  StatisticalFilter               │
│  AnnotationEngine        │       │  DEGClassifier                   │
│  PfamMapper              │       │  ORAEngine (Fisher + BH)         │
│  KEGGMapper_A            │       │  GOEnrichment                    │
│  FunctionalCategorizer   │       │  KEGGMapper_B                    │
│  TaxonomyEngine ◄── NEW  │       └──────────────────────────────────┘
│  DashboardGenerator ◄ NEW│               │
└──────────────────────────┘               │
               │                           │
               └───────────────┬───────────┘
                               ▼
┌──────────────────────────────────────────────────────────────────────────┐
│                          TIER 3 — OUTPUT LAYER                            │
│  VisualizationEngine → AIAdapter (Provider-Agnostic) → ReportGenerator  │
│  MetadataWriter → Logger → DatabaseCache                                  │
└──────────────────────────────────────────────────────────────────────────┘
```

---

## 2. UPDATED COMPLETE FOLDER STRUCTURE

```
pathoscope_ai/
│
├── main.py
├── config.yaml
├── requirements.txt
├── environment.yml
├── README.md
├── .env.example
│
├── core/
│   ├── __init__.py
│   ├── input_manager.py
│   ├── file_detector.py
│   ├── validator.py
│   ├── workflow_router.py          ← Safety gate; enforces enrichment rules
│   └── config_loader.py
│
├── workflow_a/
│   ├── __init__.py
│   ├── qc_engine.py
│   ├── sequence_cleaner.py
│   ├── orf_predictor.py
│   ├── protein_translator.py
│   ├── annotation_engine.py        ← Uses DIAMOND/BLAST, NOT custom alignment
│   ├── pfam_mapper.py
│   ├── kegg_mapper_a.py
│   ├── functional_categorizer.py
│   ├── taxonomy_engine.py          ← NEW: NCBI Taxonomy API + cache
│   └── dashboard_generator.py      ← NEW: Genome Intelligence Dashboard
│
├── workflow_b/
│   ├── __init__.py
│   ├── gene_id_normalizer.py
│   ├── expression_loader.py
│   ├── statistical_filter.py
│   ├── deg_classifier.py
│   ├── ora_engine.py
│   ├── go_enrichment.py
│   └── kegg_mapper_b.py
│
├── ai/
│   ├── __init__.py
│   ├── provider_base.py             ← NEW: Abstract AIProvider base class
│   ├── gemini_provider.py           ← NEW: Google Gemini adapter
│   ├── openai_provider.py           ← NEW: OpenAI GPT-4 adapter
│   ├── claude_provider.py           ← NEW: Anthropic Claude adapter
│   ├── provider_factory.py          ← NEW: Returns correct provider from config
│   └── prompt_templates.py
│
├── visualization/
│   ├── __init__.py
│   ├── orf_visualizer.py
│   ├── volcano_plotter.py
│   ├── enrichment_plotter.py
│   ├── domain_plotter.py
│   ├── taxonomy_tree_viz.py         ← NEW: Interactive taxonomy tree
│   ├── genome_map_viz.py            ← NEW: NCBI-style genome viewer
│   ├── similarity_dashboard.py      ← NEW: Virus similarity comparison
│   ├── executive_panel.py           ← NEW: Metric cards for dashboard
│   └── report_builder.py
│
├── databases/
│   ├── pfam_viral_subset/           ← Viral Pfam HMM profiles (~50 MB)
│   ├── pathoscope.sqlite            ← All caches and run tracking
│   ├── ensembl_hgnc_map.tsv
│   ├── entrez_hgnc_map.tsv
│   ├── gene_universe.tsv
│   ├── taxonomy_cache.sqlite        ← NEW: NCBI taxonomy lineage cache
│   └── virus_reference_db.sqlite   ← NEW: Known virus reference metadata
│
├── tests/
│   ├── __init__.py
│   ├── test_data/
│   │   ├── fasta/
│   │   │   ├── ms2_bacteriophage.fasta      (FASTA-01)
│   │   │   ├── phix174.fasta                (FASTA-02)
│   │   │   ├── sars_cov2.fasta              (FASTA-03)
│   │   │   ├── influenza_a_ha.fasta         (FASTA-04)
│   │   │   ├── invalid_chars.fasta          (ERR-01)
│   │   │   ├── high_ambiguity.fasta         (ERR-02)
│   │   │   └── duplicates.fasta             (ERR-03)
│   │   ├── fastq/
│   │   │   ├── sars_cov2_reads.fastq        (FASTQ-01)
│   │   │   ├── phix_control.fastq           (FASTQ-02)
│   │   │   └── corrupted.fastq              (FASTQ-03)
│   │   └── gene_lists/
│   │       ├── hgnc_symbols.txt             (FG-01)
│   │       ├── ensembl_ids.txt              (FG-02)
│   │       ├── entrez_ids.txt               (FG-03)
│   │       ├── mixed_ids.csv                (FG-04)
│   │       ├── expression_matrix.csv        (FG-05)
│   │       └── expression_missing_vals.csv  (FG-06)
│   │
│   ├── test_input_manager.py
│   ├── test_workflow_router.py
│   ├── test_qc_engine.py
│   ├── test_orf_predictor.py
│   ├── test_protein_translator.py
│   ├── test_gene_normalizer.py
│   ├── test_statistical_filter.py
│   ├── test_deg_classifier.py
│   ├── test_ora_engine.py
│   ├── test_taxonomy_engine.py      ← NEW
│   ├── test_dashboard_generator.py  ← NEW
│   └── test_report_generator.py
│
└── output/
    └── run_YYYYMMDD_HHMMSS/
        ├── metadata.json
        ├── run.log
        ├── config_snapshot.yaml
        ├── qc/
        ├── sequences/
        ├── annotation/
        ├── taxonomy/               ← NEW
        │   ├── taxonomy_lineage.json
        │   ├── taxonomy_tree.html
        │   └── taxonomy_report.csv
        ├── dashboard/              ← NEW
        │   ├── executive_panel.html
        │   ├── genome_map.html
        │   ├── domain_treemap.html
        │   └── similarity_dashboard.html
        ├── deg/
        ├── enrichment/
        ├── visualizations/
        ├── ai/
        └── reports/
```

---

## 3. TECHNOLOGY STACK

### Core Bioinformatics
| Library | Version | Purpose |
|---|---|---|
| BioPython | >=1.81 | FASTA/FASTQ parsing, ORF detection, sequence ops |
| pyhmmer | >=0.10 | Pfam HMM domain search (no binary dependency) |
| pandas | >=2.0 | Data handling, tables, results |
| numpy | >=1.24 | Numerical operations |
| scipy | >=1.11 | Fisher's exact test, statistical functions |
| statsmodels | >=0.14 | Benjamini-Hochberg FDR correction |
| requests | >=2.31 | API calls (KEGG, NCBI Taxonomy, Ensembl) |

**Alignment Note:** PathoScope AI uses DIAMOND (via subprocess) or BLAST (via BioPython NCBIXML) for annotation. No custom alignment algorithm is implemented. Smith-Waterman from scratch is explicitly excluded from the project scope.

### Visualization
| Library | Purpose |
|---|---|
| matplotlib | ORF maps, domain frequency charts (static PNG) |
| plotly | Volcano plots, enrichment charts, taxonomy tree, genome map (interactive HTML) |
| seaborn | Heatmaps, correlation charts |
| networkx | Taxonomy tree graph construction |
| pyvis | Interactive network/tree HTML rendering |

### AI Provider Layer (NEW — Provider-Agnostic)
| Component | Purpose |
|---|---|
| anthropic | Claude API client (ClaudeProvider adapter) |
| google-generativeai | Gemini API client (GeminiProvider adapter) |
| openai | OpenAI GPT-4 client (OpenAIProvider adapter) |

### Database and Caching
| Library | Purpose |
|---|---|
| sqlite3 | Caching, run tracking, taxonomy cache |
| PyYAML | config.yaml |
| python-dotenv | .env API key loading |

### Reporting
| Library | Purpose |
|---|---|
| jinja2 | HTML report templating |
| WeasyPrint | PDF generation |

### Testing
| Library | Purpose |
|---|---|
| pytest | Test runner |
| pytest-cov | Coverage reporting |
| responses | Mock HTTP requests in tests |
| pytest-mock | General mocking |

---

## 4. MODULE SPECIFICATIONS

### Module 4.1 — WorkflowRouter (core/workflow_router.py)

**This is a scientific validity enforcement component, not just a routing utility.**

```python
WORKFLOW_MAP = {
    'fasta':            'A',
    'fastq':            'A',
    'csv_expression':   'B',
    'tsv_expression':   'B',
    'txt_genes':        'B',
    'csv_genes':        'B',
}

WORKFLOW_A_FORBIDDEN = ['ORA', 'DEG', 'GSEA', 'FDR_ENRICHMENT']
WORKFLOW_B_FORBIDDEN = ['ORF_PREDICTION', 'PFAM_MAPPING', 'TAXONOMY']
```

**Functions:**
```
determine_workflow(detected_format: str, override: str = None) -> str
validate_workflow_capability(workflow: str, analysis_type: str) -> None
    # Raises WorkflowViolationError if analysis_type is forbidden for workflow
log_routing_decision(format, workflow, logger)
```

**WorkflowViolationError** — custom exception that MUST be raised and logged if any module attempts to apply enrichment analysis to Workflow A data. This prevents silent scientific errors.

---

### Module 4.2 — QCEngine (workflow_a/qc_engine.py)

Same specification as previous TRD. Key thresholds:
- min_sequence_length: 50 bp (configurable)
- max_ambiguity_ratio: 0.10 (configurable)
- min_fastq_quality: 20 (configurable)

---

### Module 4.3 — AnnotationEngine (workflow_a/annotation_engine.py)

**Tool selection policy:**
```
Priority 1: DIAMOND (subprocess) — fastest, preferred for large datasets
Priority 2: BioPython BLAST wrapper (via NCBIBlast+ local install if available)
Priority 3: NCBI BLAST+ web API (online fallback, rate-limited)

NO CUSTOM ALIGNMENT CODE.
```

**Functions:**
```
run_annotation(protein_fasta_path: str, config: dict, logger) -> list[AnnotationResult]
run_diamond(protein_fasta: str, db_path: str, config: dict) -> list[Hit]
run_blast_api(protein_fasta: str, config: dict) -> list[Hit]
parse_annotation_results(raw_results: list) -> list[AnnotationResult]
filter_by_evalue(results: list, evalue: float) -> list[AnnotationResult]
filter_by_identity(results: list, min_identity: float) -> list[AnnotationResult]
```

---

### Module 4.4 — TaxonomyEngine (workflow_a/taxonomy_engine.py) ← NEW

**Purpose:** After annotation identifies the closest reference genome, retrieve the complete ICTV/NCBI taxonomic lineage, cache it locally, and prepare a structured taxonomy report for the dashboard.

**Data Sources (priority):**
1. Local taxonomy_cache.sqlite (fastest, offline)
2. NCBI Taxonomy API (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/)
3. NCBI Datasets API (https://api.ncbi.nlm.nih.gov/datasets/v2/)

**Key Functions:**
```python
fetch_taxonomy_lineage(taxon_id: int, config: dict, db_conn) -> TaxonomyLineage
    # Returns full lineage from Realm to Species

classify_virus_from_annotation(annotation_results: list, config: dict, db_conn) -> VirusClassification
    # 1. Extract best annotation hit accession
    # 2. Look up taxon ID via NCBI E-utilities
    # 3. Retrieve full lineage

get_taxon_id(accession: str, config: dict) -> int
    # Query: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi
    #        db=taxonomy&term={accession}

get_lineage(taxon_id: int, config: dict) -> dict
    # Query: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
    #        db=taxonomy&id={taxon_id}&retmode=xml

cache_taxonomy(lineage: dict, taxon_id: int, db_conn)
get_cached_taxonomy(taxon_id: int, db_conn) -> dict or None

generate_taxonomy_report(classification: VirusClassification, output_dir: str) -> str
    # Writes taxonomy_lineage.json + taxonomy_report.csv
    # Returns path to JSON file

find_similar_viruses(classification: VirusClassification, db_conn) -> list[SimilarVirus]
    # Queries virus_reference_db.sqlite for viruses in the same family
    # Returns list with name, similarity %, coverage %, genome size
```

**TaxonomyLineage dataclass:**
```python
@dataclass
class TaxonomyLineage:
    taxon_id: int
    organism_name: str
    realm: str
    kingdom: str
    phylum: str
    subphylum: str
    class_: str
    order: str
    family: str
    subfamily: str
    genus: str
    species: str
    host: str           # from NCBI metadata where available
    genome_type: str    # dsDNA, ssRNA, etc.
    retrieval_method: str  # 'ncbi_api' or 'local_cache'
```

**Offline Fallback Strategy:**
The file `databases/taxonomy_cache.sqlite` must be pre-populated with lineage data for the four validation genomes (MS2, PhiX174, SARS-CoV-2, Influenza A) so that viva demonstrations work without internet.

---

### Module 4.5 — DashboardGenerator (workflow_a/dashboard_generator.py) ← NEW

**Purpose:** Orchestrate all Taxonomy Dashboard visualizations and generate the master dashboard HTML page, inspired by NCBI Virus Dashboard (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/dashboard).

**Functions:**
```python
generate_full_dashboard(
    qc_result: QCResult,
    orfs: list[ORF],
    annotation_results: list[AnnotationResult],
    taxonomy: TaxonomyLineage,
    similar_viruses: list[SimilarVirus],
    pfam_results: list[PfamResult],
    kegg_results: list[KEGGResult],
    output_dir: str,
    config: dict
) -> DashboardPaths

generate_executive_panel(taxonomy, orfs, annotation_results, pfam_results) -> str
    # HTML fragment: metric cards for key statistics

generate_taxonomy_tree(lineage: TaxonomyLineage, output_dir: str) -> str
    # Interactive tree using networkx + pyvis
    # Returns path to HTML file

generate_genome_map(orfs: list[ORF], genome_length: int, annotation_results, output_dir) -> str
    # NCBI-style horizontal genome viewer
    # Returns path to HTML file

generate_domain_treemap(pfam_results: list[PfamResult], output_dir: str) -> str
    # Plotly treemap of Pfam domain distribution
    # Returns path to HTML file

generate_similarity_dashboard(similar_viruses: list[SimilarVirus], output_dir: str) -> str
    # Horizontal bar chart: virus name vs similarity %
    # Returns path to HTML file

generate_pathway_mapping_dashboard(kegg_results: list[KEGGResult], output_dir: str) -> str
    # Sunburst / bar chart of KEGG category distribution
    # Returns path to HTML file

assemble_dashboard_html(all_panels: dict, taxonomy: TaxonomyLineage, output_dir: str) -> str
    # Renders complete taxonomy_dashboard.html using Jinja2
    # Embeds all panel HTML as iframe sections or inline
    # Returns path to final dashboard HTML
```

**DashboardPaths dataclass:**
```python
@dataclass
class DashboardPaths:
    executive_panel_html: str
    taxonomy_tree_html: str
    genome_map_html: str
    domain_treemap_html: str
    similarity_dashboard_html: str
    pathway_mapping_html: str
    complete_dashboard_html: str    # Master assembled page
```

---

### Module 4.6 — Provider-Agnostic AI Layer ← NEW

**File: ai/provider_base.py**
```python
from abc import ABC, abstractmethod

class AIProvider(ABC):
    """Abstract base class for all AI interpretation providers."""
    
    @abstractmethod
    def generate_interpretation(self, prompt: str, max_tokens: int) -> str:
        pass
    
    @abstractmethod
    def is_available(self) -> bool:
        """Check if provider API key is set and provider is reachable."""
        pass
    
    @property
    @abstractmethod
    def provider_name(self) -> str:
        pass
```

**File: ai/gemini_provider.py**
```python
import google.generativeai as genai

class GeminiProvider(AIProvider):
    def __init__(self, api_key: str, model: str = "gemini-1.5-pro"):
        genai.configure(api_key=api_key)
        self.model = genai.GenerativeModel(model)
    
    def generate_interpretation(self, prompt: str, max_tokens: int) -> str:
        response = self.model.generate_content(prompt)
        return response.text
    
    def is_available(self) -> bool:
        return bool(os.getenv('GEMINI_API_KEY'))
    
    @property
    def provider_name(self) -> str:
        return "Google Gemini"
```

**File: ai/claude_provider.py**
```python
from anthropic import Anthropic

class ClaudeProvider(AIProvider):
    def __init__(self, api_key: str, model: str = "claude-sonnet-4-20250514"):
        self.client = Anthropic(api_key=api_key)
        self.model = model
    
    def generate_interpretation(self, prompt: str, max_tokens: int) -> str:
        message = self.client.messages.create(
            model=self.model,
            max_tokens=max_tokens,
            messages=[{"role": "user", "content": prompt}]
        )
        return message.content[0].text
    
    def is_available(self) -> bool:
        return bool(os.getenv('ANTHROPIC_API_KEY'))
    
    @property
    def provider_name(self) -> str:
        return "Anthropic Claude"
```

**File: ai/openai_provider.py**
```python
from openai import OpenAI

class OpenAIProvider(AIProvider):
    def __init__(self, api_key: str, model: str = "gpt-4o"):
        self.client = OpenAI(api_key=api_key)
        self.model = model
    
    def generate_interpretation(self, prompt: str, max_tokens: int) -> str:
        response = self.client.chat.completions.create(
            model=self.model,
            max_tokens=max_tokens,
            messages=[{"role": "user", "content": prompt}]
        )
        return response.choices[0].message.content
    
    def is_available(self) -> bool:
        return bool(os.getenv('OPENAI_API_KEY'))
    
    @property
    def provider_name(self) -> str:
        return "OpenAI GPT-4"
```

**File: ai/provider_factory.py**
```python
def get_ai_provider(config: dict) -> AIProvider:
    """Returns the configured AI provider. Falls back gracefully if unavailable."""
    
    provider_name = config['ai'].get('provider', 'claude').lower()
    
    providers = {
        'claude':  (ClaudeProvider, 'ANTHROPIC_API_KEY', config['ai'].get('claude_model', 'claude-sonnet-4-20250514')),
        'gemini':  (GeminiProvider, 'GEMINI_API_KEY',    config['ai'].get('gemini_model', 'gemini-1.5-pro')),
        'openai':  (OpenAIProvider, 'OPENAI_API_KEY',    config['ai'].get('openai_model', 'gpt-4o')),
    }
    
    if provider_name not in providers:
        raise ConfigError(f"Unknown AI provider: {provider_name}. Supported: {list(providers.keys())}")
    
    ProviderClass, env_key, model = providers[provider_name]
    api_key = os.getenv(env_key)
    
    if not api_key:
        logger.warning(f"AI provider '{provider_name}' requires {env_key} env var. AI section will be skipped.")
        return NullProvider()  # Returns empty interpretation gracefully
    
    return ProviderClass(api_key=api_key, model=model)
```

---

## 5. UPDATED DATABASE SCHEMA ADDITIONS

### New Table: taxonomy_cache
```sql
CREATE TABLE taxonomy_cache (
    taxon_id        INTEGER PRIMARY KEY,
    organism_name   TEXT NOT NULL,
    realm           TEXT,
    kingdom         TEXT,
    phylum          TEXT,
    subphylum       TEXT,
    class_name      TEXT,
    order_name      TEXT,
    family          TEXT,
    subfamily       TEXT,
    genus           TEXT,
    species         TEXT,
    host_organism   TEXT,
    genome_type     TEXT,
    raw_ncbi_xml    TEXT,
    retrieved_at    TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    expires_at      TIMESTAMP NOT NULL
);
CREATE INDEX idx_taxonomy_family ON taxonomy_cache(family);
CREATE INDEX idx_taxonomy_genus  ON taxonomy_cache(genus);
```

### New Table: virus_reference_db
```sql
CREATE TABLE virus_reference_db (
    accession       TEXT PRIMARY KEY,
    virus_name      TEXT NOT NULL,
    family          TEXT,
    genus           TEXT,
    genome_length   INTEGER,
    gc_content      REAL,
    host            TEXT,
    genome_type     TEXT,
    ncbi_taxon_id   INTEGER,
    reference_seq   BOOLEAN DEFAULT FALSE
);
```

### New Table: dashboard_runs
```sql
CREATE TABLE dashboard_runs (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          TEXT NOT NULL REFERENCES runs(run_id),
    taxon_id        INTEGER,
    organism_name   TEXT,
    taxonomy_method TEXT,   -- 'ncbi_api' or 'local_cache'
    dashboard_path  TEXT,   -- path to complete dashboard HTML
    panels_generated TEXT,  -- JSON array of panel names generated
    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

---

## 6. UPDATED CONFIGURATION SYSTEM (config.yaml)

```yaml
# PathoScope AI v2.0 Configuration — Revised v2.1

project:
  name: "PathoScope AI"
  version: "2.1.0"
  output_dir: "./output"
  log_level: "INFO"

workflow:
  mode: "auto"    # auto | workflow_a | workflow_b

quality_control:
  min_sequence_length: 50
  max_ambiguity_ratio: 0.10
  min_fastq_quality: 20
  remove_duplicates: true
  strict_mode: false    # if true: reject file on any invalid sequence

orf_prediction:
  min_orf_length: 100
  genetic_code: 1
  both_strands: true
  remove_nested: false

annotation:
  # IMPORTANT: No custom alignment. Use validated tools only.
  tool: "diamond"         # Options: diamond | blast_local | blast_api
  diamond_db: "./databases/diamond/viral_proteins.dmnd"
  blast_db: "./databases/blast/viral_nr"
  evalue: 0.001
  identity: 30.0
  coverage: 50.0
  pfam_evalue: 0.01
  pfam_db: "./databases/pfam_viral_subset/Pfam-A-viral.hmm"

statistical:
  logfc_up_threshold: 1.0
  logfc_down_threshold: -1.0
  pvalue_cutoff: 0.05
  fdr_cutoff: 0.05
  fdr_method: "benjamini-hochberg"
  min_pathway_genes: 15
  max_pathway_genes: 500

taxonomy:
  enabled: true
  ncbi_taxonomy_base_url: "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  ncbi_datasets_url: "https://api.ncbi.nlm.nih.gov/datasets/v2/"
  taxonomy_db: "./databases/taxonomy_cache.sqlite"
  ttl_days: 90              # taxonomy is very stable
  generate_dashboard: true  # set false to skip dashboard generation

api:
  kegg_base_url: "https://rest.kegg.jp"
  ensembl_base_url: "https://rest.ensembl.org"
  ncbi_base_url: "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  request_timeout: 30
  retry_attempts: 3
  retry_delay: 2
  offline_mode: false

cache:
  sqlite_path: "./databases/pathoscope.sqlite"
  kegg_ttl_days: 7
  gene_id_ttl_days: 30

# AI Provider Configuration — Choose one provider
ai:
  enabled: true
  provider: "gemini"      # Options: gemini | openai | claude
  max_tokens: 2000

  # Provider-specific models (only active provider model is used)
  gemini_model: "gemini-1.5-pro"
  openai_model: "gpt-4o"
  claude_model: "claude-sonnet-4-20250514"

  # AI boundary enforcement
  allowed_tasks:
    - "biological_interpretation"
    - "pathway_explanation"
    - "protein_function_description"
  forbidden_tasks:
    - "orf_prediction"
    - "statistics_generation"
    - "taxonomy_classification"
    - "pathway_enrichment"
    - "gene_annotation"
```

---

## 7. UPDATED .env.example

```
# PathoScope AI — API Keys
# Only the key for your chosen provider (config.yaml ai.provider) is required

ANTHROPIC_API_KEY=your_claude_api_key_here
GEMINI_API_KEY=your_gemini_api_key_here
OPENAI_API_KEY=your_openai_api_key_here

NCBI_API_KEY=your_ncbi_api_key_here_optional
# NCBI key is optional but increases rate limit from 3 to 10 req/sec
```

---

## 8. TAXONOMY DASHBOARD — COMPONENT SPECIFICATIONS

### Component 8.1 — Executive Panel
Displayed at top of dashboard. Six metric cards in a responsive grid.

```
┌──────────────────────────────────────────────────────────────────┐
│  VIRUS GENOME INTELLIGENCE DASHBOARD                             │
│  Run: run_20260531_102345  │  SARS-CoV-2  │  Betacoronavirus     │
├──────────┬──────────┬──────────┬──────────┬──────────┬──────────┤
│ Genome   │ GC       │ ORFs     │ Annotated│ Pfam     │ Top      │
│ Length   │ Content  │ Found    │ Proteins │ Domains  │ Pathway  │
│          │          │          │          │          │ Category │
│ 29,903bp │ 37.97%   │   12     │ 10/12    │   8      │Replication│
└──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘
│ Closest Reference: NC_045512.2  │  Family: Coronaviridae         │
│ Genus: Betacoronavirus          │  Host: Homo sapiens           │
└──────────────────────────────────────────────────────────────────┘
```

### Component 8.2 — Taxonomy Tree
Interactive phylogenetic-style tree showing lineage from Realm to Species.

**Technology:** networkx (graph) + pyvis (HTML rendering)

```
Riboviria
  └── Orthornavirae
        └── Pisuviricota
              └── Pisoniviricetes
                    └── Nidovirales
                          └── Coronaviridae
                                └── Betacoronavirus
                                      └── ● SARS-CoV-2  ← highlighted
```

Each node: clickable, shows NCBI taxon ID and description on hover.
Highlighted path: gold color from root to query virus.

### Component 8.3 — Genome Map
NCBI-style horizontal genome viewer showing ORF positions.

**Technology:** Plotly (interactive horizontal bar chart)

```
5'━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━3'
   [━━━━━━━━━━━━━━━━━━━━━━━━━ ORF1ab (pp1ab) ━━━━━━━━━━━━━━━━━━━]
                                               [━━ S (Spike) ━━]
                                                              [E][M]
                                                                 [━ N ━]
   0    5000   10000   15000   20000   25000   29903 bp
```

Hover on each ORF: shows ORF ID, protein name, start, stop, length, annotation.
Color coding by protein function family.

### Component 8.4 — Protein Domain Treemap
Plotly treemap showing hierarchical distribution of Pfam domain families.

```
┌─────────────────────────────────────────────────────────────┐
│  Protein Domain Distribution (Pfam)                         │
│                                                             │
│  ┌──────────────────────┬──────────────────────────────┐   │
│  │                      │                              │   │
│  │  RNA-dependent       │   Spike protein binding      │   │
│  │  RNA polymerase      │   domain (PF09408)           │   │
│  │  (PF00680)           │                              │   │
│  │                      ├──────────┬───────────────────┤   │
│  │                      │ Nucleocap│   Envelope        │   │
│  └──────────────────────┘ (PF00937)│   protein         │   │
│                           └─────────┴───────────────────┘  │
└─────────────────────────────────────────────────────────────┘
```

### Component 8.5 — Virus Similarity Dashboard
Horizontal bar chart showing top 10 most similar known viruses.

```
Virus                          Similarity%   Coverage%   Accession
SARS-CoV-2 Ref (NC_045512.2)  ███████████  99.8%       100%
SARS-CoV (NC_004718.3)        ████████░░░  79.4%        95%
Bat coronavirus RaTG13         ███████░░░░  96.2%        91%
MERS-CoV (NC_019843.3)        ████░░░░░░░  51.2%        78%
...
```

### Component 8.6 — Pathway Mapping Dashboard
Sunburst chart: KEGG functional categories → subcategories → individual KO entries.

**Technology:** Plotly Sunburst

### Component 8.7 — Host Distribution (when metadata available)
Horizontal bar or pie chart of host species from virus_reference_db.sqlite for the same family.

### Component 8.8 — Download Center
```
Downloads available:
┌──────────────────────────────────────────────────────────────┐
│ 📄 Final Report (HTML)        final_report.html   [Download] │
│ 📊 Taxonomy Dashboard (HTML)  taxonomy_dashboard.html        │
│ 📋 Taxonomy Lineage (JSON)    taxonomy_lineage.json          │
│ 📊 Annotation Results (CSV)   annotation_results.csv         │
│ 📊 ORF Results (CSV)          orf_results.csv                │
│ 🧬 ORF Sequences (FASTA)      predicted_orfs.fasta           │
│ 🧬 Protein Sequences (FASTA)  predicted_proteins.fasta       │
│ 📋 GFF3 Annotation            annotation.gff3               │
│ 📋 Run Metadata (JSON)        metadata.json                  │
│ 📦 Download All (ZIP)                           [Download All]│
└──────────────────────────────────────────────────────────────┘
```

---

## 9. UPDATED TEST MODULE SPECIFICATIONS

### test_taxonomy_engine.py
```python
test_detect_taxon_id_from_accession()
    # Input: NC_045512 → Expected taxon_id: 2697049

test_fetch_sars_cov2_lineage()
    # Expected: family='Coronaviridae', genus='Betacoronavirus'

test_offline_fallback_uses_cache()
    # Pre-populate cache, set offline_mode=True
    # Expected: lineage returned from cache, no HTTP calls

test_lineage_cached_after_api_call()
    # After API call, verify record exists in taxonomy_cache table

test_similar_viruses_found_for_coronaviridae()
    # Input: Coronaviridae family
    # Expected: at least 3 similar viruses from reference_db

test_taxonomy_report_csv_generated()
test_invalid_accession_handled_gracefully()
    # Input: made-up accession
    # Expected: VirusClassification with organism_name='Unknown', no crash
```

### test_dashboard_generator.py
```python
test_executive_panel_html_generated()
test_taxonomy_tree_html_generated()
test_genome_map_html_generated()
test_domain_treemap_html_generated()
test_similarity_dashboard_html_generated()
test_full_dashboard_assembles_without_error()
test_dashboard_offline_mode_works()
    # All components generated without internet using cached data
```

---

## 10. ENRICHMENT RULES — FORMAL DOCUMENTATION

The following rules are enforced programmatically by WorkflowRouter and documented here for scientific auditability:

```
RULE ENR-01: ORA_VIRAL_FORBIDDEN
  Condition: workflow == 'A' AND analysis_type == 'ORA'
  Action: raise WorkflowViolationError
  Scientific Reason: Viral proteomes (<100 proteins) lack statistical power
                     for Over-Representation Analysis. Background universe
                     undefined for non-host genomes.

RULE ENR-02: GSEA_ALWAYS_FORBIDDEN
  Condition: analysis_type == 'GSEA'
  Action: raise WorkflowViolationError
  Scientific Reason: GSEA requires continuous ranked expression data from
                     RNA-seq differential expression analysis. Neither FASTA
                     sequences nor simple gene lists provide this data.

RULE ENR-03: DEG_VIRAL_FORBIDDEN
  Condition: workflow == 'A' AND analysis_type == 'DEG_CLASSIFICATION'
  Action: raise WorkflowViolationError
  Scientific Reason: DEG classification requires paired expression values
                     (logFC, pvalue) from RNA-seq. FASTA sequences do not
                     contain expression measurements.

RULE ENR-04: FUNCTIONAL_CATEGORIZATION_VIRAL_ONLY
  Condition: workflow == 'A'
  Action: use FunctionalCategorizer (NOT enrichment statistics)
  Scientific Reason: For small viral proteomes, reporting the distribution
                     of KEGG functional categories is scientifically valid
                     and interpretable. It describes what categories of
                     proteins are present, without claiming statistical
                     over-representation.

RULE ENR-05: ORA_BACKGROUND_REQUIRED
  Condition: workflow == 'B' AND analysis_type == 'ORA'
  Pre-condition: universe = load_gene_universe() with len(universe) >= 1000
  Action: proceed with Fisher's exact test
  Scientific Reason: ORA is valid when a proper background gene universe
                     is defined. We use the full human protein-coding gene
                     set (~19,000 genes) as background.
```

---

*End of Document 2 — TRD v2.1*
