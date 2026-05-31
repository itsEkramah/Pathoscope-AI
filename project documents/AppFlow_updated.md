# PathoScope AI v2.0 — DOCUMENT 3: APP FLOW — REVISED v2.1

---

# ═══════════════════════════════════════════════════════════════
# DOCUMENT 3 — APP FLOW DOCUMENT — REVISED v2.1
# ═══════════════════════════════════════════════════════════════

**Changes in v2.1:**
- Executive Dashboard replaces static Home Page as primary landing view
- Taxonomy tab added to results navigation
- Dashboard tab added as major results view
- AI provider shown in UI (Gemini/OpenAI/Claude indicator)
- Download Center expanded (GFF3, taxonomy JSON)
- Error states added for taxonomy classification failures
- Enrichment guard state documented (shows explanation when ORA requested on Workflow A)

---

## 1. SYSTEM INTERACTION OVERVIEW

PathoScope AI v2.0 supports two interaction modes:

**Mode 1 — CLI (Command Line Interface)**
Primary mode. User runs via terminal. Best for automated scripting and headless servers.

**Mode 2 — Web Dashboard (Flask)**
Secondary mode. User accesses via browser. Primary mode for viva demonstration and student use.

Both modes use identical backend pipeline modules. The Flask server is a thin presentation layer over the same Python functions.

---

## 2. CLI USER FLOW (Updated v2.1)

### Flow 2.1 — Standard Viral Genome Analysis (Workflow A)
```
$ python main.py --input data/sars_cov2.fasta

╔══════════════════════════════════════════╗
║         PathoScope AI v2.1               ║
║   Viral Functional Genomics Platform     ║
║   NCB, Quaid-i-Azam University           ║
╚══════════════════════════════════════════╝

[DETECT]  Format: FASTA
[ROUTE]   Workflow A (Viral Genomics)  ← Safety gate confirmed: ORA/DEG disabled
[RUN ID]  run_20260531_102345
[OUTPUT]  ./output/run_20260531_102345/

Pipeline:
  ✓ [QC]           1 sequence loaded, 0 removed
  ✓ [ORF]          12 ORFs predicted
  ✓ [TRANSLATE]    12 protein sequences
  ⟳ [ANNOTATE]    10/12 annotated (DIAMOND)...
  ✓ [PFAM]         8 Pfam domains mapped
  ✓ [KEGG]         4 KEGG categories
  ⟳ [TAXONOMY]    Fetching NCBI lineage for NC_045512.2...
  ✓ [TAXONOMY]    SARS-CoV-2 → Betacoronavirus → Coronaviridae
  ✓ [DASHBOARD]   Taxonomy dashboard generated
  ✓ [AI]          Interpretation generated (Provider: Gemini)
  ✓ [REPORT]      HTML report generated

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  ✓ Analysis complete in 47.2 seconds
  Run ID : run_20260531_102345
  Report : ./output/run_20260531_102345/reports/final_report.html
  Dashboard: ./output/run_20260531_102345/dashboard/taxonomy_dashboard.html
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Flow 2.2 — Expression Matrix Analysis (Workflow B)
```
$ python main.py --input data/expression_matrix.csv

[DETECT]  Format: CSV Expression Matrix (logFC + pvalue columns detected)
[ROUTE]   Workflow B (Functional Genomics)  ← Taxonomy/ORF/Pfam disabled
[NOTE]    ORA enrichment will use human gene universe (19,194 genes)

Pipeline:
  ✓ [LOAD]          2,000 genes loaded
  ✓ [NORMALIZE]     1,987 genes normalized to HGNC (13 unmapped, see log)
  ✓ [FDR]           BH correction applied
  ✓ [CLASSIFY]      UP: 198, DOWN: 204, NS: 1,598
  ✓ [ORA KEGG]      12 significant pathways (FDR < 0.05)
  ✓ [ORA GO-BP]     34 significant GO terms
  ✓ [VOLCANO]       Volcano plot generated
  ✓ [ENRICHMENT]    Enrichment charts generated
  ✓ [AI]            Interpretation generated (Provider: Gemini)
  ✓ [REPORT]        HTML report generated
```

### Flow 2.3 — Offline Mode
```
$ python main.py --input data/sars_cov2.fasta --offline

[MODE] Offline mode enabled
[NOTE] All API calls disabled. Using cached databases only.
[TAXONOMY] Found cached lineage for SARS-CoV-2 (taxon_id: 2697049)
...
[NOTE] AI interpretation skipped: provider unavailable in offline mode
```

### Flow 2.4 — Switch AI Provider
No code change required. Only config.yaml change:
```yaml
ai:
  provider: "gemini"   # change to: openai | claude
```
```
[AI] Provider: Gemini → OpenAI GPT-4  (config change only)
```

### Flow 2.5 — Error Flow: Invalid FASTA
```
$ python main.py --input data/invalid_chars.fasta

[DETECT]  Format: FASTA
[ROUTE]   Workflow A
[QC]      ✗ Sequence 'seq_001': invalid characters found: X, Z, J
[QC]      ✗ Sequence 'seq_002': ambiguity ratio 0.23 exceeds threshold 0.10
[RESULT]  0 sequences passed QC
[ERROR]   InputValidationError: No valid sequences remain after QC.
          Use --strict-mode=False to skip invalid sequences instead of stopping.
          Log: ./output/run_20260531_104512/run.log
```

---

## 3. WEB DASHBOARD USER FLOW (Updated v2.1)

### Page 3.0 — EXECUTIVE DASHBOARD (New Default Landing Page)

**Purpose:** Replaces the old static home page. On first visit, shows platform description and upload prompt. After any analysis, shows the most recent run's results directly as an executive dashboard. This is inspired by NCBI Virus Dashboard's always-informative landing view.

**Layout — First Visit (No Prior Analysis):**
```
┌────────────────────────────────────────────────────────────────────┐
│  🔬 PathoScope AI v2.1    [Upload] [Docs] [Settings] [About]      │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│   Viral Functional Genomics & Genome Intelligence Platform         │
│   NCB, Quaid-i-Azam University                                    │
│                                                                    │
│   ┌──────────────────────────────────────────────────────────┐   │
│   │  ↑ Drop a file here to begin analysis                    │   │
│   │  FASTA  FASTQ  FASTQ.GZ  CSV  TSV  TXT                  │   │
│   │  [Browse Files]   or   [Load Example: SARS-CoV-2]       │   │
│   └──────────────────────────────────────────────────────────┘   │
│                                                                    │
│   ── Or choose a workflow ─────────────────────────────────────   │
│   [🧬 Workflow A — Viral Genomics]  [📊 Workflow B — Func. Genomics]│
│                                                                    │
│   ── Supported Analysis Types ────────────────────────────────   │
│   Workflow A: QC → ORF → Pfam → KEGG → Taxonomy Dashboard       │
│   Workflow B: Gene ID Normalization → DEG → ORA → Enrichment     │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

**Layout — After Analysis (Workflow A — shows NCBI-inspired dashboard):**
```
┌────────────────────────────────────────────────────────────────────┐
│  🔬 PathoScope AI v2.1    [New Analysis] [Docs] [Settings]        │
│  Run: run_20260531_102345 │ sars_cov2.fasta │ 47.2s │ AI: Gemini  │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  VIRUS GENOME INTELLIGENCE DASHBOARD                               │
│  ┌────────────────────────────────────────────────────────────┐   │
│  │ Organism:  SARS-CoV-2                                     │   │
│  │ Taxonomy:  Coronaviridae → Betacoronavirus → SARS-CoV-2   │   │
│  │ Reference: NC_045512.2                                     │   │
│  └────────────────────────────────────────────────────────────┘   │
│                                                                    │
│  ┌──────────┬──────────┬──────────┬──────────┬──────────┬──────┐  │
│  │29,903 bp │ 37.97%   │    12    │  10/12   │    8     │ Rep. │  │
│  │ Genome   │    GC    │   ORFs   │Annotated │  Pfam    │ KEGG │  │
│  └──────────┴──────────┴──────────┴──────────┴──────────┴──────┘  │
│                                                                    │
│  [Summary][QC][ORF Results][Annotation][Taxonomy▲][Dashboard▲]    │
│  [Pathways][Visualizations][AI Report][Downloads]                  │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

**Success States:**
- Analysis complete → show executive panel immediately
- Recent runs listed below upload zone

**Empty State:**
- No prior analysis → show upload zone + workflow cards

---

### Page 3.1 — Upload Page

**Purpose:** File upload with optional configuration.

**Layout:**
```
┌────────────────────────────────────────────────────────────────────┐
│  ← Back    Upload File for Analysis                               │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  Step 1: Select or Drop Your File                                  │
│  ┌──────────────────────────────────────────────────────────────┐ │
│  │   📁 Drag & Drop here                                        │ │
│  │   Supported: .fasta .fastq .fastq.gz .csv .tsv .txt         │ │
│  │   [Browse Files]                                             │ │
│  │   [Load Example: SARS-CoV-2]  [Load Example: DEG Matrix]    │ │
│  └──────────────────────────────────────────────────────────────┘ │
│                                                                    │
│  Step 2: Workflow (auto-detected)                                  │
│  ● Auto (recommended)  ○ Force Workflow A  ○ Force Workflow B     │
│                                                                    │
│  Step 3: Quick Configuration                                       │
│  ▼ Advanced Options (click to expand)                             │
│    Min Sequence Length: [50 ] bp                                  │
│    Max Ambiguity Ratio: [10 ] %                                   │
│    Min ORF Length:      [100] bp                                  │
│    AI Provider: [● Gemini] [○ OpenAI] [○ Claude] [○ None]        │
│    Generate Taxonomy Dashboard: [✓]                               │
│    Offline Mode: [○]                                              │
│                                                                    │
│  [ ▶ Run Analysis ]      [ Reset ]                                │
└────────────────────────────────────────────────────────────────────┘
```

**Validation Rules:**
- File must exist and be non-empty before [Run Analysis] enables
- Unsupported extension → show warning inline
- Size > 500 MB → show size warning
- Workflow A force + CSV gene list → show warning: "This file appears to be a gene list. Auto-detect is recommended."

---

### Page 3.2 — Analysis Progress Page

**Layout:**
```
┌────────────────────────────────────────────────────────────────────┐
│  🔬 Analysis Running — Run: run_20260531_102345                   │
├────────────────────────────────────────────────────────────────────┤
│  File: sars_cov2.fasta │ Workflow A │ AI: Gemini │ Dashboard: ON  │
│                                                                    │
│  ✓ Input Validation      [██████████]  DONE                       │
│  ✓ Quality Control       [██████████]  DONE  — 1 seq, 0 removed  │
│  ✓ ORF Prediction        [██████████]  DONE  — 12 ORFs           │
│  ✓ Protein Translation   [██████████]  DONE  — 12 proteins       │
│  ⟳ Functional Annotation [████░░░░░░]  RUNNING (DIAMOND)...      │
│  ○ Pfam Domain Mapping   [░░░░░░░░░░]  WAITING                   │
│  ○ KEGG Mapping          [░░░░░░░░░░]  WAITING                   │
│  ○ Taxonomy Engine       [░░░░░░░░░░]  WAITING                   │
│  ○ Dashboard Generation  [░░░░░░░░░░]  WAITING                   │
│  ○ AI Interpretation     [░░░░░░░░░░]  WAITING  (Gemini)         │
│  ○ Report Generation     [░░░░░░░░░░]  WAITING                   │
│                                                                    │
│  Live Log ─────────────────────────────────────────────────────   │
│  [10:24:03] ORF_009 → Spike protein (identity: 98.2%)            │
│  [10:24:04] DIAMOND annotation: 10/12 proteins annotated         │
│  [10:24:05] Pfam search started...                               │
│                                                                    │
│  Estimated remaining: ~35 seconds                                 │
│  [✗ Cancel]                                                       │
└────────────────────────────────────────────────────────────────────┘
```

**Workflow B Progress (different steps):**
```
  ✓ Load Expression Matrix    — 2,000 genes
  ✓ Gene ID Normalization     — 1,987 normalized, 13 unmapped
  ✓ FDR Correction (BH)       — applied to 2,000 p-values
  ✓ DEG Classification        — UP:198, DOWN:204, NS:1,598
  ⟳ KEGG ORA Enrichment       — running Fisher's exact test...
  ○ GO Enrichment
  ○ Volcano Plot
  ○ Enrichment Charts
  ○ AI Interpretation (Gemini)
  ○ Report Generation
```

**Error State:**
```
  ✗ Taxonomy Engine   [░░░░░░░░░░]  FAILED
  ⚠ Warning: NCBI Taxonomy API unavailable. 
    Using local cache for SARS-CoV-2 (taxon_id: 2697049).
    Cached data is 12 days old.
  ⟳ Dashboard Generation  (continuing with cached taxonomy)
```

---

### Page 3.3 — Results Dashboard (Full Navigation)

**Tab Bar (11 tabs — Workflow A):**
```
[Summary] [QC Report] [ORF Results] [Annotation] 
[Taxonomy ✦NEW] [Dashboard ✦NEW] [Pathways] 
[Visualizations] [AI Report] [Downloads]
```

**Tab Bar (9 tabs — Workflow B):**
```
[Summary] [Gene Normalization] [DEG Results]
[Enrichment] [GO Analysis] [Visualizations] [AI Report] [Downloads]
```

---

### Page 3.4 — Summary Tab

**Content:**
```
Run Summary
───────────────────────────────────────────────────
Run ID:        run_20260531_102345
Input:         sars_cov2.fasta  (29,903 bp, 1 sequence)
Workflow:      A — Viral Genomics
AI Provider:   Gemini (gemini-1.5-pro)
Status:        Complete
Runtime:       47.2 seconds
───────────────────────────────────────────────────
QC:            1 sequence input, 1 passed
ORFs:          12 predicted (6 forward, 6 reverse frames checked)
Annotated:     10/12 proteins
Pfam Domains:  8 unique domains
KEGG:          4 functional categories
Taxonomy:      SARS-CoV-2 (Betacoronavirus, Coronaviridae)
               Retrieved from: NCBI Taxonomy API
```

---

### Page 3.5 — Taxonomy Tab ← NEW

**Purpose:** Show complete taxonomic classification with lineage table and lineage source.

**Layout:**
```
┌────────────────────────────────────────────────────────────────────┐
│  TAXONOMIC CLASSIFICATION                                          │
├─────────────────────┬──────────────────────────────────────────────┤
│  Realm              │  Riboviria                                   │
│  Kingdom            │  Orthornavirae                               │
│  Phylum             │  Pisuviricota                                │
│  Subphylum          │  —                                           │
│  Class              │  Pisoniviricetes                             │
│  Order              │  Nidovirales                                 │
│  Family             │  Coronaviridae                               │
│  Subfamily          │  Orthocoronavirinae                          │
│  Genus              │  Betacoronavirus                             │
│  Species            │  Severe acute respiratory syndrome-related   │
│                     │  coronavirus                                 │
│  Isolate            │  SARS-CoV-2                                  │
│  NCBI Taxon ID      │  2697049                                     │
│  Host               │  Homo sapiens                               │
│  Genome Type        │  +ssRNA, non-segmented                       │
├─────────────────────┴──────────────────────────────────────────────┤
│  Retrieved from: NCBI Taxonomy API │ Cached: Yes │ Cache age: 0d  │
├────────────────────────────────────────────────────────────────────┤
│  [Download Lineage JSON]  [Download Lineage CSV]                  │
└────────────────────────────────────────────────────────────────────┘
```

**Empty State:**
```
⚠ Taxonomy classification could not be retrieved.
  Possible reasons:
  - No annotation hits found (proteins unannotated)
  - Accession not in NCBI Taxonomy
  - Offline mode with no cached data for this virus
  
  The rest of the analysis is still valid.
  [Retry Taxonomy] [Continue Without Taxonomy]
```

---

### Page 3.6 — Dashboard Tab ← NEW (Major Feature)

**Purpose:** Full-page Virus Genome Intelligence Dashboard. Inspired by NCBI Virus Dashboard. This is the showpiece feature for viva demonstration.

**Layout:**
```
┌────────────────────────────────────────────────────────────────────┐
│  🧬 VIRUS GENOME INTELLIGENCE DASHBOARD                            │
│  SARS-CoV-2 │ Betacoronavirus │ Coronaviridae │ run_20260531      │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  ┌── EXECUTIVE SUMMARY ──────────────────────────────────────────┐ │
│  │ 29,903 bp │ 37.97% GC │ 12 ORFs │ 10 Annotated │ 8 Pfam │    │ │
│  │ Closest: NC_045512.2 (99.8% similar) │ Host: Homo sapiens    │ │
│  └───────────────────────────────────────────────────────────────┘ │
│                                                                    │
│  ┌── TAXONOMY TREE ──────────────────────────────────────────────┐ │
│  │  [Interactive Pyvis tree renders here]                        │ │
│  │  Riboviria → ... → Coronaviridae → Betacoronavirus → SARS-CoV-2│ │
│  └───────────────────────────────────────────────────────────────┘ │
│                                                                    │
│  ┌── GENOME MAP ─────────────────────────────────────────────────┐ │
│  │  [Interactive Plotly genome viewer renders here]               │ │
│  │  ORF1ab──────────────────────────┤S──┤E┤M┤N──┤               │ │
│  │  0                         15000      25000   29903bp         │ │
│  └───────────────────────────────────────────────────────────────┘ │
│                                                                    │
│  ┌── PROTEIN DOMAINS ────────────┐ ┌── VIRUS SIMILARITY ────────┐  │
│  │  [Pfam Treemap (Plotly)]       │ │  [Similarity Bar Chart]     │  │
│  │  RdRp | Spike | Nucleocapsid  │ │  NC_045512.2  ██████ 99.8%  │  │
│  └───────────────────────────────┘ │  SARS-CoV     █████ 79.4%  │  │
│                                    │  RaTG13        ████ 96.2%  │  │
│                                    └────────────────────────────┘  │
│                                                                    │
│  ┌── PATHWAY MAPPING ────────────────────────────────────────────┐ │
│  │  [KEGG Category Sunburst (Plotly)]                            │ │
│  │  Replication | Translation | Host Interaction | Assembly      │ │
│  └───────────────────────────────────────────────────────────────┘ │
│                                                                    │
│  [Open Full Dashboard in New Tab]  [Download Dashboard HTML]      │
│  [Download Dashboard PDF]          [Download All Results ZIP]     │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

**Empty State (Workflow B — Dashboard not applicable):**
```
ℹ Genome Intelligence Dashboard is only available for Workflow A 
  (Viral genomics — FASTA/FASTQ input).
  
  Your current analysis used Workflow B (Functional Genomics).
  See the Enrichment and DEG Results tabs for your results.
```

---

### Page 3.7 — QC Report Tab

Content:
- QC metrics summary table
- Sequences loaded, passed, removed (with reason breakdown)
- Average length, min/max length, GC content, ambiguity distribution
- FASTQ quality histogram (if FASTQ input)

---

### Page 3.8 — ORF Results Tab (Workflow A)

Sortable table:

| ORF ID | Start | Stop | Length (bp) | Length (aa) | Frame | Strand | Status |
|---|---|---|---|---|---|---|---|
| orf_001 | 265 | 21555 | 21291 | 7097 | +1 | Fwd | Annotated |

- Click ORF → expand row to show: DNA sequence (truncated), protein sequence, annotation hit
- Filter: strand, frame, annotation status
- Sort: any column

---

### Page 3.9 — Annotation Tab (Workflow A)

Table: protein annotation results with Pfam domains and KEGG KO numbers.

---

### Page 3.10 — Pathways Tab

**Workflow A:** KEGG Functional Category table (NOT enrichment p-values).
Label at top: "Functional categorization — not enrichment analysis. Small viral proteomes cannot support statistical enrichment."

**Workflow B:** ORA results table:
| Pathway | Gene Count | Universe Count | p-value | FDR | Significant |

---

### Page 3.11 — DEG Results Tab (Workflow B)

Classification table with filter buttons: [All] [Upregulated] [Downregulated] [Not Significant]

---

### Page 3.12 — Gene Normalization Tab (Workflow B)

Table showing per-gene normalization:
| Original ID | Detected Type | HGNC Symbol | Method | Status |

---

### Page 3.13 — Visualizations Tab

**Workflow A charts:** ORF map, domain frequency, domain treemap, KEGG category chart
**Workflow B charts:** Volcano plot, KEGG enrichment chart, GO bubble chart

All charts: download PNG + download SVG + download interactive HTML

---

### Page 3.14 — AI Report Tab

```
┌────────────────────────────────────────────────────────────────────┐
│  🤖 AI-Assisted Biological Interpretation                         │
│  Provider: Google Gemini (gemini-1.5-pro)                         │
├────────────────────────────────────────────────────────────────────┤
│  ⚠ DISCLAIMER: Generated by AI language model. All statistics    │
│  and pathway assignments are computationally derived independently.│
│  AI provides contextual biological explanation only.              │
├────────────────────────────────────────────────────────────────────┤
│  [Interpretation text appears here]                               │
├────────────────────────────────────────────────────────────────────┤
│  [📋 Copy Text]  [📄 Download TXT]                                │
└────────────────────────────────────────────────────────────────────┘
```

**Empty State (AI disabled or key missing):**
```
ℹ AI interpretation is not available.
  Reason: GEMINI_API_KEY environment variable not set.
  
  To enable AI interpretation:
  1. Edit .env file
  2. Add: GEMINI_API_KEY=your_key_here
  3. Or change ai.provider in config.yaml and set the appropriate key
  
  All other analysis results are complete and valid.
```

---

### Page 3.15 — Downloads Tab

```
Downloads — run_20260531_102345

Category: Reports
  📄 Final Report (HTML)         final_report.html            [↓]
  📄 Taxonomy Dashboard (HTML)   taxonomy_dashboard.html      [↓]
  📄 Final Report (PDF)          final_report.pdf             [↓]

Category: Sequence Data
  🧬 Predicted ORFs (FASTA)      predicted_orfs.fasta         [↓]
  🧬 Protein Sequences (FASTA)   predicted_proteins.fasta     [↓]
  📋 GFF3 Annotation             annotation.gff3              [↓]

Category: Analysis Results
  📊 Annotation Results (CSV)    annotation_results.csv       [↓]
  📊 ORF Results (CSV)           orf_results.csv              [↓]
  📊 Pathway Mapping (CSV)       pathway_results.csv          [↓]
  📊 QC Metrics (CSV)            qc_metrics.csv               [↓]

Category: Taxonomy
  📋 Taxonomy Lineage (JSON)     taxonomy_lineage.json        [↓]
  📊 Taxonomy Report (CSV)       taxonomy_report.csv          [↓]

Category: Metadata
  📋 Run Metadata (JSON)         metadata.json                [↓]
  📋 Run Log (TXT)               run.log                      [↓]
  📋 Config Snapshot (YAML)      config_snapshot.yaml         [↓]

  [📦 Download Everything as ZIP]
```

---

## 4. COMPLETE NAVIGATION FLOW DIAGRAM

```
Application Entry
├── First Visit / No Analysis
│   └── Executive Dashboard (Upload Zone)
│       ├── Upload File → Format Detection → Progress → Results
│       ├── [Load Example: SARS-CoV-2] → Auto-load → Progress → Results
│       └── Workflow Cards
│           ├── [Workflow A] → Upload (A) → Progress → Results
│           └── [Workflow B] → Upload (B) → Progress → Results
│
├── After Workflow A Analysis
│   └── Results Dashboard (11 tabs)
│       ├── Summary
│       ├── QC Report
│       ├── ORF Results
│       ├── Annotation
│       ├── Taxonomy ← NEW
│       ├── Dashboard ← NEW (major feature)
│       ├── Pathways
│       ├── Visualizations
│       ├── AI Report
│       └── Downloads (expanded)
│
└── After Workflow B Analysis
    └── Results Dashboard (9 tabs)
        ├── Summary
        ├── Gene Normalization
        ├── DEG Results
        ├── Enrichment
        ├── GO Analysis
        ├── Visualizations
        ├── AI Report
        └── Downloads
```

---

## 5. COMPLETE VALIDATION RULES TABLE

| Input | Rule | Error Message |
|---|---|---|
| File type | Must match supported extensions OR content | "Unsupported format. Supported: .fasta .fastq .csv .tsv .txt" |
| File size | Max 500 MB | "File too large. Max: 500 MB" |
| Empty file | Must not be empty | "Uploaded file is empty." |
| FASTA + ORA request | WorkflowViolationError | "ORA enrichment is not applicable to viral genomic sequences. Workflow A uses functional categorization." |
| Expression matrix | Must have logFC and pvalue columns | "Missing required columns. Need: 'log2FoldChange' and 'pvalue'" |
| Gene list | At least 1 valid ID required | "No valid gene identifiers found in input file." |
| AI provider | Config must name valid provider | "Unknown AI provider. Supported: gemini, openai, claude" |
| Taxonomy API fail | Not fatal | Warning shown in UI; dashboard generates with cached or partial data |

---

*End of Document 3 — App Flow v2.1*
