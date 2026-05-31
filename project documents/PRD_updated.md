# PathoScope AI v2.0 — Foundation Package (Revised v2.1)
## National Centre for Bioinformatics (NCB), Quaid-i-Azam University
### Functional Genomics Course — Sir Ghulam Abbas
**Team:** Muhammad Ekramah | Taha Javed | Asad Imam | Mishal Tariq
**Revision:** v2.1 — Taxonomy Dashboard | Provider-Agnostic AI | Enforced Enrichment Rules

---

# ═══════════════════════════════════════════════════════════════
# DOCUMENT 1 — PRODUCT REQUIREMENTS DOCUMENT (PRD) — REVISED v2.1
# ═══════════════════════════════════════════════════════════════

---

## REVISION LOG

| Version | Change | Reason |
|---|---|---|
| v2.0 | Initial dual-workflow design | Original architecture |
| v2.1 | AI provider made configurable (Gemini/OpenAI/Claude) | Avoid vendor lock-in |
| v2.1 | Smith-Waterman removed from requirements | Cost/benefit does not justify |
| v2.1 | Complete test dataset library added | Viva and validation coverage |
| v2.1 | Workflow Router formally documented as safety component | Scientific validity |
| v2.1 | Enrichment analysis rules enforced with biological justification | Scientific correctness |
| v2.1 | Module 10: Virus Taxonomy & Genome Intelligence Dashboard added | Major new feature |
| v2.1 | NCBI Virus Dashboard-inspired interface specified | Reference design standard |

---

## 1. EXECUTIVE SUMMARY

PathoScope AI v2.0 is an automated, dual-workflow bioinformatics platform designed for academic use in functional genomics education and research. The system accepts FASTA sequences, FASTQ reads, gene lists (TXT/CSV), and expression matrices (CSV/TSV), routing them through scientifically validated computational pipelines to produce biological insights, visualizations, and AI-assisted interpretations.

**New in v2.1:** The platform now includes a Virus Taxonomy & Genome Intelligence Dashboard — a major analytical feature inspired by the NCBI Virus Dashboard (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/dashboard). After ORF prediction and annotation are complete for Workflow A inputs, the system automatically classifies the analyzed virus, retrieves its complete taxonomic lineage from NCBI Taxonomy, and generates an interactive intelligence dashboard featuring: taxonomy trees, genome maps, protein domain heatmaps, similarity dashboards, pathway charts, and an executive summary panel.

The AI interpretation layer is now **provider-agnostic**: the system supports Gemini, OpenAI (GPT-4), and Claude through a pluggable adapter architecture. The active provider is selected via config.yaml without any code changes.

PathoScope AI intentionally avoids custom implementations of sequence alignment algorithms and instead relies on validated bioinformatics tools (DIAMOND, BLAST, BioPython alignment modules).

**Platform Type:** Command-line Python pipeline + web-based dashboard
**Primary Language:** Python 3.10+
**Target Environment:** Student laptops, university workstations, Linux/macOS/Windows (WSL)
**Deployment Mode:** Local installation, offline-first, online-enhanced

---

## 2. PRODUCT VISION

> "Give any biology student a raw genomic file or gene list and return a complete, scientifically defensible biological analysis report — including virus taxonomy classification and genome intelligence dashboard — in minutes, with zero manual tool switching."

PathoScope AI aims to be the reference student-facing platform for functional genomics education at QAU, demonstrating that publication-quality analysis can be achieved academically through automated, reproducible pipelines anchored in scientific validity.

---

## 3. PROBLEM STATEMENT

### Problem 3.1 — Workflow Fragmentation
Current functional genomics analysis requires students to manually chain FASTQC, Trimmomatic, BLAST, KEGG REST, Bioconductor, and custom scripts — time-consuming, error-prone, and irreproducible.

### Problem 3.2 — Scientific Inconsistency (Previous Architecture)
The original PathoScope AI applied GSEA and ORA enrichment directly to viral FASTA files. This is biologically invalid:
- GSEA requires pre-ranked gene lists from differential expression experiments.
- ORA requires a meaningful background gene universe; a viral genome with 4–30 proteins cannot provide one.
- Enrichment p-values computed from <50 proteins are statistically meaningless.

### Problem 3.3 — Missing Functional Genomics Support
Gene list processing, expression matrices, Gene ID normalization (HGNC/Ensembl/Entrez), DEG classification, and ORA were absent from the previous design.

### Problem 3.4 — Missing Virus Context and Intelligence
After annotation, there was no way to understand where the analyzed virus sits in the tree of life, what family it belongs to, what hosts it infects, or how it compares to known reference viruses.

### Problem 3.5 — AI Provider Lock-In
The original design hardcoded Claude API. This creates a dependency on a single commercial provider. A configurable adapter design allows switching providers without code changes.

### Problem 3.6 — Insufficient Validation Strategy
No structured test dataset library existed. Testing was ad hoc and insufficient for viva defense.

---

## 4. TARGET USERS

### Primary Users
- Bioinformatics MSc/BSc students at NCB, QAU
- Functional genomics course participants

### Secondary Users
- Research supervisors evaluating project quality
- Virology researchers needing rapid annotation and classification
- Course instructors demonstrating end-to-end genomics workflows

### User Personas

**Persona A — Beginner Student (Amna)**
- Uploads SARS-CoV-2 FASTA
- Expects: clean visual dashboard, virus name, taxonomy tree, plain-English AI interpretation
- Cannot interpret raw BLAST output

**Persona B — Advanced Student (Bilal)**
- Adjusts config.yaml thresholds
- Wants volcano plots, pathway enrichment tables, and reproducible metadata.json
- Will present Taxonomy Dashboard in viva

**Persona C — Professor / Evaluator (Sir Ghulam Abbas)**
- Checks scientific validity of enrichment strategy
- Verifies workflow routing (no ORA on viral genomes)
- Evaluates reproducibility, documentation, and viva demonstration quality

---

## 5. WORKFLOW SAFETY: THE WORKFLOW ROUTER

### 5.1 Scientific Justification

The Workflow Router is not merely a technical routing component — it is a **scientific validity enforcement layer**. Without it, a student could accidentally upload a viral FASTA file and receive statistically meaningless ORA enrichment results, presenting fake p-values in their viva.

The router enforces these biological rules:
- **FASTA input → Workflow A only.** Viral sequences are processed through ORF prediction, annotation, domain mapping, taxonomy classification, and functional categorization. ORA, DEG, and GSEA are never applied.
- **FASTQ input → Workflow A only.** Same as FASTA after quality-based filtering.
- **Gene List (TXT/CSV) → Workflow B only.** Gene IDs are normalized and routed to ORA enrichment against a human gene universe.
- **Expression Matrix (CSV/TSV with logFC + pvalue columns) → Workflow B only.** Full DEG + ORA pipeline.
- **Unknown format → Error.** Pipeline stops with informative message.

### 5.2 Why Enrichment Is Invalid for Small Viral Proteomes

ORA (Fisher's exact test) requires:
1. A meaningful DEG gene list (tens to hundreds of genes)
2. A background gene universe (typically the full genome's protein-coding genes, ~19,000 for Homo sapiens)
3. Pathway gene sets with sufficient overlap to generate statistical power

A typical viral genome provides:
- MS2 Bacteriophage: 4 proteins
- PhiX174: 11 proteins
- SARS-CoV-2: ~12 canonical proteins
- Influenza A (single segment): 2–3 proteins

Testing 4 proteins against a pathway database of thousands of pathways, each with 15–500 genes, produces p-values that are either 0 (pathway genes > DEG list, Fisher's test has no power) or 1 (no genes in common), with no biological meaning. Any "significant" pathway found under these conditions is a statistical artifact, not a biological finding.

**Therefore: Workflow A uses Functional Categorization, not enrichment statistics.**

---

## 6. FUNCTIONAL REQUIREMENTS

### FR-01 — Input Handling
- System SHALL accept FASTA (.fa, .fasta, .fna)
- System SHALL accept FASTQ (.fq, .fastq, .fastq.gz)
- System SHALL accept gene lists (.txt, .csv) with HGNC/Ensembl/Entrez IDs
- System SHALL accept expression matrices (.csv, .tsv) with gene IDs and expression values
- System SHALL auto-detect file format from extension and content
- System SHALL reject malformed files with descriptive error messages

### FR-02 — Workflow Router (Safety Component)
- System SHALL route FASTA/FASTQ exclusively to Workflow A
- System SHALL route gene lists and expression matrices exclusively to Workflow B
- System SHALL log routing decision with reason
- System SHALL allow manual override only via explicit CLI flag
- System SHALL NEVER apply ORA, GSEA, or DEG analysis to Workflow A inputs
- System SHALL NEVER apply ORF prediction or Pfam mapping to Workflow B inputs

### FR-03 — Workflow A: Viral Genomics Pipeline
- System SHALL perform QC: format validation, duplicate removal, ambiguity filtering, length filtering
- System SHALL perform six-frame ORF prediction (BioPython)
- System SHALL translate ORFs to protein sequences (BioPython, standard genetic code)
- System SHALL perform functional annotation via DIAMOND or BLAST (not custom alignment)
- System SHALL map proteins to Pfam domains via pyhmmer
- System SHALL map proteins to KEGG Orthology categories
- System SHALL perform Functional Categorization (NOT enrichment statistics)
- System SHALL route to Taxonomy Engine after annotation
- System SHALL generate Virus Taxonomy & Genome Intelligence Dashboard
- System SHALL generate visualizations: ORF map, domain frequency, KEGG category chart
- System SHALL send validated results to AI Provider for interpretation
- System SHALL generate HTML report with embedded dashboard

### FR-03.5 — Virus Taxonomy & Genome Intelligence Dashboard (NEW)
- System SHALL classify analyzed virus by querying NCBI Taxonomy API
- System SHALL retrieve complete taxonomic lineage: Realm → Kingdom → Phylum → Class → Order → Family → Genus → Species
- System SHALL display an interactive taxonomy tree visualization
- System SHALL generate an Executive Summary card panel showing: Virus Name, Closest Reference, Genome Length, GC Content, ORF Count, Annotated Proteins, Conserved Domains, Top Pathway Category
- System SHALL generate a Genome Map showing ORF positions and protein locations
- System SHALL generate a Protein Domain Dashboard with treemap and bar charts
- System SHALL generate a Virus Similarity Dashboard comparing to closest NCBI reference genomes
- System SHALL display Host Distribution information when metadata is available
- System SHALL generate a Pathway Mapping Dashboard showing KEGG category distribution
- System SHALL work in offline mode using locally cached taxonomy data
- System SHALL export dashboard as standalone HTML + PDF

### FR-04 — Workflow B: Functional Genomics Pipeline
- System SHALL normalize gene identifiers to HGNC symbols (HGNC/Ensembl/Entrez supported)
- System SHALL apply statistical filtering: logFC threshold, p-value cutoff, BH FDR correction
- System SHALL classify genes: Upregulated (logFC > 1, padj < 0.05), Downregulated (logFC < -1, padj < 0.05), Not Significant
- System SHALL perform ORA enrichment using Fisher's exact test against human gene universe
- System SHALL perform GO Biological Process enrichment
- System SHALL generate volcano plots, enrichment bar charts, GO bubble plots
- System SHALL send validated results to AI Provider for interpretation
- System SHALL generate HTML report

### FR-05 — Statistical Analysis
- System SHALL apply Benjamini-Hochberg FDR to all multiple hypothesis tests
- System SHALL use Fisher's exact test for ORA (Workflow B only)
- System SHALL use configurable thresholds from config.yaml
- System SHALL NEVER compute or report enrichment statistics for Workflow A inputs
- System SHALL NEVER fabricate statistical values

### FR-06 — Provider-Agnostic AI Interpretation
- System SHALL support multiple AI providers: Gemini, OpenAI, Claude
- Active provider SHALL be selected via config.yaml (ai.provider field)
- System SHALL implement AIProvider adapter pattern (GeminiProvider, OpenAIProvider, ClaudeProvider)
- AI SHALL only interpret validated computational outputs
- AI SHALL NOT generate statistics, p-values, pathway assignments, ORFs, or taxonomic classifications
- AI outputs SHALL include a mandatory disclaimer
- System SHALL continue without AI if provider is unavailable (log warning, skip section)

### FR-07 — Reporting and Export
- System SHALL generate HTML report for every run
- System SHALL generate CSV tables for all results
- System SHALL generate PDF report (WeasyPrint)
- System SHALL generate GFF3 file for ORF annotation results (Workflow A)
- System SHALL generate JSON metadata for every run (metadata.json)
- System SHALL save all outputs to timestamped run directory

### FR-08 — Offline/Online Mode
- System SHALL work fully offline using locally cached databases and taxonomy data
- System SHALL cache all API responses in SQLite
- System SHALL log data source (online/offline) for every lookup

### FR-09 — Reproducibility
- System SHALL read all parameters from config.yaml
- System SHALL write metadata.json with: input MD5, all parameters, tool versions, database versions
- System SHALL support rerunning any prior analysis from its run ID

---

## 7. NON-FUNCTIONAL REQUIREMENTS

| ID | Requirement | Target |
|---|---|---|
| NFR-01 | Runtime on student laptop (4-core) | < 15 min for SARS-CoV-2 end-to-end |
| NFR-02 | Cross-platform | Python 3.10+ on Linux, macOS, Windows WSL |
| NFR-03 | Dependency specification | requirements.txt + environment.yml |
| NFR-04 | Offline operation | All core features work without internet |
| NFR-05 | Database version tracking | Recorded in metadata.json |
| NFR-06 | Modularity | One function per biological task, max 50 lines |
| NFR-07 | Test coverage | > 90% unit test pass rate |
| NFR-08 | Chart quality | All PNG outputs at 300 DPI minimum |
| NFR-09 | AI agnosticism | Provider swappable via config only |
| NFR-10 | No custom alignment | Use DIAMOND/BLAST/BioPython only |

---

## 8. SCIENTIFIC POLICY STATEMENT

> PathoScope AI v2.0 intentionally avoids custom implementations of sequence alignment algorithms and instead relies on validated bioinformatics tools (DIAMOND, BLAST, BioPython alignment modules). Custom alignment implementation (e.g. Smith-Waterman from scratch) has been evaluated and excluded because its engineering cost exceeds its educational value in the context of this functional genomics project.

> AI components (Gemini, OpenAI, Claude) are used exclusively for natural language interpretation of validated computational results. No AI component performs ORF prediction, sequence annotation, statistical analysis, pathway enrichment, or taxonomic classification. These remain deterministic, auditable, scientifically defensible computational modules.

---

## 9. COMPLETE SCIENTIFIC VALIDATION DATASET LIBRARY

### 9.1 Viral FASTA Datasets

**Dataset FASTA-01: MS2 Bacteriophage**
- NCBI Accession: NC_001417
- Genome size: 3,569 bp
- Expected ORFs: 4 canonical (maturation, coat, lysis, replicase)
- Purpose: Minimal genome baseline; fast test case; well-characterized
- Expected QC: 1 sequence, passes all filters
- Expected Annotation: coat protein (Pfam PF01819), replicase
- Expected Taxonomy: Leviviricetes → Norzivirales → Fiersviridae → Levivirus → MS2
- Validation Criteria: ≥4 ORFs, ≥1 Pfam domain, correct taxonomy lineage retrieved

**Dataset FASTA-02: PhiX174**
- NCBI Accession: NC_001422
- Genome size: 5,386 bp
- Expected ORFs: 11 canonical
- Purpose: Classic benchmark genome; widely used as sequencing control
- Expected QC: 1 sequence, passes all filters
- Expected Taxonomy: Phiacceleraviricetes → Petitvirales → Microviridae → Bullavirinae
- Validation Criteria: 10–14 ORFs, taxonomy retrieved, capsid protein domain detected

**Dataset FASTA-03: SARS-CoV-2**
- NCBI Accession: NC_045512.2
- Genome size: 29,903 bp
- Expected ORFs: ~12 canonical (ORF1ab, S, E, M, N + accessory ORFs)
- Purpose: Real viral genome; publicly known annotation for validation
- Expected Taxonomy: Riboviria → Orthornavirae → Pisuviricota → Pisoniviricetes → Nidovirales → Coronaviridae → Betacoronavirus → SARS-CoV-2
- Validation Criteria: Spike protein detected (>3000 bp ORF), Pfam RdRp domain detected, full taxonomy lineage retrieved, dashboard renders correctly

**Dataset FASTA-04: Influenza A (H1N1)**
- NCBI Accession: CY121680 (segment 4, HA gene) or full genome set
- Genome size: segmented (~13,500 bp total, 8 segments)
- Expected ORFs: 2–3 per segment
- Purpose: Test segmented genome handling
- Expected Taxonomy: Insthoviricetes → Articulavirales → Orthomyxoviridae → Alphainfluenzavirus
- Validation Criteria: HA protein detected, neuraminidase domain found, taxonomy correct

### 9.2 Viral FASTQ Datasets

**Dataset FASTQ-01: Valid SARS-CoV-2 Short Reads**
- Source: SRA sample ERR5556343 (amplicon reads, synthetic)
- Purpose: Validate FASTQ QC pipeline, quality filtering, Phred score handling
- Expected: Reads pass QC after quality filtering at Phred ≥ 20
- Validation Criteria: QC summary CSV generated, low-quality reads flagged

**Dataset FASTQ-02: PhiX Control Reads**
- Source: Illumina PhiX control reads
- Purpose: Baseline sequencing quality control test
- Expected: High-quality reads, minimal filtering required
- Validation Criteria: >95% reads pass Phred 20 threshold

**Dataset FASTQ-03: Corrupted FASTQ**
- Manually constructed: truncated quality lines, mismatched lengths
- Purpose: Test error handling
- Expected: FileFormatError raised with clear message identifying corrupt records
- Validation Criteria: System does not crash; error message names the problematic record and line number

### 9.3 Functional Genomics Datasets

**Dataset FG-01: HGNC Gene Symbol List**
- Content: 100 well-characterized human gene symbols (TP53, BRCA1, BRCA2, MYC, EGFR, PTEN, RB1, APC, VHL, MLH1, and 90 more)
- Purpose: Test gene ID normalization for HGNC input (identity conversion)
- Expected: 100/100 recognized as HGNC, method='identity'
- Validation Criteria: No API calls made for recognized HGNC symbols

**Dataset FG-02: Ensembl ID List**
- Content: 100 Ensembl IDs (ENSG format) mapping to same genes as FG-01
- Purpose: Test Ensembl→HGNC conversion via local table
- Expected: 100/100 converted via local table (no API required)
- Validation Criteria: method='local_table' for all; confirm TP53 = ENSG00000141510

**Dataset FG-03: Entrez ID List**
- Content: 100 Entrez IDs mapping to same genes
- Purpose: Test Entrez→HGNC conversion
- Expected: 100/100 converted via local table
- Validation Criteria: method='local_table'; confirm TP53 Entrez=7157

**Dataset FG-04: Mixed ID List**
- Content: 50 HGNC + 25 Ensembl + 25 Entrez IDs (shuffled)
- Purpose: Test multi-format auto-detection in a single file
- Expected: All 100 normalized to HGNC, auto-type detection correct per row
- Validation Criteria: Per-row type detection logged in normalization CSV

**Dataset FG-05: Differential Expression Matrix**
- Content: 2,000-gene CSV with columns: gene_id, log2FoldChange, pvalue
  - 200 upregulated (logFC > 1, pvalue < 0.05)
  - 200 downregulated (logFC < -1, pvalue < 0.05)
  - 1,600 not significant
- Purpose: Full DEG classification + ORA validation
- Expected: TP53 pathway genes enriched; Cell cycle pathway FDR < 0.05
- Validation Criteria: DEG counts match expected; volcano plot visually correct

**Dataset FG-06: DEG Matrix with Missing Values**
- Content: Expression matrix with 50 rows containing NaN in pvalue column
- Purpose: Test missing value handling
- Expected: NaN rows flagged and removed with count logged
- Validation Criteria: System does not crash; validation CSV shows removed rows

### 9.4 Error/Edge Case Datasets

**Dataset ERR-01: Invalid Characters in FASTA**
- Content: FASTA with sequences containing X, Z, J characters
- Expected: Sequences rejected; error message lists each invalid character found; run continues with valid sequences unless --strict flag set

**Dataset ERR-02: High Ambiguity FASTA**
- Content: FASTA with sequences where >10% bases are N
- Expected: Sequences filtered; count logged in QC report

**Dataset ERR-03: Duplicate Sequences**
- Content: FASTA with 5 identical copies of the same sequence
- Expected: 1 sequence retained, 4 removed; QC report shows removed_duplicates=4

**Dataset ERR-04: Empty File**
- Expected: InputValidationError immediately; no pipeline steps run

**Dataset ERR-05: Wrong Workflow File Type**
- Content: FASTA file with misleading .csv extension
- Expected: Content-based detection overrides extension; FASTA format detected; routed to Workflow A

---

## 10. MVP SCOPE (Version 1.0 / Submission Build)

**In Scope:**
- Workflow A: FASTA/FASTQ → QC → ORF → Translation → Pfam → KEGG → Functional Categorization → Taxonomy Dashboard → AI Interpretation → Report
- Workflow B: Gene List/Expression Matrix → Gene ID Normalization → DEG Classification → KEGG ORA → GO ORA → Volcano Plot → Enrichment Charts → AI Interpretation → Report
- Taxonomy Dashboard (Module 10): Taxonomy Tree, Executive Panel, Genome Map, Protein Domain Dashboard, Virus Similarity
- Provider-agnostic AI (Gemini, OpenAI, Claude adapters)
- Offline mode with SQLite caching
- config.yaml system + metadata.json reproducibility

**Explicitly Out of Scope for v1.0:**
- GSEA (biologically inappropriate for available input types)
- Custom Smith-Waterman implementation (excluded by scientific policy)
- D3.js custom visualizations (complexity not justified; Plotly sufficient)
- Multi-omics integration
- Cloud deployment
- Real-time collaborative analysis
- Geographic distribution map (optional advanced feature, deferred)
- Multi-sample comparative genomics across uploaded files (deferred to v2)

---

## 11. SUCCESS METRICS

| Metric | Target |
|---|---|
| SARS-CoV-2 FASTA end-to-end | < 10 minutes including dashboard |
| Taxonomy retrieved for SARS-CoV-2 | Correct lineage to species level |
| Gene list normalized (500 genes) | < 2 minutes |
| Offline mode functional | All core features work without internet |
| Report + Dashboard generated per run | Yes — HTML always generated |
| Tests passing on all 14 test datasets | > 90% pass rate |
| Reproducibility | Same input → identical outputs |
| Viva live demo | Completes on student laptop |
| AI provider switch (Gemini ↔ OpenAI ↔ Claude) | Config change only, no code change |

---

## 12. USER STORIES

**US-01:** As a student, I upload a SARS-CoV-2 FASTA file and receive a complete analysis report showing the virus taxonomy, genome map, annotated proteins, domain distribution, and plain-language AI interpretation.

**US-02:** As a researcher, I upload a differential expression matrix and receive a volcano plot, enriched pathways with FDR values, classified DEG lists, and AI-generated biological commentary.

**US-03:** As a supervisor, I review the metadata.json and confirm the exact thresholds used, tool versions, and database versions — ensuring the analysis is reproducible.

**US-04:** As a student demonstrating in viva, I show the Taxonomy Dashboard for SARS-CoV-2 displaying its lineage tree, protein domains, and similarity to related coronaviruses without internet access (offline mode).

**US-05:** As an advanced user, I change ai.provider from 'gemini' to 'claude' in config.yaml and rerun without modifying any code — the AI section of the report now uses Claude instead.

**US-06:** As a student, I upload a gene list with mixed Ensembl, Entrez, and HGNC identifiers and receive a normalized output where all genes are represented by HGNC symbols with conversion method documented per gene.

---

## 13. RISKS AND LIMITATIONS

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| NCBI Taxonomy API unavailable | Medium | Taxonomy Dashboard degraded | Local taxonomy cache in SQLite; fallback to cached lineage for known viruses |
| KEGG API rate limiting | High | Annotation delayed | Rate limiter + SQLite cache with 7-day TTL |
| AI provider API key not set | High | AI section empty | Graceful degradation; report generated without AI section |
| Ensembl API downtime | Medium | Gene normalization delayed | Local mapping table covers >95% of common human genes |
| Small laptop performance | Medium | Timeout on large FASTQ | Progress bar with estimated time; FASTQ size warning at upload |
| pyhmmer Pfam database size | Low | Disk space on student PC | Provide viral Pfam subset (~50 MB) instead of full Pfam-A (6 GB) |

---

## 14. FEATURES EXCLUDED FROM v1.0 WITH JUSTIFICATION

| Feature | Reason Excluded |
|---|---|
| GSEA | Requires pre-ranked DEG lists from RNA-seq; not applicable to FASTA or small gene lists without expression data |
| Smith-Waterman from scratch | Engineering cost >> educational value; use BLAST/DIAMOND |
| D3.js custom charts | JavaScript complexity not justified; Plotly achieves same interactive results in Python |
| Cloud deployment | Out of scope for academic project; local deployment sufficient for viva |
| Geographic distribution map | Requires sample metadata not available in input files; deferred to v2 |
| Multi-sample comparative analysis | Requires batch upload UI; deferred to v2 |

---

*End of Document 1 — PRD v2.1*
