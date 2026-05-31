# PathoScope AI v2.0
# DOCUMENT 3 — APP FLOW DOCUMENT
# DOCUMENT 4 — UI/UX DESIGN BRIEF

---

# ═══════════════════════════════════════════════════
# DOCUMENT 3 — APP FLOW DOCUMENT
# ═══════════════════════════════════════════════════

---

## 1. SYSTEM INTERACTION OVERVIEW

PathoScope AI v2.0 supports two interaction modes:

**Mode 1 — Command Line Interface (CLI)**
Primary mode. User runs pipeline via terminal. Suitable for advanced users and automated scripting.

**Mode 2 — Web Dashboard (Flask)**
Secondary mode. User uploads files via browser, views results as interactive HTML. Suitable for demo, viva, and non-technical users.

Both modes use the same backend pipeline modules.

---

## 2. CLI USER FLOW

### Flow 2.1 — Standard Single File Analysis
```
Step 1: User runs:
        python main.py --input data/sars_cov2.fasta

Step 2: System detects file format (FASTA)

Step 3: System displays:
        [PathoScope AI v2.0]
        Input: sars_cov2.fasta
        Format: FASTA
        Routing to: Workflow A (Viral Genomics)
        Press ENTER to begin or Ctrl+C to cancel.

Step 4: System runs full pipeline
        Progress bar shows:
        [QC]        ████████░░ 80%
        [ORF]       ██████████ 100%
        [Pfam]      ████████░░ 80%
        [KEGG]      ██████████ 100%
        [AI]        ██████████ 100%
        [Report]    ██████████ 100%

Step 5: System outputs:
        ✓ Analysis complete
        Run ID: run_20260531_102345
        Output directory: ./output/run_20260531_102345/
        HTML Report: ./output/run_20260531_102345/final_report.html
        Runtime: 24.7 seconds
```

### Flow 2.2 — Gene List Analysis
```
Step 1: User runs:
        python main.py --input data/gene_list.txt

Step 2: System detects format (TXT gene list)

Step 3: System displays:
        Format: Gene List (TXT)
        Routing to: Workflow B (Functional Genomics)

Step 4: System runs:
        [GeneNorm]  ██████████ 100%
        [Filter]    ██████████ 100%
        [ORA]       ████████░░ 80%
        [Visualize] ██████████ 100%
        [AI]        ██████████ 100%
        [Report]    ██████████ 100%
```

### Flow 2.3 — Expression Matrix Analysis
```
Step 1: User runs:
        python main.py --input data/expression_matrix.csv

Step 2: System detects: CSV with logFC and pvalue columns

Step 3: Routes to Workflow B

Step 4: Runs:
        [Load Matrix]  ██████████ 100%
        [GeneNorm]     ██████████ 100%
        [FDR Correct]  ██████████ 100%
        [DEG Classify] ██████████ 100%
        [ORA KEGG]     ██████████ 100%
        [ORA GO]       ██████████ 100%
        [Volcano Plot] ██████████ 100%
        [AI Interpret] ██████████ 100%
        [Report]       ██████████ 100%
```

### Flow 2.4 — Override Workflow
```
python main.py --input data/file.fasta --workflow A
python main.py --input data/file.fasta --workflow B
python main.py --input data/file.csv --offline
python main.py --input data/file.fasta --config custom_config.yaml
```

### Flow 2.5 — Error Flow (Invalid File)
```
Step 1: User runs:
        python main.py --input data/invalid.fasta

Step 2: System detects:
        ERROR: Sequence 'seq_001' contains invalid characters: X, Z
        ERROR: 3 sequences exceed ambiguity threshold (N > 10%)

Step 3: System outputs:
        ✗ Validation failed. 
        Log file: ./output/run_20260531_104512/run.log
        Fix errors and rerun, or set --strict-mode=False to skip invalid sequences.
```

---

## 3. WEB DASHBOARD USER FLOW

### Page 3.1 — Dashboard / Home Page

**Purpose:** Entry point. User selects analysis type and uploads file.

**Layout:**
```
┌────────────────────────────────────────────────────────┐
│  🔬 PathoScope AI v2.0                    [About] [Docs]│
├────────────────────────────────────────────────────────┤
│                                                        │
│        Select Your Analysis Type                       │
│                                                        │
│  ┌────────────────────┐  ┌────────────────────────┐   │
│  │  🧬 Workflow A      │  │  📊 Workflow B          │   │
│  │  Viral Genomics     │  │  Functional Genomics    │   │
│  │                     │  │                        │   │
│  │  Input: FASTA/FASTQ │  │  Input: Gene List,     │   │
│  │  Analysis: ORF,     │  │  Expression Matrix     │   │
│  │  Pfam, KEGG         │  │  Analysis: DEG,        │   │
│  │                     │  │  ORA, Pathway          │   │
│  │  [Select →]         │  │  [Select →]            │   │
│  └────────────────────┘  └────────────────────────┘   │
│                                                        │
│  ─── OR ────────────────────────────────────────────  │
│                                                        │
│  [ ↑ Upload Any Supported File — Auto-detect Mode ]   │
│                                                        │
│  Supported: .fasta .fastq .fastq.gz .csv .tsv .txt    │
│                                                        │
│  ─── Recent Runs ───────────────────────────────────  │
│  run_20260531_102345  sars_cov2.fasta  [View Results] │
│  run_20260530_091200  gene_list.txt    [View Results] │
└────────────────────────────────────────────────────────┘
```

**Buttons:**
- [Select →] — opens upload page for selected workflow
- [↑ Upload Any Supported File] — opens generic upload with auto-detect
- [View Results] — opens previous run report

**Validation Rules:**
- File type must match supported extensions
- File size limit: 500MB (configurable)
- Empty files rejected immediately

---

### Page 3.2 — Upload Page

**Purpose:** File upload with configuration options.

**Layout:**
```
┌────────────────────────────────────────────────────────┐
│  ← Back to Home       Workflow A — Viral Genomics      │
├────────────────────────────────────────────────────────┤
│                                                        │
│  Step 1: Upload Your Sequence File                     │
│                                                        │
│  ┌──────────────────────────────────────────────────┐ │
│  │                                                  │ │
│  │   Drag and drop your FASTA or FASTQ file here   │ │
│  │                                                  │ │
│  │            OR   [ Browse Files ]                │ │
│  │                                                  │ │
│  │  Supported: .fasta .fastq .fastq.gz              │ │
│  └──────────────────────────────────────────────────┘ │
│                                                        │
│  Step 2: Configure Parameters (Optional)               │
│                                                        │
│  Min Sequence Length:  [50      ] bp                   │
│  Max Ambiguity Ratio:  [10      ] %                    │
│  Min ORF Length:       [100     ] bp                   │
│  AI Interpretation:    [✓ Enabled]                    │
│  Offline Mode:         [○ Disabled]                   │
│                                                        │
│  Step 3: Start Analysis                                │
│                                                        │
│  [ ▶ Run Analysis ]     [ Reset to Defaults ]         │
│                                                        │
└────────────────────────────────────────────────────────┘
```

**Inputs:** File upload, numeric parameters, toggle switches
**Outputs:** Redirect to Analysis Progress page
**Validation Rules:**
- File must be uploaded before clicking Run
- Min sequence length must be positive integer
- Ambiguity ratio must be 0–100%

---

### Page 3.3 — Analysis Progress Page

**Purpose:** Real-time display of pipeline execution.

**Layout:**
```
┌────────────────────────────────────────────────────────┐
│  🔬 PathoScope AI — Analysis Running                   │
├────────────────────────────────────────────────────────┤
│                                                        │
│  Run ID: run_20260531_102345                           │
│  File: sars_cov2.fasta                                 │
│  Workflow: A (Viral Genomics)                          │
│                                                        │
│  Pipeline Progress:                                    │
│                                                        │
│  ✓ Input Validation         [ ██████████ ] COMPLETE    │
│  ✓ Quality Control          [ ██████████ ] COMPLETE    │
│  ✓ ORF Prediction           [ ██████████ ] COMPLETE    │
│  ✓ Protein Translation      [ ██████████ ] COMPLETE    │
│  ⟳ Pfam Domain Mapping      [ ████░░░░░░ ] RUNNING...  │
│  ○ KEGG Mapping             [ ░░░░░░░░░░ ] WAITING     │
│  ○ Visualization            [ ░░░░░░░░░░ ] WAITING     │
│  ○ AI Interpretation        [ ░░░░░░░░░░ ] WAITING     │
│  ○ Report Generation        [ ░░░░░░░░░░ ] WAITING     │
│                                                        │
│  Live Log:                                             │
│  ┌──────────────────────────────────────────────────┐ │
│  │ [10:23:47] ORF Prediction: 12 ORFs found         │ │
│  │ [10:23:48] Protein Translation: 12 proteins      │ │
│  │ [10:23:50] Pfam search started...                │ │
│  └──────────────────────────────────────────────────┘ │
│                                                        │
│  Estimated time remaining: ~45 seconds                 │
│                                                        │
│  [✗ Cancel Analysis]                                  │
│                                                        │
└────────────────────────────────────────────────────────┘
```

**Success State:** Redirects automatically to Results Dashboard when complete.
**Failure State:** Shows error message, link to log file, option to retry.

---

### Page 3.4 — Results Dashboard

**Purpose:** Display all analysis results in tabbed layout.

**Layout:**
```
┌────────────────────────────────────────────────────────┐
│  ← New Analysis    Run: run_20260531_102345            │
│  sars_cov2.fasta — Workflow A — 24.7 seconds           │
├────────────────────────────────────────────────────────┤
│  [Summary] [QC Report] [ORF Results] [Annotation]      │
│  [Pathways] [Visualizations] [AI Report] [Downloads]   │
├────────────────────────────────────────────────────────┤
│                                                        │
│  ╔═══════════════ SUMMARY TAB ══════════════════════╗  │
│  ║                                                  ║  │
│  ║  Input:     SARS-CoV-2 Genome (29,903 bp)        ║  │
│  ║  Sequences: 1                                    ║  │
│  ║  ORFs Found: 12                                  ║  │
│  ║  Proteins Annotated: 10 / 12                     ║  │
│  ║  Pfam Domains: 8 unique domains                  ║  │
│  ║  KEGG Pathways: 4 pathways mapped                ║  │
│  ║                                                  ║  │
│  ║  ┌─────────────────────────────────────────────┐ ║  │
│  ║  │ [Bar Chart: ORF size distribution]          │ ║  │
│  ║  └─────────────────────────────────────────────┘ ║  │
│  ║                                                  ║  │
│  ╚══════════════════════════════════════════════════╝  │
└────────────────────────────────────────────────────────┘
```

**Tabs:**
- **Summary** — statistics overview, key charts
- **QC Report** — sequence quality metrics table
- **ORF Results** — sortable table of all predicted ORFs
- **Annotation** — protein annotation results table
- **Pathways** — KEGG/GO pathway mapping results
- **Visualizations** — full-size charts (ORF map, domain chart, enrichment chart)
- **AI Report** — plain language biological interpretation
- **Downloads** — download links for all output files

---

### Page 3.5 — QC Report Tab

**Content:**
- Sequences loaded: N
- Sequences after QC: N
- Removed (ambiguity): N
- Removed (duplicates): N
- Removed (length): N
- Average length: X bp
- GC content: X%
- Ambiguity distribution: table + chart
- FASTQ quality score histogram (if FASTQ input)

---

### Page 3.6 — ORF Results Tab

**Content:** Sortable, filterable table:

| ORF ID | Start | Stop | Length | Frame | Strand | Status |
|---|---|---|---|---|---|---|
| ORF_001 | 265 | 21555 | 21291 | +1 | Forward | Annotated |
| ORF_002 | 21563 | 25384 | 3822 | +1 | Forward | Annotated |

- Filter by frame, strand, minimum length
- Click ORF → see protein sequence, annotation details

---

### Page 3.7 — Annotation Tab (Workflow B: DEG Results Tab)

**Workflow A content:** Protein annotation table:
| Protein | Best Match | Identity | E-value | Pfam Domain | KEGG |
|---|---|---|---|---|---|

**Workflow B content:** DEG classification table:
| Gene | logFC | p-value | padj | Classification |
|---|---|---|---|---|
| TP53 | 2.3 | 0.001 | 0.008 | Upregulated |
| BRCA1 | -1.8 | 0.002 | 0.010 | Downregulated |

Filter by: Upregulated | Downregulated | All
Sort by: logFC | padj

---

### Page 3.8 — Visualizations Tab

**Workflow A Charts:**
- ORF Position Map (Matplotlib): Shows all ORFs on genome coordinates
- Domain Frequency Chart (Plotly bar): Top 10 Pfam domains
- KEGG Category Chart (Plotly bar): Functional category distribution

**Workflow B Charts:**
- Volcano Plot (Plotly interactive): logFC vs -log10(padj), colored by classification
- KEGG Enrichment Chart (Plotly bar): Top 20 enriched pathways
- GO Bubble Chart (Plotly): GO terms (x=GeneRatio, y=Term, size=Count, color=padj)

All charts include:
- Download as PNG button
- Download as SVG button
- Hover tooltips with gene/pathway details (Plotly)

---

### Page 3.9 — AI Interpretation Tab

**Layout:**
```
┌────────────────────────────────────────────────────────┐
│  🤖 AI-Assisted Biological Interpretation              │
│     Generated by Claude (claude-sonnet-4-20250514)     │
├────────────────────────────────────────────────────────┤
│  ⚠ DISCLAIMER: This interpretation is generated by    │
│  an AI language model. All statistical values and      │
│  pathway assignments are computationally derived       │
│  and validated independently. AI provides contextual  │
│  biological explanation only.                         │
├────────────────────────────────────────────────────────┤
│                                                        │
│  Interpretation:                                       │
│                                                        │
│  [Plain language biological interpretation             │
│   text from Claude API appears here,                   │
│   explaining the key findings in context               │
│   of the validated computational results]              │
│                                                        │
│  [ 📋 Copy Text ] [ 📄 Download as TXT ]              │
│                                                        │
└────────────────────────────────────────────────────────┘
```

---

### Page 3.10 — Downloads Tab

**Content:**
```
Available Downloads for run_20260531_102345:

📄 Final Report (HTML)          final_report.html          [Download]
📊 Annotation Results (CSV)     annotation_results.csv     [Download]
📊 Pathway Results (CSV)        pathway_results.csv        [Download]
🖼 ORF Map (PNG)               orf_map.png                [Download]
🖼 Enrichment Chart (PNG)      enrichment_chart.png       [Download]
🖼 Volcano Plot (PNG)          volcano_plot.png           [Download]
🗂 Run Log (TXT)               run.log                    [Download]
📋 Metadata (JSON)             metadata.json              [Download]
📦 Download All (ZIP)                                     [Download All]
```

---

## 4. NAVIGATION FLOW DIAGRAM (Text)

```
Home Page
├── Select Workflow A → Upload Page (Workflow A) → Progress → Results Dashboard
├── Select Workflow B → Upload Page (Workflow B) → Progress → Results Dashboard
├── Auto Upload → Format Detection → Correct Upload Page → Progress → Results
└── View Previous Run → Results Dashboard (from cache)

Results Dashboard
├── Summary Tab
├── QC Report Tab
├── ORF / DEG Results Tab
├── Annotation Tab
├── Pathways Tab
├── Visualizations Tab
├── AI Report Tab
└── Downloads Tab
    └── All downloadable outputs
```

---

## 5. CRITICAL VALIDATION RULES (All Pages)

| Input | Validation | Error Message |
|---|---|---|
| File type | Must match supported extensions | "Unsupported file type. Please upload FASTA, FASTQ, CSV, or TXT." |
| File size | Max 500MB | "File too large. Maximum 500MB." |
| Empty file | File must not be empty | "Uploaded file is empty." |
| Min ORF length | Must be integer > 0 | "Min ORF length must be a positive integer." |
| Ambiguity ratio | Must be 0–100 | "Ambiguity ratio must be between 0 and 100." |
| Expression matrix | Must have logFC and pvalue columns | "Expression matrix must contain 'log2FoldChange' and 'pvalue' columns." |
| Gene list | Must contain at least 1 valid gene ID | "No valid gene identifiers found in input file." |

---
---

# ═══════════════════════════════════════════════════
# DOCUMENT 4 — UI/UX DESIGN BRIEF
# ═══════════════════════════════════════════════════

---

## 1. DESIGN PHILOSOPHY

PathoScope AI v2.0 follows the **Scientific Clarity** design philosophy:

> Every design decision must serve the researcher's goal of understanding biological results — not impress with visual complexity.

**Core Principles:**
1. **Data First:** Charts and tables are primary. Decorative elements are secondary.
2. **Information Density over Whitespace:** Scientific software users expect dense, rich information displays. Avoid oversimplified layouts.
3. **No Ambiguity:** Every button, label, and status indicator must be self-explanatory without tooltips.
4. **Reproducibility Visible:** Run IDs, timestamps, and parameter snapshots must be visible at all times.
5. **Error Transparency:** Errors must be shown completely, not hidden behind generic messages.

---

## 2. COLOR PALETTE

```
Primary Background:    #0F1117  (deep dark navy — scientific instrument aesthetic)
Secondary Background:  #1A1D2E  (slightly lighter panel background)
Card Background:       #252842  (elevated card surface)
Border:                #3A3F6E  (subtle border for panels)

Primary Accent:        #4E9FFF  (blue — primary actions, links, workflow A)
Secondary Accent:      #7C63FF  (purple — workflow B, functional genomics)
Success:               #2DD4A0  (teal green — completed steps, upregulated genes)
Warning:               #FFB84E  (amber — warnings, not significant genes)
Danger:                #FF5252  (red — errors, downregulated genes)
Info:                  #63D4FF  (light blue — informational highlights)

Text Primary:          #E8EAF6  (near white — headings, primary labels)
Text Secondary:        #9094B8  (muted blue-gray — secondary labels, metadata)
Text Disabled:         #4A4F72  (dark gray — disabled states)

Chart Colors:
  Upregulated:         #2DD4A0  (green)
  Downregulated:       #FF5252  (red)
  Not Significant:     #4A4F72  (dark gray)
  Pathway Bar:         #4E9FFF  (blue gradient)
  Domain Bar:          #7C63FF  (purple gradient)
```

---

## 3. TYPOGRAPHY

```
Display Font:    JetBrains Mono (monospace, technical, code-like)
  Use for:       Run IDs, gene IDs, sequence identifiers, metrics
  Source:        Google Fonts

Body Font:       IBM Plex Sans
  Use for:       Paragraphs, descriptions, labels, table content
  Source:        Google Fonts

Heading Sizes:
  H1 (Page title):     28px, weight 700, Primary Accent color
  H2 (Section title):  22px, weight 600, Text Primary
  H3 (Card title):     18px, weight 600, Text Primary
  Label:               14px, weight 500, Text Secondary
  Body:                15px, weight 400, Text Primary
  Caption:             12px, weight 400, Text Secondary
  Code/ID:             13px, JetBrains Mono, Info color
```

---

## 4. LAYOUT GUIDELINES

**Grid System:** 12-column CSS Grid
**Max Content Width:** 1400px (centered)
**Sidebar Width:** 240px (when present)
**Gutter:** 16px between columns
**Section Spacing:** 32px between major sections

**Card Component:**
```css
background: var(--card-bg);
border: 1px solid var(--border);
border-radius: 12px;
padding: 24px;
box-shadow: 0 4px 24px rgba(0,0,0,0.3);
```

**Panel System:**
- Left sidebar: navigation between tabs (64px icons + labels)
- Main content: takes remaining width
- Right panel (optional): run metadata, quick actions

---

## 5. COMPONENT SPECIFICATIONS

### 5.1 Navigation Bar
```
Height: 56px
Background: #0F1117 with border-bottom: 1px solid #3A3F6E
Content: [Logo + "PathoScope AI"] [Version badge] [spacer] [Mode indicator] [Help]
Logo: Simple DNA helix SVG icon in Primary Accent blue
```

### 5.2 Workflow Selection Card
```
Width: calc(50% - 16px)
Height: 240px
Hover state: border changes to Primary Accent, slight scale transform (1.02)
Workflow A card: left border accent in #4E9FFF
Workflow B card: left border accent in #7C63FF
Icon: Large (48px) icon in card center
```

### 5.3 Upload Zone
```
Background: #1A1D2E
Border: 2px dashed #3A3F6E
Border-radius: 12px
Height: 180px
Drag active state: border-color becomes Primary Accent, background lightens
```

### 5.4 Progress Bar
```
Background track: #252842
Fill: Primary Accent (#4E9FFF) with shimmer animation when running
Height: 8px, border-radius: 4px
Completed state: fill becomes Success green (#2DD4A0)
Failed state: fill becomes Danger red (#FF5252)
```

### 5.5 Data Tables
```
Header row: background #252842, text Text Secondary, uppercase 12px
Data rows: alternating background (transparent / #1A1D2E)
Hover: background #2D3150
Border: 1px solid #3A3F6E (horizontal only)
Sortable columns: show sort arrow icon in header
Pagination: show 25 rows by default, options: 25/50/100/All
```

### 5.6 Status Badges
```
Upregulated:      background #2DD4A0/15%, text #2DD4A0, border #2DD4A0/30%
Downregulated:    background #FF5252/15%, text #FF5252, border #FF5252/30%
Not Significant:  background #4A4F72/30%, text #9094B8, border #4A4F72/50%
Complete:         background #2DD4A0/15%, text #2DD4A0
Running:          background #4E9FFF/15%, text #4E9FFF (pulse animation)
Failed:           background #FF5252/15%, text #FF5252
```

### 5.7 Buttons
```
Primary:   background #4E9FFF, text #0F1117, border none, hover: #6BAFFF
Secondary: background transparent, text #4E9FFF, border 1px solid #4E9FFF
Danger:    background #FF5252, text white, hover: #FF7070
Disabled:  background #252842, text #4A4F72, cursor not-allowed
Border-radius: 8px
Padding: 10px 20px
Font: IBM Plex Sans, 14px, weight 600
```

### 5.8 Visualization Panels
```
All charts: dark background (#1A1D2E), matching color palette above
Chart title: H3 typography above chart
Chart border: Card component spec
Download buttons: appear on hover over chart panel
Plotly theme: custom dark theme with matching palette
Matplotlib: custom rcParams with matching colors
```

---

## 6. SCREEN DESIGNS BY PAGE

### Screen 6.1 — Home Page (Dashboard)
```
Full dark background
Top navigation bar (fixed)
Hero section: large welcome text + description (centered)
Two workflow cards side-by-side (responsive: stacked on mobile)
Universal upload zone below workflow cards
Recent runs list: compact table at bottom
Footer: version info, team, university
```

### Screen 6.2 — Upload Page
```
Breadcrumb: Home > Workflow A > Upload
Large upload zone (top half)
Configuration accordion below (collapsed by default, expand for options)
Run button: full-width Primary button at bottom
```

### Screen 6.3 — Progress Page
```
Run metadata card at top (Run ID, file, workflow, start time)
Step-by-step progress table (icon + step name + progress bar + status)
Scrollable live log terminal below (dark background, monospace font, auto-scroll)
Cancel button (bottom right, Danger style)
```

### Screen 6.4 — Results Dashboard
```
Run metadata header bar (always visible, collapsible)
Tab navigation (horizontal pill tabs)
Main content area (full width)
Each tab: own layout described in App Flow Document
Sticky download bar at bottom: quick access to key output files
```

---

## 7. RESPONSIVE DESIGN

**Breakpoints:**
- Desktop (≥1280px): Full two-column layout
- Tablet (768–1279px): Single column, sidebar collapses to bottom navigation
- Mobile (< 768px): Stacked layout, simplified navigation, charts replace with summary tables

**Mobile Considerations:**
- Workflow cards stack vertically
- Data tables become horizontally scrollable
- Charts use simplified Plotly mobile layout
- Upload zone uses native file picker (no drag-drop on mobile)

---

## 8. ACCESSIBILITY STANDARDS

- All interactive elements: minimum 44×44px touch target
- Color: all text must meet WCAG AA contrast ratio (4.5:1 minimum)
- Tables: proper ARIA labels, scope attributes
- Charts: alt text descriptions, data available as accessible table
- Error messages: role="alert" for screen readers
- Form inputs: associated labels, required fields marked
- Keyboard navigation: all actions accessible without mouse (Tab, Enter, Arrow keys)

---

## 9. ICON SYSTEM

Use **Heroicons** (MIT license) for all UI icons.
- DNA helix: custom SVG in brand color
- Upload: ArrowUpTrayIcon
- Play/Run: PlayIcon
- Download: ArrowDownTrayIcon
- Check: CheckCircleIcon
- Error: XCircleIcon
- Warning: ExclamationTriangleIcon
- Info: InformationCircleIcon
- Settings: Cog6ToothIcon

---

## 10. LOADING AND ANIMATION

- Page load: fade-in (opacity 0→1, 200ms)
- Progress steps: slide-in from left (150ms stagger between steps)
- Charts: Plotly built-in animation on render
- Table rows: fade-in as data loads (staggered by row)
- Hover states: transition 150ms ease
- Progress bar fill: transition 300ms linear
- No animations that cannot be paused (prefers-reduced-motion respected)

---

*End of Document 4 — UI/UX Design Brief*
