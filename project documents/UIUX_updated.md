# PathoScope AI v2.0 — DOCUMENT 4: UI/UX DESIGN BRIEF — REVISED v2.1

---

# ═══════════════════════════════════════════════════════════════
# DOCUMENT 4 — UI/UX DESIGN BRIEF — REVISED v2.1
# ═══════════════════════════════════════════════════════════════

**Changes in v2.1:**
- NCBI Virus Dashboard design language added as reference standard
- Dark/Light mode system added
- Taxonomy Dashboard visual component specifications added
- Genome Map component specified
- Taxonomy Tree component specified
- Protein Domain Treemap specified
- Virus Similarity Dashboard specified
- Dashboard layout grid specified
- Scientific color palette extended for taxonomy visualization
- Metric card component system added

---

## 1. DESIGN PHILOSOPHY

PathoScope AI v2.1 follows two complementary design principles:

**Principle 1 — Scientific Clarity**
Every visual element exists to communicate biological results. Nothing is decorative without purpose. The interface should look like a professional bioinformatics tool, not a startup product.

**Principle 2 — NCBI Intelligence Dashboard Aesthetic**
The Virus Taxonomy & Genome Intelligence Dashboard takes visual inspiration from the NCBI Virus Dashboard (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/dashboard). Key characteristics of that design that we adopt:
- Clean card-based information panels
- Metric cards for immediate quantitative overview
- Interactive zoomable charts
- Scientific color palette (not consumer/marketing colors)
- Taxonomy tree prominently featured
- Genome map as central visual element
- Publication-quality data density

**What to Avoid:**
- Purple gradient backgrounds (clichéd AI product aesthetic)
- Floating animations unrelated to data
- Rounded buttons that look like mobile apps
- Overly large whitespace that wastes screen space on scientific data
- Flashy loading animations for biological analysis steps

---

## 2. DUAL THEME SYSTEM (Dark + Light)

PathoScope AI v2.1 provides both dark and light modes. Default is **dark mode** (preferred for bioinformatics software, matches NCBI Virus Dashboard's dark panels). Users toggle via the nav bar.

### Dark Theme (Default)
```
--bg-primary:        #0D1117    (deep dark — GitHub/NCBI inspired)
--bg-secondary:      #161B22    (panel background)
--bg-card:           #1C2128    (elevated card)
--bg-code:           #0A0E13    (terminal / log panel)
--border:            #30363D    (subtle border)
--border-active:     #484F58    (active/hover border)

--text-primary:      #E6EDF3    (high contrast near-white)
--text-secondary:    #8B949E    (muted description text)
--text-disabled:     #484F58    (disabled elements)
--text-code:         #79C0FF    (monospace code text, gene IDs)

--accent-blue:       #388BFD    (primary actions, links, Workflow A)
--accent-purple:     #BC8CFF    (Workflow B, functional genomics)
--accent-green:      #3FB950    (success, upregulated, complete)
--accent-red:        #F85149    (error, downregulated, danger)
--accent-orange:     #D29922    (warning, borderline results)
--accent-teal:       #56D364    (taxonomy highlighting)
--accent-gold:       #E3B341    (active taxonomy tree path)

--chart-up:          #3FB950    (upregulated genes)
--chart-down:        #F85149    (downregulated genes)
--chart-neutral:     #30363D    (not significant)
--chart-1:           #388BFD    (primary chart series)
--chart-2:           #BC8CFF    (secondary chart series)
--chart-3:           #56D364    (tertiary chart series)
--chart-4:           #D29922    (quaternary chart series)
```

### Light Theme
```
--bg-primary:        #FFFFFF
--bg-secondary:      #F6F8FA
--bg-card:           #FFFFFF
--bg-code:           #F6F8FA
--border:            #D0D7DE
--text-primary:      #1F2328
--text-secondary:    #656D76
--accent-blue:       #0969DA
--accent-green:      #1A7F37
--accent-red:        #CF222E
(all other accents adjusted proportionally for light background)
```

---

## 3. TYPOGRAPHY

```
Monospace Font: JetBrains Mono
  Use for:  Run IDs, gene IDs, sequence IDs, accession numbers,
            DNA/protein sequences, metric values, code in report
  Sizes:    Run ID labels: 11px / Gene IDs: 13px
  Source:   Google Fonts

Display Font: Inter Variable
  Use for:  Page titles, section headings, navigation labels
  Sizes:    Page Title (H1): 24px, weight 600
            Section Heading (H2): 18px, weight 600
            Card Heading (H3): 15px, weight 600
  Source:   Google Fonts (Inter Variable — all weights in one file)

Body Font: Inter Variable (regular weight)
  Use for:  Body text, descriptions, table content, labels
  Sizes:    Body: 14px, weight 400
            Caption: 12px, weight 400
            Table: 13px, weight 400

Scientific Data:
  All numerical values (p-values, fold changes, lengths, 
  accession numbers): JetBrains Mono, 13px
```

---

## 4. LAYOUT SYSTEM

**Grid:** 12-column CSS Grid
**Container max-width:** 1440px (left-centered for wide monitors)
**Content max-width:** 1200px
**Sidebar width:** 220px (collapsed: 56px icon-only)
**Main content:** fluid, fills remaining width
**Card gutter:** 12px
**Section spacing:** 24px

**Layout Zones:**
```
┌──────────────────────────────────────────────────────────────────┐
│  TOP NAV BAR (56px fixed)                                        │
├──────────┬───────────────────────────────────────────────────────┤
│  SIDE    │  MAIN CONTENT AREA                                    │
│  NAV     │                                                        │
│  (220px) │  [Tab Bar / Page Header]                              │
│          │                                                        │
│  ●       │  [Main Content Grid]                                  │
│  ●       │                                                        │
│  ●       │                                                        │
│          │                                                        │
│  ●       │  [Footer / Run Info Bar]                              │
└──────────┴───────────────────────────────────────────────────────┘
```

---

## 5. COMPONENT SPECIFICATIONS

### 5.1 Navigation Bar (Top)
```
Height:      56px
Background:  var(--bg-secondary)
Border:      0 0 1px 0 solid var(--border)
Left:        [🔬 PathoScope AI logo] [version badge "v2.1"]
Center:      [Run ID badge if analysis active]
Right:       [AI provider indicator] [🌙 dark/light toggle] [⚙ settings] [? help]

AI Provider Indicator:
  Background: var(--bg-card)
  Border: 1px solid var(--border)
  Shows: "AI: Gemini" | "AI: OpenAI" | "AI: Claude" | "AI: Off"
  Color: accent-teal for active provider
```

### 5.2 Side Navigation
```
Width:      220px (expanded), 56px (collapsed)
Background: var(--bg-secondary)
Border:     0 1px 0 0 solid var(--border)

Items (icons + labels):
  🏠 Dashboard
  ↑ New Analysis
  📁 Recent Runs
  ─────────── (divider)
  [Active run tabs appear here when run is selected]
  ─────────── (divider)
  📚 Documentation
  ⚙ Settings
```

### 5.3 Metric Cards (Executive Panel)
Six cards in a responsive 6-column grid (2-column on tablet, 3x2 on mobile).

```
Card dimensions: min-width 140px, height 80px
Background: var(--bg-card)
Border: 1px solid var(--border)
Border-radius: 8px
Hover: border-color: var(--border-active), slight translateY(-2px)

Layout inside card:
  Top: label text (12px, --text-secondary, uppercase, JetBrains Mono)
  Middle: value (22px, --text-primary, JetBrains Mono, weight 600)
  Bottom: sublabel (11px, --text-secondary)

Example card:
┌──────────────────┐
│  GENOME LENGTH   │  ← label (12px, muted)
│  29,903 bp       │  ← value (22px, primary, monospace)
│  Single segment  │  ← sublabel (11px, muted)
└──────────────────┘
```

### 5.4 Taxonomy Tree Component
```
Container: full width, min-height 300px
Background: var(--bg-secondary)
Border-radius: 8px
Border: 1px solid var(--border)

Technology: PyVis generates self-contained HTML
  iframe embedded in dashboard OR injected as inline HTML

Node styles:
  Normal node:     border-radius 50%, background var(--bg-card),
                   border var(--border), text var(--text-secondary)
  Active path:     border var(--accent-gold), 
                   background rgba(227,179,65,0.15) — gold highlight
  Query node (leaf): larger size, background var(--accent-blue)/20%,
                     border var(--accent-blue), label bold

NetworkX graph construction:
  nodes = taxonomy levels (Realm, Kingdom, ..., Species, Query)
  edges = parent → child (directed)
  layout = hierarchical top-down (sugiyama / dot layout)
  
PyVis options:
  physics: disabled (hierarchical layout is static)
  node font: JetBrains Mono 12px
  edge arrows: show direction
  zoom: enabled
  pan: enabled
  tooltip on hover: shows NCBI taxon ID + description
```

### 5.5 Genome Map Component (NCBI-style)
```
Technology: Plotly horizontal bar chart with custom styling

Layout:
  X-axis: genomic position (0 to genome_length)
  Y-axis: protein name labels (ORF_ID → annotation_name if available)
  Bars: horizontal, colored by protein function category

Color coding (protein categories):
  Replication proteins:  var(--accent-blue)
  Structural proteins:   var(--accent-green)
  Accessory proteins:    var(--accent-orange)
  Unknown function:      var(--text-disabled)

Bar styling:
  Height: 0.6 (Plotly bar width)
  Border-radius: 3px (rounded bars)
  No gaps between bars on same track

Annotations on bars:
  Label each bar with ORF ID (if bar is wide enough > 500 bp)
  Font: JetBrains Mono 10px, color white

Genome ruler:
  X-axis ticks every 5,000 bp
  Major ticks at multiples of 10,000 bp
  Label format: "10 kb", "20 kb", "29.9 kb"

Hover tooltip:
  ORF ID | Protein Name | Start-Stop | Length | Identity% | E-value

Background: var(--bg-secondary)
Paper background: var(--bg-secondary)
Font color: var(--text-primary)
```

### 5.6 Protein Domain Treemap
```
Technology: Plotly Treemap

Data structure:
  Level 1: Pfam family category (e.g., Viral core, RNA synthesis)
  Level 2: Individual Pfam domain (e.g., PF00680 RdRp_1)
  Size: number of proteins with this domain
  Color: E-value quality (low e-value = dark accent-blue, high = muted)

Styling:
  Background: var(--bg-secondary)
  Text: JetBrains Mono for accessions, Inter for domain names
  Hover: domain description, accession, count, best e-value
```

### 5.7 Virus Similarity Dashboard
```
Technology: Plotly horizontal bar chart

Data: top 10 most similar reference viruses from virus_reference_db

Layout:
  Y-axis: virus name + accession (JetBrains Mono)
  X-axis (top): Similarity % (0-100)
  X-axis (bottom): Coverage % (0-100)
  Bars: split (similarity = solid, coverage = outlined)
  
  First bar always: query virus vs itself (100%) — shown as reference

Color: gradient from --accent-green (>90% similar) to --accent-orange (50-70%) to --text-disabled (<50%)
Hover: full virus name, accession, genome size, host, similarity score
```

### 5.8 Volcano Plot Component
```
Technology: Plotly scatter

Points:
  Upregulated (logFC>1, padj<0.05):   --chart-up (#3FB950), opacity 0.8
  Downregulated (logFC<-1, padj<0.05): --chart-down (#F85149), opacity 0.8
  Not significant:                      --chart-neutral, opacity 0.4, size 4

Threshold lines:
  Vertical x=1 and x=-1: dashed, --text-secondary, opacity 0.6
  Horizontal -log10(0.05)=~1.3: dashed, --text-secondary, opacity 0.6

Top gene labels:
  Top 10 most significant genes labeled with gene symbol
  Label font: JetBrains Mono 10px
  Arrows pointing to labeled points

Axis titles:
  X: "log₂ Fold Change" (subscript 2)
  Y: "-log₁₀(adjusted p-value)"

Legend: bottom-right corner, shows counts per category
```

### 5.9 KEGG Enrichment Bar Chart
```
Technology: Plotly horizontal bar

Data: Top 20 significant pathways (sorted by padj ascending = most significant on top)
X-axis: -log10(padj)
Y-axis: Pathway name (truncated at 50 chars, full name in tooltip)
Bar color: encodes gene_count using sequential colorscale (light to dark --accent-blue)
Color bar legend: "Gene Count"

Hover: pathway name, gene_count, universe_count, pvalue, padj, gene_list (truncated)

Dashed line at x = -log10(0.05) = 1.3: "FDR threshold"
```

### 5.10 Pathway Mapping Sunburst (Workflow A — Functional Categorization)
```
Technology: Plotly Sunburst

Inner ring: KEGG top-level categories (Metabolism, Genetic Info Processing, etc.)
Outer ring: Sub-categories
Segment size: protein count in each category

Hover: category name, protein count, example protein
Center label: "KEGG Functional Categories"

Note panel below chart: 
"Functional categorization (not enrichment statistics) — 
 appropriate for small viral proteomes."
```

### 5.11 GO Bubble Chart
```
Technology: Plotly scatter

X-axis: GeneRatio (genes_in_term / total_degs)
Y-axis: GO term name (top 20, sorted by padj)
Bubble size: gene count
Bubble color: padj colorscale (red=significant, blue=less significant)

Legend: colorbar for padj + size legend for gene count
```

### 5.12 Data Tables
```
Header: background var(--bg-secondary), text var(--text-secondary),
        uppercase 11px, JetBrains Mono
Rows:   alternating transparent / var(--bg-secondary) at 40% opacity
Hover:  background var(--bg-card)
Border: 1px solid var(--border) (horizontal dividers only)

Sortable columns: chevron icon appears on hover over sortable headers
Active sort: accent color chevron

Pagination: bottom of table
  Show: [25] [50] [100] [All] rows per page

Compact mode for large tables: row height 32px (default 40px)
Full mode for annotation tables: row height 48px (for multi-line content)
```

### 5.13 Status Badges
```
Complete:       bg rgba(63,185,80,0.15), text #3FB950, border rgba(63,185,80,0.3)
Running:        bg rgba(56,139,253,0.15), text #388BFD, animate pulse
Failed:         bg rgba(248,81,73,0.15),  text #F85149, border rgba(248,81,73,0.3)
Warning:        bg rgba(210,153,34,0.15), text #D29922, border rgba(210,153,34,0.3)
Upregulated:    bg rgba(63,185,80,0.12),  text #3FB950
Downregulated:  bg rgba(248,81,73,0.12),  text #F85149
Not Significant:bg rgba(48,54,61,0.5),    text #8B949E
```

### 5.14 Progress Steps
```
Bar track:  height 4px, background var(--bg-card), border-radius 2px
Fill:       var(--accent-blue) with shimmer animation (running)
             var(--accent-green) (complete)
             var(--accent-red) (failed)
             var(--accent-orange) (warning — completed with issues)

Step icon:
  Waiting:   ○ (open circle, --text-disabled)
  Running:   ⟳ (spinning, --accent-blue)
  Complete:  ✓ (filled, --accent-green)
  Failed:    ✗ (filled, --accent-red)
  Warning:   ⚠ (filled, --accent-orange)
```

### 5.15 Upload Zone
```
Background:    var(--bg-secondary)
Border:        2px dashed var(--border)
Border-radius: 12px
Height:        160px

Drag active state:
  border-color: var(--accent-blue)
  background:   rgba(56,139,253,0.05)
  transition:   150ms ease

File accepted state:
  border-color: var(--accent-green)
  Shows: filename, size, detected format badge
  
Error state:
  border-color: var(--accent-red)
  Shows: error message with specific reason
```

### 5.16 Buttons
```
Primary:    bg var(--accent-blue), text #0D1117 (dark), weight 600
            hover: brightness(1.1), transition 150ms
Secondary:  bg transparent, text var(--accent-blue), 
            border 1px solid var(--accent-blue)
Danger:     bg var(--accent-red), text white
Disabled:   bg var(--bg-card), text var(--text-disabled), cursor not-allowed

All buttons:
  border-radius: 6px
  padding: 8px 16px
  font: Inter 13px weight 600
  transition: 150ms ease
```

---

## 6. FULL SCREEN SPECIFICATIONS

### Screen 6.1 — Executive Dashboard (Landing, after analysis)
```
NAV BAR (56px)
├── EXECUTIVE HEADER SECTION
│   ├── Organism name + taxonomy summary (2 rows)
│   └── 6 METRIC CARDS (responsive grid)
│
├── MAIN DASHBOARD GRID (2 columns on desktop)
│   ├── LEFT (60% width)
│   │   ├── TAXONOMY TREE (full height PyVis panel)
│   │   └── GENOME MAP (Plotly, full width)
│   │
│   └── RIGHT (40% width)
│       ├── PROTEIN DOMAINS TREEMAP
│       ├── VIRUS SIMILARITY CHART
│       └── PATHWAY SUNBURST
│
└── DOWNLOAD CENTER (footer bar, 48px)
```

### Screen 6.2 — Upload Page
```
NAV BAR
BREADCRUMB: Home > New Analysis
UPLOAD ZONE (centered, 60% width max 600px)
WORKFLOW AUTO-DETECT PREVIEW (shown after file selected)
CONFIGURATION ACCORDION (collapsed by default)
RUN BUTTON (full width, Primary style)
```

### Screen 6.3 — Progress Page
```
NAV BAR
RUN HEADER: Run ID, file, workflow, AI provider
PROGRESS STEPS TABLE (icon + name + bar + status + detail)
LIVE LOG TERMINAL (dark background, JetBrains Mono, auto-scroll, max height 200px)
ESTIMATED TIME REMAINING
CANCEL BUTTON (bottom right)
```

### Screen 6.4 — Results Pages (Tabbed)
```
NAV BAR
RUN INFO HEADER (collapsible, 40px): run_id | file | workflow | AI | runtime
TAB BAR (horizontal scrollable)
TAB CONTENT (full remaining height, scrollable)
QUICK ACTION BAR (bottom, 40px sticky): [Download All] [New Analysis] [Share Run ID]
```

---

## 7. RESPONSIVE DESIGN

| Breakpoint | Layout |
|---|---|
| ≥1440px (Desktop XL) | Full sidebar + 2-col dashboard + all panels visible |
| 1280–1439px (Desktop) | Full sidebar + 2-col dashboard |
| 1024–1279px (Laptop) | Collapsed sidebar (icon only) + 2-col dashboard |
| 768–1023px (Tablet) | Hidden sidebar (hamburger menu) + 1-col stacked layout |
| <768px (Mobile) | Bottom nav tabs + cards stack vertically + simplified charts |

**Mobile-specific adaptations:**
- Dashboard metric cards: 2x3 grid instead of 6x1
- Taxonomy tree: replaced by collapsible text lineage table (PyVis does not render well on mobile)
- Charts: Plotly mobile-optimized layout
- Tables: horizontal scroll enabled
- Upload: native file picker (no drag-and-drop)

---

## 8. ACCESSIBILITY STANDARDS

- All text: WCAG AA minimum contrast ratio (4.5:1)
- Interactive elements: 44×44px minimum touch target
- Charts: described with alt-text summary; downloadable as CSV for screen readers
- Tables: proper ARIA scope, role="grid" for sortable tables
- Color: never used as sole indicator (always accompanied by icon or text)
- Keyboard: all actions accessible via Tab/Enter/Arrow keys
- Focus states: visible 2px var(--accent-blue) outline on all interactive elements
- Reduced motion: all CSS animations respect prefers-reduced-motion
- Language: lang="en" on html element

---

## 9. ICON SYSTEM

**Icon library:** Heroicons (MIT license)
**Icon size:** 16px (inline text), 20px (buttons), 24px (nav items), 32px (empty states)

Key icon assignments:
```
🔬 Platform Logo:  custom SVG — DNA double helix in --accent-blue
🧬 FASTA/Genomics: Heroicons BeakerIcon
📊 Expression:     Heroicons ChartBarIcon
🌳 Taxonomy:       Heroicons Square3Stack3DIcon (hierarchical)
🔍 Annotation:     Heroicons MagnifyingGlassIcon
📄 Report:         Heroicons DocumentTextIcon
⬇ Download:       Heroicons ArrowDownTrayIcon
▶ Run:             Heroicons PlayIcon
✓ Complete:        Heroicons CheckCircleIcon
✗ Error:           Heroicons XCircleIcon
⚠ Warning:        Heroicons ExclamationTriangleIcon
🤖 AI:             Heroicons SparklesIcon
⚙ Settings:       Heroicons Cog6ToothIcon
🌙/☀ Theme:       Heroicons MoonIcon / SunIcon
```

---

## 10. ANIMATIONS AND TRANSITIONS

**Policy:** Animations must serve a purpose. No decorative animations.

| Trigger | Animation | Duration |
|---|---|---|
| Page load | Fade-in content area (opacity 0→1) | 200ms |
| Analysis step completes | Step row flashes accent-green once | 400ms |
| Chart renders | Plotly built-in enter animation | 300ms |
| Card hover | translateY(-2px) + border brighten | 150ms |
| Tab switch | Fade content (opacity 0→1) | 150ms |
| Progress bar fill | Width transition linear | 300ms |
| Error message | Slide in from top + shake once | 250ms |
| Upload drop zone | Border + background transition | 150ms |
| Taxonomy tree load | PyVis built-in stabilization | 500ms max |

**Reduced motion:** When prefers-reduced-motion is set, all transitions are instant (duration: 0ms) except Plotly animations which are suppressed via `transition: false` in Plotly layout.

---

## 11. EMPTY STATES (Every Panel)

Each data panel has a specific empty state message rather than a generic "No data" label.

| Panel | Empty State Message |
|---|---|
| Taxonomy Tree | "Taxonomy unavailable. Possible reason: no annotation hits, or offline mode without cache. Other analysis is complete." |
| Genome Map | "No ORFs passed the minimum length filter. Try reducing min_orf_length in config.yaml." |
| Volcano Plot | "No genes in expression matrix. Check input file format." |
| Enrichment Chart | "No significant pathways found (FDR < 0.05). Possible reason: DEG list too small or no overlap with pathway gene sets." |
| AI Report | "AI interpretation unavailable. Set [PROVIDER]_API_KEY in .env file and re-run." |
| Downloads | "Run in progress or failed. Downloads available when analysis completes." |

---

## 12. PRINT / EXPORT STYLES

A separate CSS media query `@media print` ensures the HTML report is print-friendly:
- Hide navigation bars, tab bars, and action buttons
- Charts: display at full width without interactive overlays
- Tables: break across pages cleanly (page-break-inside: avoid on rows)
- Colors: light background for printing (force white background, dark text)
- Font: reduced to 11px body, 13px headings
- Page headers/footers: show run ID, page number, timestamp

---

## 13. FULL CHANGE LOG FOR THIS REVISION

| Section Changed | Change | Reason |
|---|---|---|
| Color palette | Extended with --accent-gold, --accent-teal for taxonomy | Taxonomy tree needs distinct colors for lineage path |
| Theme system | Added light mode alongside dark mode | Accessibility; some users prefer light mode for printing |
| Executive Dashboard layout | New 2-column grid with 6 metric cards | NCBI Virus Dashboard reference design |
| Taxonomy Tree component | Full specification added | New Module 10 requirement |
| Genome Map component | Full specification added | New dashboard feature |
| Protein Domain Treemap | Full specification added | New dashboard feature |
| Virus Similarity Dashboard | Full specification added | New dashboard feature |
| Nav bar | Added AI provider indicator | User needs to know which AI provider is active |
| Upload zone | Added [Load Example: SARS-CoV-2] button | Quick demo for viva without file upload |
| Progress page | Added taxonomy and dashboard steps | New pipeline steps |
| Downloads tab | Added GFF3, taxonomy JSON, dashboard HTML | New output types |
| Metric card component | Formally specified | Reusable across dashboard |
| Empty states | Documented per panel | Better UX on failures |

---

*End of Document 4 — UI/UX Design Brief v2.1*
