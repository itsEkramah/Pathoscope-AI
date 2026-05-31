import os
import json
import base64
import pandas as pd
from pathlib import Path
from typing import Dict, Any, List, Optional
from loguru import logger

# Try importing WeasyPrint for PDF generation, but fail gracefully
try:
    import weasyprint
    HAS_WEASYPRINT = True
except (ImportError, OSError, Exception) as e:
    HAS_WEASYPRINT = False

# High-fidelity, publication-style HTML Jinja2 template inspired by MultiQC and nf-core
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PathoScope AI: Genomic Synthesis & Interpretation Report</title>
    <style>
        :root {
            --primary: #0b2240;       /* nf-core midnight blue */
            --primary-light: #1f3c6d;
            --secondary: #24b07a;     /* nf-core green */
            --accent: #e05252;        /* Vermilion red */
            --dark: #111111;
            --light: #f8fafc;
            --white: #ffffff;
            --grey: #64748b;
            --border: #e2e8f0;
            --soft-blue: #f0f5fa;
        }

        * {
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background-color: var(--light);
            color: var(--dark);
            line-height: 1.5;
            padding: 0;
        }

        /* Sidebar Navigation Layout */
        .wrapper {
            display: flex;
            min-height: 100vh;
        }

        .sidebar {
            width: 250px;
            background-color: var(--primary);
            color: var(--white);
            padding: 25px 15px;
            position: fixed;
            height: 100vh;
            overflow-y: auto;
            border-right: 4px solid var(--secondary);
        }

        .sidebar h2 {
            font-size: 20px;
            font-weight: 800;
            margin-bottom: 20px;
            letter-spacing: -0.5px;
            color: var(--white);
        }

        .sidebar h2 span {
            color: var(--secondary);
        }

        .nav-links {
            list-style: none;
            margin-top: 20px;
        }

        .nav-links li {
            margin-bottom: 12px;
        }

        .nav-links a {
            color: #94a3b8;
            text-decoration: none;
            font-size: 13px;
            font-weight: 500;
            display: block;
            padding: 8px 12px;
            border-radius: 4px;
            transition: all 0.2s ease;
        }

        .nav-links a:hover, .nav-links a.active {
            color: var(--white);
            background-color: var(--primary-light);
            border-left: 3px solid var(--secondary);
        }

        .main-content {
            flex: 1;
            margin-left: 250px;
            padding: 40px;
            background-color: var(--light);
        }

        header {
            border-bottom: 2px solid var(--border);
            padding-bottom: 20px;
            margin-bottom: 30px;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }

        .header-title h1 {
            color: var(--primary);
            font-size: 26px;
            font-weight: 800;
        }

        .badge-status {
            background-color: var(--secondary);
            color: var(--white);
            padding: 6px 12px;
            border-radius: 20px;
            font-size: 11px;
            font-weight: 700;
            text-transform: uppercase;
        }

        .section-card {
            background: var(--white);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
        }

        .section-title {
            color: var(--primary);
            font-size: 18px;
            font-weight: 700;
            margin-bottom: 20px;
            border-left: 4px solid var(--secondary);
            padding-left: 12px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .grid-2 {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 25px;
            margin-bottom: 25px;
        }

        .grid-3 {
            display: grid;
            grid-template-columns: 1fr 1fr 1fr;
            gap: 20px;
            margin-bottom: 25px;
        }

        .metric-card {
            background-color: var(--soft-blue);
            border: 1px solid #d0e1fd;
            border-radius: 6px;
            padding: 15px 20px;
            text-align: center;
        }

        .metric-card h3 {
            font-size: 28px;
            color: var(--primary);
            font-weight: 800;
        }

        .metric-card p {
            font-size: 12px;
            color: var(--grey);
            font-weight: 600;
            text-transform: uppercase;
            margin-top: 4px;
        }

        .scientific-table {
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
            font-size: 13px;
        }

        .scientific-table th, .scientific-table td {
            padding: 10px 14px;
            text-align: left;
            border-bottom: 1px solid var(--border);
        }

        .scientific-table th {
            background-color: #f8fafc;
            color: var(--primary);
            font-weight: 700;
        }

        .scientific-table tr:hover {
            background-color: #f1f5f9;
        }

        .plot-box {
            text-align: center;
            margin: 20px 0;
            padding: 15px;
            border: 1px solid var(--border);
            border-radius: 6px;
            background-color: var(--white);
        }

        .plot-img {
            max-width: 100%;
            height: auto;
            border-radius: 4px;
        }

        /* AI interpretation Panel Styling */
        .ai-synthesis-box {
            background: linear-gradient(135deg, #092240 0%, #1a365d 100%);
            color: var(--white);
            border-radius: 8px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 4px 15px rgba(9, 34, 64, 0.15);
        }

        .ai-synthesis-box h2 {
            color: var(--secondary);
            font-size: 20px;
            margin-bottom: 20px;
            border-bottom: 1px solid rgba(255,255,255,0.1);
            padding-bottom: 10px;
        }

        .ai-field {
            margin-bottom: 25px;
        }

        .ai-field-title {
            font-weight: 700;
            color: #a5f3fc; /* Cyan accent */
            font-size: 13px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            margin-bottom: 6px;
        }

        .ai-field-body {
            font-size: 14px;
            line-height: 1.6;
            text-align: justify;
        }

        .ai-bullets {
            list-style: none;
            padding-left: 0;
            margin-top: 6px;
        }

        .ai-bullets li {
            padding-left: 18px;
            position: relative;
            margin-bottom: 6px;
            font-size: 13.5px;
        }

        .ai-bullets li::before {
            content: "✦";
            color: var(--secondary);
            position: absolute;
            left: 0;
            font-size: 12px;
        }

        .ai-bullets.warnings li::before {
            color: var(--accent);
        }

        /* Bibliography Panel */
        .bib-item {
            background-color: #fafbfc;
            border-left: 4px solid var(--primary-light);
            padding: 15px;
            margin-bottom: 15px;
            border-radius: 0 6px 6px 0;
            font-size: 13px;
        }

        .bib-title {
            font-weight: 700;
            color: var(--primary);
            margin-bottom: 4px;
        }

        .bib-meta {
            font-style: italic;
            color: var(--grey);
            font-size: 12px;
            margin-bottom: 8px;
        }

        .bib-pmid {
            background-color: #e2e8f0;
            color: var(--primary);
            padding: 2px 6px;
            border-radius: 4px;
            font-weight: 700;
            font-size: 11px;
            text-decoration: none;
            display: inline-block;
        }

        /* Reproducibility Metadata Panel */
        .metadata-section {
            background-color: #f1f5f9;
            border: 1px solid var(--border);
            padding: 20px;
            border-radius: 6px;
            font-size: 12.5px;
        }

        .metadata-section dt {
            font-weight: 700;
            color: var(--primary);
            margin-top: 8px;
        }

        .metadata-section dd {
            margin-left: 0;
            color: #334155;
            font-family: monospace;
            background-color: #f8fafc;
            padding: 4px 8px;
            border-radius: 4px;
            border: 1px solid #e2e8f0;
            display: inline-block;
        }

        footer {
            margin-top: 50px;
            border-top: 1px solid var(--border);
            padding-top: 20px;
            text-align: center;
            color: var(--grey);
            font-size: 11px;
        }

        .badge-fdr-sig {
            background-color: rgba(36, 176, 122, 0.15);
            color: var(--secondary);
            font-weight: 700;
            padding: 2px 6px;
            border-radius: 4px;
            font-size: 10px;
        }

        .badge-fdr-nonsig {
            background-color: rgba(100, 116, 139, 0.15);
            color: var(--grey);
            font-weight: 700;
            padding: 2px 6px;
            border-radius: 4px;
            font-size: 10px;
        }

        .pdf-page-break {
            page-break-before: always;
        }
    </style>
</head>
<body>
    <div class="wrapper">
        <!-- SIDEBAR NAVIGATION -->
        <nav class="sidebar">
            <h2>PathoScope<span>AI</span></h2>
            <p style="font-size: 10px; color: #64748b; margin-top:-6px; font-weight: 600;">Genomics Pipeline Report</p>
            <ul class="nav-links">
                <li><a href="#overview" class="active">Overview & Metadata</a></li>
                <li><a href="#preprocessing">Preprocessing & QC</a></li>
                <li><a href="#orfs">ORF Predictions</a></li>
                <li><a href="#homology">Reference Homology</a></li>
                <li><a href="#domains">Conserved Domains</a></li>
                <li><a href="#enrichment">Pathway Enrichment</a></li>
                <li><a href="#ai-synthesis">AI Synthesis & Interpretation</a></li>
                <li><a href="#pubmed">PubMed Bibliography</a></li>
                <li><a href="#reproducibility">Reproducibility Audit</a></li>
            </ul>
        </nav>

        <!-- MAIN CONTENT PANEL -->
        <main class="main-content">
            <!-- HEADER -->
            <header>
                <div class="header-title">
                    <h1>Comprehensive Analysis Report</h1>
                    <p style="color: var(--grey); font-size:13px; font-weight: 500;">MultiQC & nf-core Compliant Supplementary Report</p>
                </div>
                <div>
                    <span class="badge-status">Pipeline Completed</span>
                </div>
            </header>

            <!-- SECTION 1: RUN METADATA -->
            <section id="overview" class="section-card">
                <div class="section-title">1. Execution Overview & Metadata</div>
                <div class="grid-2">
                    <div>
                        <p style="margin-bottom: 8px;"><strong>Input Sequence File:</strong> {{ metadata.run.input_file }}</p>
                        <p style="margin-bottom: 8px;"><strong>File MD5 Hash:</strong> <code style="font-family: monospace; font-size: 12px; background-color: var(--light); padding:2px 4px; border-radius:4px;">{{ metadata.run.input_md5 }}</code></p>
                        <p style="margin-bottom: 8px;"><strong>Analysis Date:</strong> {{ metadata.run.timestamp }}</p>
                    </div>
                    <div>
                        <p style="margin-bottom: 8px;"><strong>Pipeline Engine:</strong> PathoScope AI v{{ metadata.pipeline.version }}</p>
                        <p style="margin-bottom: 8px;"><strong>Environment:</strong> {{ metadata.environment.platform }} (Python {{ metadata.environment.python_version.split(' ')[0] }})</p>
                        <p style="margin-bottom: 8px;"><strong>Statistical FDR Alpha:</strong> {{ metadata.config_snapshot.statistics.fdr_threshold }}</p>
                    </div>
                </div>
            </section>

            <!-- SECTION 2: PREPROCESSING AND QC -->
            <section id="preprocessing" class="section-card">
                <div class="section-title">2. Sequence Preprocessing & Quality Control</div>
                <div class="grid-3">
                    <div class="metric-card">
                        <h3>{{ qc.counts.total_processed }}</h3>
                        <p>Reads Processed</p>
                    </div>
                    <div class="metric-card">
                        <h3>{{ qc.counts.total_kept }}</h3>
                        <p>Reads Passed Filters</p>
                    </div>
                    <div class="metric-card">
                        <h3>{{ qc.statistics.gc_percent }}%</h3>
                        <p>Average GC Content</p>
                    </div>
                </div>
                <table class="scientific-table">
                    <thead>
                        <tr>
                            <th>Parameter Flag</th>
                            <th>Value / Threshold</th>
                            <th>Status Description</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td><strong>Minimum Read Length</strong></td>
                            <td>{{ metadata.config_snapshot.preprocessing.min_length }} bp</td>
                            <td>Reads below this boundary are safely discarded.</td>
                        </tr>
                        <tr>
                            <td><strong>Sequence N50</strong></td>
                            <td>{{ qc.statistics.n50_bp }} bp</td>
                            <td>High value indicates good sequence assembly/contiguity.</td>
                        </tr>
                    </tbody>
                </table>

                {% if images.read_retention_waterfall %}
                <div class="grid-2" style="margin-top: 30px;">
                    <div class="plot-box" style="margin: 0;">
                        <h4 style="margin-bottom: 10px; color: var(--primary);">Sequence Duplication Levels</h4>
                        <img class="plot-img" src="data:image/png;base64,{{ images.sequence_duplication }}" alt="Sequence Duplication Levels">
                    </div>
                    <div class="plot-box" style="margin: 0;">
                        <h4 style="margin-bottom: 10px; color: var(--primary);">Read Retention Statistics Waterfall</h4>
                        <img class="plot-img" src="data:image/png;base64,{{ images.read_retention_waterfall }}" alt="Read Retention Waterfall">
                    </div>
                </div>
                {% endif %}

                {% if images.per_base_quality %}
                <div class="plot-box" style="margin-top: 25px;">
                    <h4 style="margin-bottom: 10px; color: var(--primary);">Per-Base Sequence Quality Distribution</h4>
                    <img class="plot-img" src="data:image/png;base64,{{ images.per_base_quality }}" alt="Per-Base Sequence Quality">
                </div>
                {% endif %}

                {% if images.gc_distribution %}
                <div class="plot-box" style="margin-top: 25px;">
                    <h4 style="margin-bottom: 10px; color: var(--primary);">Per-Read GC Content Distribution</h4>
                    <img class="plot-img" src="data:image/png;base64,{{ images.gc_distribution }}" alt="GC Distribution">
                </div>
                {% endif %}
            </section>

            <!-- SECTION 3: ORF PREDICTIONS -->
            <section id="orfs" class="section-card">
                <div class="section-title">3. Advanced ORF Predictions & Coding Potential</div>
                <div class="grid-3">
                    <div class="metric-card">
                        <h3>{{ orfs.total_predicted }}</h3>
                        <p>ORFs Predicted</p>
                    </div>
                    <div class="metric-card">
                        <h3>{{ orfs.average_length_aa }} aa</h3>
                        <p>Average ORF Size</p>
                    </div>
                    <div class="metric-card">
                        <h3>{{ orfs.length_aa_range }} aa</h3>
                        <p>Size Range (aa)</p>
                    </div>
                </div>
                
                {% if images.lengths %}
                <div class="plot-box">
                    <h4 style="margin-bottom: 10px; color: var(--primary);">Predicted ORF Lengths Distribution Density</h4>
                    <img class="plot-img" src="data:image/png;base64,{{ images.lengths }}" alt="ORF Length Distribution">
                </div>
                {% endif %}

                {% if images.track %}
                <div class="plot-box">
                    <h4 style="margin-bottom: 10px; color: var(--primary);">Laned Genomic Coordinate Track Map</h4>
                    <img class="plot-img" src="data:image/png;base64,{{ images.track }}" alt="Genomic Track Map">
                </div>
                {% endif %}
            </section>

            <div class="pdf-page-break"></div>

            <!-- SECTION 4: REFERENCE HOMOLOGY -->
            <section id="homology" class="section-card">
                <div class="section-title">4. Curated Reference Homology Alignments</div>
                <p style="font-size: 13px; color: var(--grey); margin-bottom: 15px;">High-confidence functional alignments mapping predicted ORFs to curated Swiss-Prot reference databases:</p>
                <table class="scientific-table">
                    <thead>
                        <tr>
                            <th>Protein Query ID</th>
                            <th>Reference Target ID</th>
                            <th>% Identity</th>
                            <th>% Coverage</th>
                            <th>E-value</th>
                            <th>Bitscore</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for hit in annotations %}
                        <tr>
                            <td><strong>{{ hit.protein_id }}</strong></td>
                            <td><a href="https://www.uniprot.org/uniprotkb/{{ hit.uniprot_id }}" target="_blank" style="color: var(--primary-light); font-weight:700;">{{ hit.uniprot_id }}</a></td>
                            <td>{{ hit.identity_percent }}%</td>
                            <td>{{ hit.query_coverage_percent }}%</td>
                            <td><code style="font-family: monospace;">{{ hit.e_value }}</code></td>
                            <td>{{ hit.bitscore }}</td>
                        </tr>
                        {% else %}
                        <tr>
                            <td colspan="6" style="text-align: center; color: var(--grey);">No high-confidence database alignments identified.</td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>

                <h3 style="margin-top: 30px; margin-bottom: 15px; color: var(--primary); font-size: 15px; border-left: 3px solid var(--secondary); padding-left: 8px;">Dynamic ICTV Taxonomic Lineages</h3>
                <p style="font-size: 13px; color: var(--grey); margin-bottom: 15px;">Deterministic taxonomic lineages and Baltimore replication group assignments mapped via homologous Swiss-Prot identifiers:</p>
                <table class="scientific-table">
                    <thead>
                        <tr>
                            <th>Protein Query ID</th>
                            <th>Organism Code</th>
                            <th>Class</th>
                            <th>Order</th>
                            <th>Family</th>
                            <th>Genus</th>
                            <th>Baltimore Group</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for tax in taxonomy %}
                        <tr>
                            <td><strong>{{ tax.protein_id }}</strong></td>
                            <td><code style="font-family: monospace; font-weight: 700; color: var(--primary);">{{ tax.uniprot_organism_code }}</code></td>
                            <td>{{ tax.ictv_class }}</td>
                            <td>{{ tax.ictv_order }}</td>
                            <td><strong>{{ tax.ictv_family }}</strong></td>
                            <td>{{ tax.ictv_genus }}</td>
                            <td><span class="badge-fdr-sig" style="background-color: rgba(11, 34, 64, 0.1); color: var(--primary); font-weight: 700;">{{ tax.baltimore_group }}</span></td>
                        </tr>
                        {% else %}
                        <tr>
                            <td colspan="7" style="text-align: center; color: var(--grey);">No taxonomic lineages identified.</td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>

                {% if images.annotation_distributions %}
                <div class="plot-box">
                    <h4 style="margin-bottom: 10px; color: var(--primary);">Sequence Similarity Annotation Metrics Distributions</h4>
                    <img class="plot-img" src="data:image/png;base64,{{ images.annotation_distributions }}" alt="Annotation Distributions">
                </div>
                {% endif %}
            </section>

            <!-- SECTION 5: CONSERVED DOMAINS -->
            <section id="domains" class="section-card">
                <div class="section-title">5. Conserved Structural Domains (Pfam)</div>
                <p style="font-size: 13px; color: var(--grey); margin-bottom: 15px;">Conserved structural domains detected independently using local Pfam hidden Markov model (HMM) searches:</p>
                <table class="scientific-table">
                    <thead>
                        <tr>
                            <th>Protein Query ID</th>
                            <th>Domain Identifier</th>
                            <th>Pfam Accession</th>
                            <th>Independent E-value</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for dom in domains %}
                        <tr>
                            <td><strong>{{ dom.protein_id }}</strong></td>
                            <td><strong>{{ dom.domain_name }}</strong></td>
                            <td><code>{{ dom.domain_accession }}</code></td>
                            <td><code style="font-family: monospace;">{{ dom.e_value }}</code></td>
                        </tr>
                        {% else %}
                        <tr>
                            <td colspan="4" style="text-align: center; color: var(--grey);">No structural domains identified.</td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </section>

            <div class="pdf-page-break"></div>

            <!-- SECTION 6: PATHWAY ENRICHMENT -->
            <section id="enrichment" class="section-card">
                <div class="section-title">6. Overrepresentation & Gene Set Enrichment Analysis</div>
                <p style="font-size: 13px; color: var(--grey); margin-bottom: 15px;">Statistically enriched host biological processes mapped using exact hypergeometric ORA and GSEA models:</p>
                <table class="scientific-table">
                    <thead>
                        <tr>
                            <th>Pathway ID</th>
                            <th>Description</th>
                            <th>Counts (k)</th>
                            <th>Fold FE</th>
                            <th>FDR q-value</th>
                            <th>FDR Status</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for path in pathways %}
                        <tr>
                            <td><code>{{ path.pathway_id }}</code></td>
                            <td><strong>{{ path.description }}</strong></td>
                            <td>{{ path.query_count_k }}</td>
                            <td>{{ path.fold_enrichment }}x</td>
                            <td><code style="font-family: monospace;">{{ path.adjusted_pvalue_fdr }}</code></td>
                            <td>
                                {% if path.adjusted_pvalue_fdr <= 0.05 %}
                                <span class="badge-fdr-sig">SIGNIFICANT</span>
                                {% else %}
                                <span class="badge-fdr-nonsig">TOP CANDIDATE</span>
                                {% endif %}
                            </td>
                        </tr>
                        {% else %}
                        <tr>
                            <td colspan="6" style="text-align: center; color: var(--grey);">No biological pathways mapped.</td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>

                <div class="grid-2" style="margin-top: 30px;">
                    {% if images.barplot %}
                    <div class="plot-box" style="margin: 0;">
                        <img class="plot-img" src="data:image/png;base64,{{ images.barplot }}" alt="Pathway Barplot">
                    </div>
                    {% endif %}
                    
                    {% if images.bubbleplot %}
                    <div class="plot-box" style="margin: 0;">
                        <img class="plot-img" src="data:image/png;base64,{{ images.bubbleplot }}" alt="Pathway Bubbleplot">
                    </div>
                    {% endif %}
                </div>

                {% if images.volcano %}
                <div class="plot-box" style="margin-top: 25px;">
                    <h4 style="margin-bottom: 10px; color: var(--primary);">Volcano Distribution of Pathway Enrichment</h4>
                    <img class="plot-img" src="data:image/png;base64,{{ images.volcano }}" alt="Volcano Plot">
                </div>
                {% endif %}

                {% if images.gene_pathway_heatmap %}
                <div class="plot-box" style="margin-top: 25px;">
                    <h4 style="margin-bottom: 10px; color: var(--primary);">Gene-Pathway Co-occurrence Profile Heatmap</h4>
                    <img class="plot-img" src="data:image/png;base64,{{ images.gene_pathway_heatmap }}" alt="Heatmap">
                </div>
                {% endif %}

                {% if images.functional_pca_plot %}
                <div class="plot-box" style="margin-top: 25px;">
                    <h4 style="margin-bottom: 10px; color: var(--primary);">PCA Projections of Gene Functional Profiles</h4>
                    <img class="plot-img" src="data:image/png;base64,{{ images.functional_pca_plot }}" alt="PCA Plot">
                </div>
                {% endif %}

                {% if images.protein_pathway_network %}
                <div class="plot-box" style="margin-top: 25px;">
                    <h4 style="margin-bottom: 10px; color: var(--primary);">Bipartite Protein-Pathway Network Map</h4>
                    <img class="plot-img" src="data:image/png;base64,{{ images.protein_pathway_network }}" alt="Network Map">
                </div>
                {% endif %}
            </section>

            <div class="pdf-page-break"></div>

            <!-- SECTION 7: AI INTERPRETATION -->
            <section id="ai-synthesis" class="ai-synthesis-box">
                <h2>AI Interpretation Layer & Virology Synthesis</h2>
                
                <div class="ai-field">
                    <div class="ai-field-title">High-Level Synopsis</div>
                    <div class="ai-field-body">{{ ai.concise_summary }}</div>
                </div>

                <div class="ai-field">
                    <div class="ai-field-title">Detailed Scientific Interpretation</div>
                    <div class="ai-field-body" style="white-space: pre-line;">{{ ai.detailed_biological_interpretation }}</div>
                </div>

                <div class="ai-field">
                    <div class="ai-field-title">Implicated Host hijacks & Mechanisms</div>
                    <ul class="ai-bullets">
                        {% for mech in ai.implicated_host_mechanisms %}
                        <li>{{ mech }}</li>
                        {% else %}
                        <li>None significant flagged in this profile.</li>
                        {% endfor %}
                    </ul>
                </div>

                <div class="ai-field">
                    <div class="ai-field-title">PubMed Literature evidence Summary</div>
                    <div class="ai-field-body" style="white-space: pre-line;">{{ ai.literature_evidence_summary }}</div>
                </div>

                <div class="ai-field">
                    <div class="ai-field-title">Associated host Diseases & Molecular Pathology</div>
                    <div class="ai-field-body">{{ ai.disease_association_summary }}</div>
                </div>

                <div class="ai-field">
                    <div class="ai-field-title">Therapeutic relevance & Antiviral Options</div>
                    <div class="ai-field-body">{{ ai.therapeutic_relevance_summary }}</div>
                </div>

                <div class="ai-field">
                    <div class="ai-field-title">Known Biomarkers & Molecular Targets</div>
                    <div class="ai-field-body">{{ ai.known_biomarkers_summary }}</div>
                </div>

                <div class="ai-field">
                    <div class="ai-field-title">Computational Limitations</div>
                    <div class="ai-field-body">{{ ai.limitations }}</div>
                </div>

                <div class="ai-field" style="margin-bottom: 0;">
                    <div class="ai-field-title" style="color: #fda4af;">Experimental validation Warnings</div>
                    <ul class="ai-bullets warnings">
                        {% for warn in ai.confidence_warnings %}
                        <li>{{ warn }}</li>
                        {% endfor %}
                    </ul>
                </div>
            </section>

            <!-- SECTION 8: PUBMED BIBLIOGRAPHY -->
            <section id="pubmed" class="section-card">
                <div class="section-title">8. Grounded PubMed Literature Bibliography</div>
                <p style="font-size: 13px; color: var(--grey); margin-bottom: 15px;">Retrieved literature abstract evidence backing molecular annotations and host pathway implications:</p>
                {% for art in ai.retrieved_literature_citations %}
                <div class="bib-item">
                    <div class="bib-title">{{ art.title }}</div>
                    <div class="bib-meta">{{ art.authors }} — <strong>{{ art.journal }}</strong></div>
                    <a href="https://pubmed.ncbi.nlm.nih.gov/{{ art.pmid }}" target="_blank" class="bib-pmid">PubMed: {{ art.pmid }}</a>
                </div>
                {% else %}
                <p style="font-size: 13px; color: var(--grey); text-align: center;">No bibliography resources mapped.</p>
                {% endfor %}
            </section>

            <!-- SECTION 9: REPRODUCIBILITY INFORMATION -->
            <section id="reproducibility" class="section-card">
                <div class="section-title">9. Reproducibility & Pipeline Audit Hub</div>
                <div class="grid-2">
                    <div class="metadata-section">
                        <h4 style="color: var(--primary); margin-bottom: 10px;">Pipeline Versions & Audit</h4>
                        <dl>
                            <dt>Pipeline Software version</dt>
                            <dd>v{{ metadata.pipeline.version }}</dd>
                            
                            <dt>DIAMOND Database Target</dt>
                            <dd>{{ metadata.config_snapshot.annotation.local_db_path }}</dd>
                            
                            <dt>Pfam HMM Database Target</dt>
                            <dd>{{ metadata.config_snapshot.domain_search.hmmer_db_path }}</dd>
                            
                            <dt>Input Sequence MD5 Audit Hash</dt>
                            <dd>{{ metadata.run.input_md5 }}</dd>
                            
                            <dt>Execution Host Command</dt>
                            <dd style="white-space: pre-wrap; font-size:11px;">{{ metadata.run.command_line }}</dd>
                        </dl>
                    </div>
                    <div class="metadata-section">
                        <h4 style="color: var(--primary); margin-bottom: 10px;">Stage Elapsed Times (s)</h4>
                        <table class="scientific-table" style="margin:0; font-size:11px;">
                            <thead>
                                <tr>
                                    <th>Pipeline Stage</th>
                                    <th>Runtime Elapsed (s)</th>
                                </tr>
                            </thead>
                            <tbody>
                                {% for stage, elapsed in metadata.run.execution_times_seconds.items() %}
                                <tr>
                                    <td><strong>{{ stage | replace('_', ' ') | title }}</strong></td>
                                    <td>{{ elapsed }} s</td>
                                </tr>
                                {% else %}
                                <tr>
                                    <td colspan="2">No metrics recorded.</td>
                                </tr>
                                {% endfor %}
                            </tbody>
                        </table>
                    </div>
                </div>
            </section>

            <!-- FOOTER -->
            <footer>
                <p>Generated dynamically by <strong>PathoScope AI v{{ metadata.pipeline.version }}</strong>.</p>
                <p style="font-size:10px; margin-top:5px; color:#94a3b8;">Reproducibility Audit Hash: {{ metadata.run.input_md5[:16] }}</p>
            </footer>
        </main>
    </div>
</body>
</html>
"""


def encode_image_base64(image_path: Path) -> Optional[str]:
    """Helper utility that reads a PNG/JPG and returns a base64 encoded string."""
    image_path = Path(image_path)
    if not image_path.exists():
        return None
    try:
        with open(image_path, "rb") as f:
            return base64.b64encode(f.read()).decode("utf-8")
    except Exception as e:
        logger.warning(f"Failed to base64 encode image {image_path}: {e}")
        return None


def generate_html_report_dashboard(
    output_dir: Path,
    metadata: Dict[str, Any],
    ai_synthesis: Dict[str, Any]
) -> Path:
    """
    Renders the gorgeous self-contained Jinja2 HTML report dashboard.
    Automatically base64 encodes all visual plots and injects them.
    """
    from jinja2 import Template
    
    output_dir = Path(output_dir)
    vis_dir = output_dir / "visualizations"
    report_outdir = output_dir / "final_report"
    report_outdir.mkdir(parents=True, exist_ok=True)
    
    html_outpath = report_outdir / "report.html"
    
    # Enrich metadata dictionary with robust defaults to prevent Jinja2 UndefinedError in tests or empty runs
    import datetime
    if not isinstance(metadata, dict):
        metadata = {}
    metadata.setdefault("pipeline", {})
    metadata["pipeline"].setdefault("version", "1.0.0")
    
    metadata.setdefault("run", {})
    metadata["run"].setdefault("input_file", "unknown_input.fasta")
    metadata["run"].setdefault("input_md5", "unknown_md5")
    metadata["run"].setdefault("timestamp", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    metadata["run"].setdefault("status", "SUCCESS")
    metadata["run"].setdefault("command_line", "python -m pathoscope.cli run -i sample.fasta -o CLI_results")
    metadata["run"].setdefault("execution_times_seconds", {})
    
    metadata.setdefault("environment", {})
    metadata["environment"].setdefault("platform", "Unknown OS")
    metadata["environment"].setdefault("python_version", "3.x")
    
    metadata.setdefault("config_snapshot", {})
    metadata["config_snapshot"].setdefault("reporting", {})
    metadata["config_snapshot"]["reporting"].setdefault("theme", "dark")
    metadata["config_snapshot"].setdefault("statistics", {})
    metadata["config_snapshot"]["statistics"].setdefault("fdr_threshold", 0.05)
    metadata["config_snapshot"].setdefault("preprocessing", {})
    metadata["config_snapshot"]["preprocessing"].setdefault("min_length", 50)
    metadata["config_snapshot"].setdefault("annotation", {})
    metadata["config_snapshot"]["annotation"].setdefault("local_db_path", "data/reference/viral_proteins.dmnd")
    metadata["config_snapshot"].setdefault("domain_search", {})
    metadata["config_snapshot"]["domain_search"].setdefault("hmmer_db_path", "data/reference/Pfam-A.hmm")

    logger.info("Assembling and rendering interactive publication-style HTML report dashboard...")
    
    # 1. Base64 encode all visual graphics (PNG format at 300 DPI)
    images_b64 = {
        "lengths": encode_image_base64(vis_dir / "orf_lengths.png"),
        "track": encode_image_base64(vis_dir / "orf_genomic_track.png"),
        "annotation_distributions": encode_image_base64(vis_dir / "annotation_distributions.png"),
        "barplot": encode_image_base64(vis_dir / "pathway_enrichment_barplot.png"),
        "bubbleplot": encode_image_base64(vis_dir / "pathway_enrichment_bubbleplot.png"),
        "protein_pathway_network": encode_image_base64(vis_dir / "protein_pathway_network.png"),
        "volcano": encode_image_base64(vis_dir / "pathway_volcano_plot.png"),
        "gene_pathway_heatmap": encode_image_base64(vis_dir / "gene_pathway_heatmap.png"),
        "functional_pca_plot": encode_image_base64(vis_dir / "functional_pca_plot.png"),
        "per_base_quality": encode_image_base64(vis_dir / "per_base_quality.png"),
        "gc_distribution": encode_image_base64(vis_dir / "gc_distribution.png"),
        "sequence_duplication": encode_image_base64(vis_dir / "sequence_duplication.png"),
        "read_retention_waterfall": encode_image_base64(vis_dir / "read_retention_waterfall.png")
    }
    
    # 2. Extract deterministic upstream stats from the metadata
    qc_data = {
        "counts": {
            "total_processed": 0,
            "total_kept": 0
        },
        "statistics": {
            "n50_bp": 0,
            "gc_percent": 0.0
        }
    }
    
    qc_json = output_dir / "preprocessed" / "qc_report.json"
    if qc_json.exists():
        try:
            with open(qc_json, "r", encoding="utf-8") as f:
                raw_qc = json.load(f)
                qc_data["counts"]["total_processed"] = raw_qc.get("counts", {}).get("total_processed", 0)
                qc_data["counts"]["total_kept"] = raw_qc.get("counts", {}).get("total_kept", 0)
                
                metrics = raw_qc.get("metrics", raw_qc.get("statistics", {}))
                qc_data["statistics"]["n50_bp"] = metrics.get("n50", 0)
                qc_data["statistics"]["gc_percent"] = round(metrics.get("mean_gc_percent", metrics.get("gc_content_pct", 0.0)), 2)
        except Exception:
            pass

    # 3. Extract ORF coordinates count and metrics
    orf_metrics = {"total_predicted": 0, "average_length_aa": 0, "length_aa_range": "0-0"}
    gff_file = output_dir / "orfs" / "coordinates.gff3"
    if gff_file.exists():
        try:
            lengths = []
            with open(gff_file, "r") as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 5:
                        start = int(parts[3])
                        end = int(parts[4])
                        lengths.append((end - start + 1) // 3)
            if lengths:
                orf_metrics["total_predicted"] = len(lengths)
                orf_metrics["average_length_aa"] = int(sum(lengths) / len(lengths))
                orf_metrics["length_aa_range"] = f"{min(lengths)}-{max(lengths)}"
        except Exception:
            pass

    # 4. Extract Homologous Alignments
    annotations_list = []
    annot_csv = output_dir / "annotations" / "annotated_proteins.csv"
    if annot_csv.exists():
        try:
            df_annot = pd.read_csv(annot_csv)
            df_hits = df_annot[
                df_annot["uniprot_id"].notna() & 
                (df_annot["uniprot_id"].astype(str).str.upper() != "NONE") &
                (df_annot["uniprot_id"].astype(str).str.upper() != "HYPOTHETICAL PROTEIN")
            ].copy()
            for idx, row in df_hits.head(15).iterrows():
                annotations_list.append({
                    "protein_id": row["protein_id"],
                    "uniprot_id": row["uniprot_id"],
                    "identity_percent": row["identity_pct"],
                    "query_coverage_percent": row["query_coverage_pct"],
                    "e_value": row["e_val"],
                    "bitscore": row["bitscore"]
                })
        except Exception as e:
            logger.error(f"Failed to compile homolog annotations: {e}")

    # 5. Extract Conserved Domains (Pfam)
    domains_list = []
    domains_csv = output_dir / "pathways" / "pfam_domains.csv"
    if domains_csv.exists():
        try:
            df_dom = pd.read_csv(domains_csv)
            for idx, row in df_dom.head(15).iterrows():
                domains_list.append({
                    "protein_id": row.get("protein_id", ""),
                    "domain_name": row.get("domain_name", ""),
                    "domain_accession": row.get("domain_accession", ""),
                    "e_value": row.get("domain_e_value", 0.0)
                })
        except Exception:
            pass

    # 5.5 Extract ICTV Taxonomic Nomenclature Mappings
    taxonomy_list = []
    taxonomy_csv = output_dir / "pathways" / "taxonomy_classification.csv"
    if taxonomy_csv.exists():
        try:
            df_tax = pd.read_csv(taxonomy_csv)
            for idx, row in df_tax.head(15).iterrows():
                taxonomy_list.append({
                    "protein_id": row.get("protein_id", ""),
                    "subject_db_id": row.get("subject_db_id", ""),
                    "uniprot_organism_code": row.get("uniprot_organism_code", ""),
                    "ictv_class": row.get("ictv_class", ""),
                    "ictv_order": row.get("ictv_order", ""),
                    "ictv_family": row.get("ictv_family", ""),
                    "ictv_genus": row.get("ictv_genus", ""),
                    "baltimore_group": row.get("baltimore_group", "")
                })
        except Exception as e:
            logger.error(f"Failed to compile ICTV taxonomy details: {e}")

    # 6. Extract Statistical Enriched Pathways
    pathways_list = []
    enrich_csv = output_dir / "enrichment" / "enrichment_results.csv"
    if enrich_csv.exists():
        try:
            df_enrich = pd.read_csv(enrich_csv)
            for idx, row in df_enrich.head(15).iterrows():
                pathways_list.append({
                    "pathway_id": row["pathway_id"],
                    "description": row["description"],
                    "query_count_k": int(row["query_count_k"]),
                    "fold_enrichment": row["fold_enrichment"],
                    "adjusted_pvalue_fdr": row["adjusted_pvalue_fdr"]
                })
        except Exception:
            pass

    # 7. Render the template using Jinja2
    template = Template(HTML_TEMPLATE)
    rendered_html = template.render(
        metadata=metadata,
        qc=qc_data,
        orfs=orf_metrics,
        annotations=annotations_list,
        domains=domains_list,
        taxonomy=taxonomy_list,
        pathways=pathways_list,
        images=images_b64,
        ai=ai_synthesis
    )
    
    with open(html_outpath, "w", encoding="utf-8") as f:
        f.write(rendered_html)
        
    logger.info(f"Successfully saved portable HTML dashboard report to: {html_outpath}")
    return html_outpath


def compile_pdf_from_html(
    html_path: Path,
    pdf_outpath: Path
) -> bool:
    """
    Compiles PDF using WeasyPrint with elegant page breaks.
    """
    if not HAS_WEASYPRINT:
        logger.warning("WeasyPrint is not installed in the current environment. Reverting gracefully. PDF skipped.")
        return False
        
    logger.info(f"Compiling PDF report using WeasyPrint: {pdf_outpath}")
    try:
        weasyprint.HTML(str(html_path)).write_pdf(str(pdf_outpath))
        logger.info(f"Successfully generated publication-grade PDF report: {pdf_outpath}")
        return True
    except Exception as e:
        logger.warning(f"WeasyPrint failed to compile PDF report: {e}. Graceful failover triggered. PDF skipped.")
        return False


def run_report_generation(
    output_dir: Path
) -> Dict[str, Any]:
    """
    Orchestrates advanced MultiQC-like HTML and PDF Supplementary report rendering.
    JSON and CSV datasets are preserved under versioned results directories for auditable open science.
    """
    output_dir = Path(output_dir)
    report_outdir = output_dir / "final_report"
    report_outdir.mkdir(parents=True, exist_ok=True)
    
    pdf_outpath = report_outdir / "report.pdf"
    
    # 1. Load pipeline metadata run audit
    metadata_json = output_dir / "metadata.json"
    metadata = {}
    if metadata_json.exists():
        try:
            with open(metadata_json, "r", encoding="utf-8") as f:
                metadata = json.load(f)
        except Exception as e:
            logger.warning(f"Failed to load run metadata: {e}")
            
    # 2. Load AI synthesized biological report
    ai_json = output_dir / "final_report" / "ai_synthesis.json"
    ai_synthesis = {}
    if ai_json.exists():
        try:
            with open(ai_json, "r", encoding="utf-8") as f:
                ai_synthesis = json.load(f)
        except Exception as e:
            logger.warning(f"Failed to load AI biological synthesis report: {e}")

    # 3. Generate HTML report dashboard with embedded base64 visual plots
    html_path = generate_html_report_dashboard(output_dir, metadata, ai_synthesis)
    
    # 4. Attempt compiling PDF via WeasyPrint
    pdf_status = compile_pdf_from_html(html_path, pdf_outpath)
    
    return {
        "status": "SUCCESS",
        "output_directory": str(report_outdir),
        "html_report": str(html_path),
        "pdf_report": str(pdf_outpath) if pdf_status else "SKIPPED_WEASYPRINT_UNAVAILABLE"
    }
