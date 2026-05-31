# pyrefly: ignore [missing-import]
import pytest
from pathlib import Path
import tempfile
import sqlite3
import pandas as pd
from unittest.mock import patch, MagicMock
from Bio import SeqIO
from pathoscope.core.pathway_mapper import (
    PfamDomain,
    PathwayCacheDB,
    parse_hmmscan_domtblout,
    query_kegg_pathways,
    query_reactome_pathways,
    get_fallback_pathways,
    parse_uniprot_id,
    process_pathway_and_domain_mapping
)
from pathoscope.utils.config_loader import AppConfig

# 1. Test Pfam domtblout line parsing
def test_pfam_domain_parsing():
    # Tabular domtblout line layout: whitespace separated
    # target name, target accession, tlen, query name, query accession, qlen, full seq evalue, full seq score, full seq bias, domain idx, domain total, c-Evalue, i-Evalue, domain score, domain bias, hmm start, hmm end, ali start, ali end, env start, env end, acc, description
    line = (
        "RNA_pol_N\tPF00561.21\t150\tORF_1\t-\t200\t1.2e-20\t80.5\t0.1\t1\t1"
        "\t2.4e-22\t1.8e-22\t85.0\t0.1\t10\t140\t15\t145\t12\t148\t0.95\tRNA-directed RNA polymerase N-terminal"
    )
    fields = line.split("\t")
    dom = PfamDomain(fields)
    
    assert dom.domain_name == "RNA_pol_N"
    assert dom.accession == "PF00561"  # strip version
    assert dom.protein_id == "ORF_1"
    assert dom.i_evalue == 1.8e-22
    assert dom.score == 85.0
    assert dom.hmm_start == 10
    assert dom.hmm_end == 140
    assert dom.ali_start == 15
    assert dom.ali_end == 145
    assert "RNA-directed RNA polymerase" in dom.description
    
    d = dom.to_dict()
    assert d["domain_name"] == "RNA_pol_N"
    assert d["alignment_start"] == 15

# 2. Test Pfam tabular parser
def test_parse_hmmscan_domtblout():
    with tempfile.NamedTemporaryFile(suffix=".domtblout", delete=False, mode="w") as f:
        # domtblout columns
        f.write("# comment line should be skipped\n")
        f.write("Pfam_domain\tPF99999.1\t100\tORF_1\t-\t100\t1e-15\t50.0\t0\t1\t1\t1e-18\t1e-18\t55.0\t0\t1\t99\t1\t99\t1\t99\t0.9\tDescription details\n")
        temp_path = Path(f.name)
        
    try:
        hits = parse_hmmscan_domtblout(temp_path)
        assert len(hits) == 1
        assert hits[0].domain_name == "Pfam_domain"
        assert hits[0].accession == "PF99999"
    finally:
        temp_path.unlink()

# 3. Test SQLite Pathway Caching
def test_sqlite_pathway_cache():
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
        temp_db_path = Path(f.name)
        
    try:
        cache = PathwayCacheDB(temp_db_path)
        
        # Verify empty get
        assert cache.get("uniprot:P12345", "kegg") is None
        
        # Set cache
        mock_data = [{"pathway_id": "map03010", "description": "Ribosome", "database": "KEGG"}]
        cache.set("uniprot:P12345", "kegg", mock_data)
        
        # Verify get hit
        retrieved = cache.get("uniprot:P12345", "kegg")
        assert retrieved is not None
        assert retrieved[0]["pathway_id"] == "map03010"
        assert retrieved[0]["description"] == "Ribosome"
        
    finally:
        if temp_db_path.exists():
            try:
                temp_db_path.unlink()
            except OSError:
                pass

# 4. Test dynamic API query and local fallbacks
@patch("requests.get")
def test_query_kegg_pathways_network_hit(mock_get):
    # Mock successful HTTP request returning KEGG pathway maps
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.text = "up:P12345\tpath:ko03010\n"
    mock_get.return_value = mock_response
    
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
        temp_db = Path(f.name)
        
    try:
        cache = PathwayCacheDB(temp_db)
        
        # Mock pathway description query as well
        with patch("pathoscope.core.pathway_mapper.fetch_kegg_pathway_description") as mock_desc:
            mock_desc.return_value = "Ribosome"
            
            paths = query_kegg_pathways("P12345", cache)
            
            assert len(paths) == 1
            assert paths[0]["pathway_id"] == "ko03010"
            assert paths[0]["description"] == "Ribosome"
            
            # Assert cache is populated
            assert cache.get("uniprot:P12345", "kegg") is not None
            
    finally:
        if temp_db.exists():
            try:
                temp_db.unlink()
            except OSError:
                pass

def test_fallback_pathways():
    # Test polymerase triggers polymerase pathway fallback offline
    paths = get_fallback_pathways("sp|P12345|RDRP_VIRUS")
    assert len(paths) == 2
    assert paths[0]["pathway_id"] == "map03010"
    
    # Non-matching subject gives empty fallback list
    assert len(get_fallback_pathways("sp|P12345|HYPOTHETICAL")) == 0

def test_parse_uniprot_id():
    assert parse_uniprot_id("sp|P12345|POL_VIRUS") == "P12345"
    assert parse_uniprot_id("P12345") == "P12345"
    assert parse_uniprot_id("None") == ""
    assert parse_uniprot_id(None) == ""

# 5. Test complete process_pathway_and_domain_mapping orchestrator
def test_process_pathway_and_domain_mapping_workflow():
    # Setup mock input files
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
        f.write(">ORF_1\nMAPKRILV\n")
        temp_proteins = Path(f.name)
        
    # Write mock annotated proteins CSV
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w") as f:
        f.write("protein_id,length_aa,annotation_status,subject_db_id,description,evalue,bitscore\n")
        f.write("ORF_1,8,Annotated,sp|P12345|POL_VIRUS,Homolog polymerase,1e-15,50.0\n")
        temp_anno = Path(f.name)
        
    temp_outdir = Path(tempfile.mkdtemp())
    temp_db_path = temp_outdir / "pathways_cache.db"
    
    # Write mock Pfam cache
    cached_dom = temp_outdir / "pfam_domains_cache.domtblout"
    with open(cached_dom, "w", encoding="utf-8") as f:
        f.write("Pfam_pol\tPF00561.21\t100\tORF_1\t-\t100\t1e-15\t50.0\t0\t1\t1\t1e-18\t1e-18\t55.0\t0\t1\t99\t1\t99\t1\t99\t0.9\tRNA-directed RNA polymerase\n")
        
    config = AppConfig()
    config.domain_search.eval_threshold = 1e-4
    config.pathway_mapping.db_cache_path = str(temp_db_path)
    config.pathway_mapping.use_kegg_api = False # disable active network calls in testing
    
    try:
        # The run will trigger fallbacks since use_kegg_api is False
        stats = process_pathway_and_domain_mapping(temp_proteins, temp_anno, temp_outdir, config)
        
        # Verify stats
        assert stats["domains"]["total_domains_detected"] == 1
        assert stats["pathways"]["total_proteins_mapped"] == 1
        
        # Verify files created
        assert (temp_outdir / "pfam_domains.csv").exists()
        assert (temp_outdir / "mapped_pathways.csv").exists()
        assert (temp_outdir / "pathway_report.json").exists()
        
        # Verify CSV table headers
        df_paths = pd.read_csv(temp_outdir / "mapped_pathways.csv")
        assert len(df_paths) >= 1
        assert df_paths.iloc[0]["protein_id"] == "ORF_1"
        assert df_paths.iloc[0]["uniprot_id"] == "P12345"
        
    finally:
        if temp_proteins.exists():
            try:
                temp_proteins.unlink()
            except OSError:
                pass
        if temp_anno.exists():
            try:
                temp_anno.unlink()
            except OSError:
                pass
        for p in temp_outdir.iterdir():
            try:
                p.unlink()
            except OSError:
                pass
        try:
            temp_outdir.rmdir()
        except OSError:
            pass


# 6. Test SQLite Cache integrity and schema checking
def test_sqlite_cache_integrity_validation():
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
        temp_db_path = Path(f.name)
        
    try:
        # Corrupt the database file
        with open(temp_db_path, "w") as f:
            f.write("corrupted garbage content")
            
        # Creating CacheDB should detect corruption and gracefully rebuild or revert to memory
        cache = PathwayCacheDB(temp_db_path)
        
        # Verify it still functions
        assert cache.get("uniprot:P12345", "kegg") is None
        cache.set("uniprot:P12345", "kegg", [{"pathway_id": "map03010", "description": "Ribosome", "database": "KEGG"}])
        assert cache.get("uniprot:P12345", "kegg") is not None
        
    finally:
        import gc
        cache = None
        gc.collect()
        if temp_db_path.exists():
            try:
                temp_db_path.unlink()
            except Exception:
                pass

# 7. Test API query retry and backoff mechanisms
@patch("requests.get")
def test_api_query_retry_mechanism(mock_get):
    # Mock a server 503 error for first 2 attempts, then success 200
    mock_503 = MagicMock()
    mock_503.status_code = 503
    
    mock_200 = MagicMock()
    mock_200.status_code = 200
    mock_200.text = "success"
    
    mock_get.side_effect = [mock_503, mock_503, mock_200]
    
    from pathoscope.core.pathway_mapper import execute_api_query_with_retry
    res = execute_api_query_with_retry("http://mock-api.com", retries=3, backoff_factor=0.01)
    
    assert res is not None
    assert res.status_code == 200
    assert mock_get.call_count == 3

# 8. Test Reactome parent stable ID traversal
@patch("pathoscope.core.pathway_mapper.execute_api_query_with_retry")
def test_reactome_hierarchy_traversal(mock_query):
    mock_res = MagicMock()
    mock_res.status_code = 200
    # Reactome parent structure matches hierarchical branches list of lists
    mock_res.json.return_value = [
        [
            {"dbId": 1, "displayName": "Top Level Category", "stId": "R-HSA-1"},
            {"dbId": 2, "displayName": "Parent Pathway", "stId": "R-HSA-2"},
            {"dbId": 3, "displayName": "Target Pathway", "stId": "R-HSA-3"}
        ]
    ]
    mock_query.return_value = mock_res
    
    temp_db = Path(":memory:")
    cache = PathwayCacheDB(temp_db)
    from pathoscope.core.pathway_mapper import query_reactome_parents
    parents = query_reactome_parents("R-HSA-3", cache)
    
    assert len(parents) == 1
    assert parents[0]["parent_id"] == "R-HSA-2"
    assert parents[0]["parent_description"] == "Parent Pathway"

# 9. Test KEGG CLASS category hierarchy parsing
@patch("pathoscope.core.pathway_mapper.execute_api_query_with_retry")
def test_kegg_class_hierarchy_parsing(mock_query):
    mock_res = MagicMock()
    mock_res.status_code = 200
    mock_res.text = (
        "ENTRY       ko03010             Pathway\n"
        "NAME        Ribosome\n"
        "CLASS       Genetic Information Processing; Translation; Ribosome\n"
    )
    mock_query.return_value = mock_res
    
    temp_db = Path(":memory:")
    cache = PathwayCacheDB(temp_db)
    from pathoscope.core.pathway_mapper import fetch_kegg_hierarchy
    cats = fetch_kegg_hierarchy("ko03010", cache)
    
    assert len(cats) == 3
    assert cats[0] == "Genetic Information Processing"
    assert cats[1] == "Translation"
    assert cats[2] == "Ribosome"

# 10. Test pathway mapping confidence scoring and Pfam boosts
def test_pathway_confidence_scoring():
    from pathoscope.core.pathway_mapper import calculate_pathway_confidence
    # Baseline average confidence without Pfam domain
    conf_base = calculate_pathway_confidence([0.8], has_pfam_validation=False)
    assert conf_base == 0.8
    
    # Pfam validated boosts confidence by 10%
    conf_boost = calculate_pathway_confidence([0.8], has_pfam_validation=True)
    assert conf_boost == 0.9
    
    # Caps score at 1.0 maximum
    conf_max = calculate_pathway_confidence([0.95], has_pfam_validation=True)
    assert conf_max == 1.0

# 11. Test pathway co-occurrence relationship network graph creation
def test_pathway_relationship_graph_creation():
    from pathoscope.core.pathway_mapper import generate_pathway_relationship_graph
    
    pathway_rows = [
        {"protein_id": "ORF_1", "pathway_id": "map03010", "pathway_description": "Ribosome"},
        {"protein_id": "ORF_1", "pathway_id": "R-HSA-2", "pathway_description": "Parent Pathway"},
        {"protein_id": "ORF_2", "pathway_id": "map03010", "pathway_description": "Ribosome"}
    ]
    
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        temp_png = Path(f.name)
        
    try:
        temp_png.unlink() # remove so we verify it is created
        generate_pathway_relationship_graph(pathway_rows, temp_png)
        assert temp_png.exists()
    finally:
        if temp_png.exists():
            temp_png.unlink()
