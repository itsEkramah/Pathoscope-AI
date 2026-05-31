"""
id_normalizer.py

Enforces HGNC-compliant gene symbol normalization. Maps Ensembl Gene IDs,
Entrez Gene IDs, and gene synonyms to official HGNC Gene Symbols.
Integrated with a local pre-populated lookup dict and an offline SQLite cache registry
with remote HGNC API fallbacks.
"""

import os
import sqlite3
import requests
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from loguru import logger
from pathoscope.core.exceptions import NormalizationError

# High-fidelity local curated mapping dictionary for core immunovirology and reference genes
# Pre-packaged to guarantee instant offline execution during vivas and test runs.
CURATED_ID_MAPS = {
    # Symbols / Synonyms mapping to HGNC official symbol
    "IL-6": "IL6", "IL_6": "IL6", "INTERLEUKIN-6": "IL6",
    "IFN-B1": "IFNB1", "IFN_B1": "IFNB1", "IFN-BETA": "IFNB1",
    "TP53": "TP53", "P53": "TP53",
    "BRCA1": "BRCA1",
    "STAT1": "STAT1",
    "MX1": "MX1", "MXA": "MX1",
    "ISG15": "ISG15", "G1P2": "ISG15",
    "CXCL10": "CXCL10", "IP-10": "CXCL10", "IP10": "CXCL10",
    "TNF": "TNF", "TNFA": "TNF", "TNF-ALPHA": "TNF",
    "IL1B": "IL1B", "IL-1B": "IL1B", "IL-1BETA": "IL1B",
    
    # Ensembl Gene IDs mapping to HGNC official symbol
    "ENSG00000136244": "IL6",
    "ENSG00000171855": "IFNB1",
    "ENSG00000141510": "TP53",
    "ENSG00000012048": "BRCA1",
    "ENSG00000115415": "STAT1",
    "ENSG00000157601": "MX1",
    "ENSG00000187608": "ISG15",
    "ENSG00000169245": "CXCL10",
    "ENSG00000232810": "TNF",
    "ENSG00000125538": "IL1B",

    # Entrez Gene IDs mapping to HGNC official symbol
    "3569": "IL6",
    "3456": "IFNB1",
    "7157": "TP53",
    "672": "BRCA1",
    "6772": "STAT1",
    "4599": "MX1",
    "9636": "ISG15",
    "3627": "CXCL10",
    "7124": "TNF",
    "3553": "IL1B"
}

class GeneIDNormalizer:
    """
    Handles translation and normalization of mixed gene identifiers (Ensembl, Entrez,
    official symbols, and aliases) to standardized HGNC Official Gene Symbols.
    """
    def __init__(self, db_path: Optional[Path] = None):
        if db_path is None:
            # Set default database location inside config or data directory
            self.db_path = Path("data/gene_registry.db")
        else:
            self.db_path = Path(db_path)
            
        self._initialize_database()

    def _initialize_database(self):
        """Initializes the local SQLite gene mapping database and seeds it with curated sets."""
        try:
            self.db_path.parent.mkdir(parents=True, exist_ok=True)
            conn = sqlite3.connect(str(self.db_path))
            cursor = conn.cursor()
            
            # Create schema for ID mappings
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS gene_mappings (
                    input_id TEXT PRIMARY KEY,
                    official_symbol TEXT NOT NULL
                )
            """)
            
            # Seed local mapping registry from our curated dictionary to ensure full offline robustness
            for input_id, official in CURATED_ID_MAPS.items():
                cursor.execute(
                    "INSERT OR IGNORE INTO gene_mappings (input_id, official_symbol) VALUES (?, ?)",
                    (input_id.upper(), official.upper())
                )
                
            conn.commit()
            conn.close()
            logger.info("Local Gene ID registry SQLite database initialized successfully.")
        except Exception as e:
            logger.warning(f"Failed to initialize SQLite gene registry database: {e}. Falling back to in-memory mappings.")

    def lookup_local_mapping(self, gene_id: str) -> Optional[str]:
        """Queries the local registry database or internal dictionary for the mapping."""
        clean_id = gene_id.strip().upper()
        
        # 1. Check curated internal dictionary first
        if clean_id in CURATED_ID_MAPS:
            return CURATED_ID_MAPS[clean_id]
            
        # 2. Query SQLite cache
        if self.db_path.exists():
            try:
                conn = sqlite3.connect(str(self.db_path))
                cursor = conn.cursor()
                cursor.execute("SELECT official_symbol FROM gene_mappings WHERE input_id = ?", (clean_id,))
                row = cursor.fetchone()
                conn.close()
                if row:
                    return row[0]
            except Exception as e:
                logger.debug(f"Local mapping SQLite cache query failed: {e}")
                
        return None

    def query_hgnc_api(self, gene_id: str) -> Optional[str]:
        """
        Dynamically queries the remote HGNC REST API as a fallback connection when online.
        """
        clean_id = gene_id.strip()
        headers = {"Accept": "application/json"}
        
        # Detect Ensembl ID
        if clean_id.upper().startswith("ENSG"):
            url = f"https://rest.genenames.org/fetch/ensembl_gene_id/{clean_id}"
        # Detect Entrez ID (numeric only)
        elif clean_id.isdigit():
            url = f"https://rest.genenames.org/fetch/entrez_id/{clean_id}"
        # Default to symbol or synonym alias search
        else:
            url = f"https://rest.genenames.org/search/alias/{clean_id}"
            
        try:
            logger.info(f"Querying HGNC REST API for identifier: {clean_id}")
            response = requests.get(url, headers=headers, timeout=5)
            if response.status_code == 200:
                data = response.json()
                docs = data.get("response", {}).get("docs", [])
                if docs:
                    official_symbol = docs[0].get("symbol")
                    if official_symbol:
                        logger.info(f"HGNC API matched '{clean_id}' to official symbol '{official_symbol}'")
                        
                        # Cache the newly discovered mapping locally
                        self.cache_mapping(clean_id, official_symbol)
                        return official_symbol
        except Exception as e:
            logger.warning(f"HGNC REST API request failed for '{clean_id}': {e}. Using offline mapping boundaries.")
            
        return None

    def cache_mapping(self, input_id: str, official_symbol: str):
        """Saves a newly discovered mapping inside the local SQLite database."""
        try:
            conn = sqlite3.connect(str(self.db_path))
            cursor = conn.cursor()
            cursor.execute(
                "INSERT OR REPLACE INTO gene_mappings (input_id, official_symbol) VALUES (?, ?)",
                (input_id.strip().upper(), official_symbol.strip().upper())
            )
            conn.commit()
            conn.close()
        except Exception as e:
            logger.debug(f"Failed to cache mapping '{input_id} -> {official_symbol}' in SQLite database: {e}")

    def normalize_id(self, gene_id: str, online_fallback: bool = True) -> str:
        """
        Main interface to resolve an unknown gene identifier.
        Checks curated database first, SQLite registry second, and HGNC API last.
        """
        if not gene_id or not str(gene_id).strip():
            return "UNKNOWN"
            
        clean_id = str(gene_id).strip()
        
        # 1. Local checks
        matched = self.lookup_local_mapping(clean_id)
        if matched:
            return matched
            
        # 2. Remote check if configured
        if online_fallback:
            matched = self.query_hgnc_api(clean_id)
            if matched:
                return matched
                
        # If unresolved, return uppercase stripped ID as clean fallback symbol
        return clean_id.upper()

    def normalize_gene_list(self, gene_list: List[str], online_fallback: bool = True) -> List[str]:
        """Normalizes an array of mixed IDs, deduplicating the output list."""
        normalized = []
        seen = set()
        
        for g_id in gene_list:
            norm = self.normalize_id(g_id, online_fallback=online_fallback)
            if norm and norm != "UNKNOWN" and norm not in seen:
                normalized.append(norm)
                seen.add(norm)
                
        logger.info(f"Normalized {len(gene_list)} mixed gene IDs to {len(normalized)} unique official HGNC symbols.")
        return normalized
