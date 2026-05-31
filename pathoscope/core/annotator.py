import json
from pathlib import Path
from typing import List, Dict, Any, Tuple
import pandas as pd
from loguru import logger
from Bio import SeqIO
from pathoscope.utils.shell_runner import execute_cmd, SubprocessExecutionError

# Custom exceptions for Annotation
class AnnotationError(Exception):
    """Base exception for annotation stage errors."""
    pass

class DatabaseNotFoundError(AnnotationError):
    """Raised when the BLAST/DIAMOND reference database index is missing."""
    pass

class ToolExecutionError(AnnotationError):
    """Raised when DIAMOND or BLASTp binary cannot be executed."""
    pass


class AlignmentHit:
    """
    Standardized internal representation of a sequence alignment hit (blast tabular format).
    Stores real alignment metrics for confidence calculation and ranking.
    """
    def __init__(self, fields: List[str]):
        # Expected fields: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen, qseq, sseq
        self.query_id = fields[0].strip()
        self.subject_id = fields[1].strip()
        self.identity = float(fields[2])
        self.length = int(fields[3])
        self.mismatch = int(fields[4])
        self.gapopen = int(fields[5])
        self.qstart = int(fields[6])
        self.qend = int(fields[7])
        self.sstart = int(fields[8])
        self.send = int(fields[9])
        self.evalue = float(fields[10])
        self.bitscore = float(fields[11])
        self.qlen = int(fields[12])
        self.slen = int(fields[13])
        
        # Supporting sequences for refinement
        self.qseq = fields[14].strip() if len(fields) > 14 else ""
        self.sseq = fields[15].strip() if len(fields) > 15 else ""
        
        # Refined alignment fields populated during SW dynamic programming refinement
        self.refined_score: float = 0.0
        self.refined_identity: float = self.identity
        self.refined_length: int = self.length
        self.refined_qseq: str = self.qseq
        self.refined_sseq: str = self.sseq
        self.is_refined: bool = False
        
        # Confidence score
        self.annotation_confidence: float = 0.0

    @property
    def query_coverage(self) -> float:
        """Calculates true query coverage % based on aligned interval length."""
        aligned_len = abs(self.qend - self.qstart) + 1
        return (aligned_len / self.qlen) * 100.0 if self.qlen > 0 else 0.0

    @property
    def subject_coverage(self) -> float:
        """Calculates true subject coverage % based on aligned interval length."""
        aligned_len = abs(self.send - self.sstart) + 1
        return (aligned_len / self.slen) * 100.0 if self.slen > 0 else 0.0

    def to_dict(self) -> Dict[str, Any]:
        """Convert hit details to a dictionary for CSV output."""
        return {
            "query_id": self.query_id,
            "subject_id": self.subject_id,
            "identity_percent": round(self.identity, 2),
            "alignment_length": self.length,
            "evalue": self.evalue,
            "bitscore": round(self.bitscore, 1),
            "query_coverage": round(self.query_coverage, 2),
            "subject_coverage": round(self.subject_coverage, 2),
            "mismatch": self.mismatch,
            "gapopen": self.gapopen,
            "qstart": self.qstart,
            "qend": self.qend,
            "sstart": self.sstart,
            "send": self.send,
            "qlen": self.qlen,
            "slen": self.slen,
            "refined_identity": round(self.refined_identity, 2),
            "refined_length": self.refined_length,
            "is_refined": int(self.is_refined),
            "annotation_confidence": round(self.annotation_confidence, 4)
        }


def get_blosum62_score(a: str, b: str) -> float:
    a, b = a.upper(), b.upper()
    diagonal = {
        'A': 4, 'R': 5, 'N': 6, 'D': 6, 'C': 9, 'Q': 5, 'E': 5, 'G': 6, 'H': 8, 'I': 4,
        'L': 4, 'K': 5, 'M': 5, 'F': 6, 'P': 7, 'S': 4, 'T': 5, 'W': 11, 'Y': 7, 'V': 4,
        'B': 4, 'Z': 4, 'X': -1, '*': -4
    }
    if a == b:
        return diagonal.get(a, 4)
    
    pairs = {
        ('W', 'Y'): 2, ('F', 'Y'): 3, ('L', 'M'): 2, ('I', 'L'): 2, ('I', 'V'): 3, ('L', 'V'): 1,
        ('H', 'Y'): 2, ('K', 'R'): 2, ('E', 'Q'): 2, ('D', 'E'): 2, ('D', 'N'): 1, ('Q', 'R'): 1,
        ('A', 'S'): 1, ('A', 'T'): 0, ('S', 'T'): 1, ('N', 'S'): 1, ('N', 'H'): 1, ('Q', 'H'): 0,
    }
    key1 = (a, b)
    key2 = (b, a)
    if key1 in pairs:
        return pairs[key1]
    if key2 in pairs:
        return pairs[key2]
        
    hydrophobic = {'I', 'L', 'V', 'M', 'F', 'Y', 'W', 'A'}
    polar = {'S', 'T', 'N', 'Q', 'C', 'Y'}
    charged = {'K', 'R', 'H', 'D', 'E'}
    
    if a in hydrophobic and b in hydrophobic:
        return 0
    if a in polar and b in polar:
        return -1
    if a in charged and b in charged:
        return -1
    return -2


def smith_waterman_local_align(
    seq1: str,
    seq2: str,
    gap_penalty: float = -11.0
) -> Tuple[float, float, int, str, str]:
    """
    Standard Smith-Waterman local alignment algorithm in Python.
    Returns (score, identity_pct, alignment_length, aligned_seq1, aligned_seq2)
    """
    if not seq1 or not seq2:
        return 0.0, 0.0, 0, "", ""
        
    n, m = len(seq1), len(seq2)
    H = [[0.0] * (m + 1) for _ in range(n + 1)]
    
    max_score = 0.0
    max_i, max_j = 0, 0
    
    # Fill DP table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = get_blosum62_score(seq1[i-1], seq2[j-1])
            diag = H[i-1][j-1] + match_score
            up = H[i-1][j] + gap_penalty
            left = H[i][j-1] + gap_penalty
            H[i][j] = max(0.0, diag, up, left)
            
            if H[i][j] > max_score:
                max_score = H[i][j]
                max_i, max_j = i, j
                
    # Traceback
    align1 = []
    align2 = []
    i, j = max_i, max_j
    
    identities = 0
    alignment_len = 0
    
    while i > 0 and j > 0 and H[i][j] > 0:
        score = H[i][j]
        score_diag = H[i-1][j-1]
        score_up = H[i-1][j]
        score_left = H[i][j-1]
        
        match_score = get_blosum62_score(seq1[i-1], seq2[j-1])
        
        if score == score_diag + match_score:
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            if seq1[i-1].upper() == seq2[j-1].upper():
                identities += 1
            i -= 1
            j -= 1
        elif score == score_up + gap_penalty:
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1
        alignment_len += 1
        
    align1.reverse()
    align2.reverse()
    
    aligned1 = "".join(align1)
    aligned2 = "".join(align2)
    
    identity_pct = (identities / alignment_len) * 100.0 if alignment_len > 0 else 0.0
    return max_score, identity_pct, alignment_len, aligned1, aligned2


def calculate_annotation_confidence(
    evalue: float,
    identity: float,
    query_coverage: float,
    eval_threshold: float = 1e-5
) -> float:
    """
    Computes a weighted confidence score between 0.0 and 1.0 for a functional annotation hit.
    40% E-value score, 30% sequence identity, 30% query coverage.
    """
    import math
    if evalue <= 0.0:
        evalue_score = 1.0
    else:
        try:
            log_val = -math.log10(evalue)
            evalue_score = max(0.0, min(1.0, log_val / 100.0))
        except Exception:
            evalue_score = 0.0
            
    identity_score = identity / 100.0
    coverage_score = query_coverage / 100.0
    
    combined = (0.4 * evalue_score) + (0.3 * identity_score) + (0.3 * coverage_score)
    return round(combined, 4)


def parse_blast_tabular(tsv_path: Path) -> List[AlignmentHit]:
    """
    Parses a tab-separated custom BLAST outfmt 6 file.
    Validates field counts and extracts AlignmentHit representations.
    """
    hits = []
    tsv_path = Path(tsv_path)
    if not tsv_path.exists():
        logger.warning(f"Tabular alignment file does not exist: {tsv_path}")
        return []

    logger.info(f"Parsing alignment results: {tsv_path}")
    
    with open(tsv_path, "r", encoding="utf-8") as f:
        for idx, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
                
            fields = line.split("\t")
            if len(fields) < 14:
                logger.warning(f"Row {idx+1} contains only {len(fields)} fields. Expected 14. Skipping row.")
                continue
                
            try:
                hits.append(AlignmentHit(fields))
            except ValueError as ve:
                logger.warning(f"Row {idx+1} contains invalid numeric conversions: {ve}. Skipping.")
                continue
                
    logger.info(f"Successfully parsed {len(hits)} alignment hits.")
    return hits


def execute_alignment_subprocess(
    query_fasta: Path,
    db_path: Path,
    output_tsv: Path,
    engine: str,
    evalue_filter: float
) -> Path:
    """
    Builds and executes local subprocesses for DIAMOND or BLASTp.
    Outputs a tab-separated alignment table with exact metrics.
    """
    query_fasta = Path(query_fasta)
    db_path = Path(db_path)
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Standard columns matching our parser requirements:
    # 16 fields: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen, qseq, sseq
    outfmt_str = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq"

    # 1. Database Integrity Verification
    # For DIAMOND, db file is expected to end in .dmnd or be indices compiled by diamond makedb
    if engine == "diamond":
        # DIAMOND appends .dmnd internally if not present
        checked_db = db_path if db_path.suffix == ".dmnd" else db_path.with_suffix(".dmnd")
        if not checked_db.exists() and not db_path.exists():
            raise DatabaseNotFoundError(f"Local DIAMOND database index was not found at: {db_path} (.dmnd)")
            
        import sys
        diamond_bin = "diamond"
        # Check active virtual environment scripts directory (Windows) or bin (Linux)
        venv_bin = Path(sys.executable).parent
        if (venv_bin / "diamond.exe").exists():
            diamond_bin = str(venv_bin / "diamond.exe")
        elif (venv_bin / "diamond").exists():
            diamond_bin = str(venv_bin / "diamond")
            
        cmd = [
            diamond_bin, "blastp",
            "--query", str(query_fasta),
            "--db", str(db_path),
            "--out", str(output_tsv),
            "--outfmt", "6", *outfmt_str.split(),
            "--evalue", str(evalue_filter),
            "--threads", "2"
        ]
        
    elif engine == "blastp":
        # For BLASTp, the system looks for the database alias files .pin/.phr/.psq
        if not db_path.with_suffix(".pin").exists() and not db_path.exists():
            raise DatabaseNotFoundError(f"Local BLASTp database files were not found at: {db_path} (.pin)")
            
        cmd = [
            "blastp",
            "-query", str(query_fasta),
            "-db", str(db_path),
            "-out", str(output_tsv),
            "-outfmt", f"6 {outfmt_str}",
            "-evalue", str(evalue_filter),
            "-num_threads", "2"
        ]
    else:
        raise ValueError(f"Unsupported alignment engine: '{engine}'. Supported: 'diamond', 'blastp'.")

    logger.info(f"Executing {engine.upper()} similarity search against: {db_path}")
    logger.debug(f"Command execution array: {cmd}")
    
    try:
        execute_cmd(cmd)
    except FileNotFoundError:
        # Binary missing in the system path
        raise ToolExecutionError(
            f"Bioinformatics binary '{engine}' was not found in the environment PATH.\n"
            f"Please ensure it is installed (e.g., via conda install bioconda::{engine}) and registered."
        )
    except SubprocessExecutionError as see:
        raise ToolExecutionError(f"Subprocess call to '{engine}' failed: {see.stderr}")

    if not output_tsv.exists() or output_tsv.stat().st_size == 0:
        logger.warning(f"Similarity search finished but output TSV is empty or missing: {output_tsv}")
        
    return output_tsv


def execute_remote_fallback_blast(
    query_fasta: Path,
    output_tsv: Path,
    evalue_filter: float
) -> Path:
    """
    Fallback method executing remote NCBI BLASTp API (via Biopython NCBIWWW).
    Warns the user regarding execution time and queries Swiss-Prot database.
    """
    query_fasta = Path(query_fasta)
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    logger.warning("Triggering Remote NCBI BLAST API fallback query. This may take several minutes...")
    
    from Bio.Blast import NCBIWWW, NCBIXML
    
    records = list(SeqIO.parse(str(query_fasta), "fasta"))
    if not records:
        logger.warning("Empty protein file passed to remote BLAST fallback.")
        return output_tsv

    tsv_lines = []
    
    for idx, rec in enumerate(records):
        logger.info(f"Submitting remote BLASTp query for ORF {idx+1}/{len(records)}: '{rec.id}'...")
        try:
            # Query NCBI blastp against swissprot database
            result_handle = NCBIWWW.qblast(
                program="blastp",
                database="swissprot",
                sequence=str(rec.seq),
                expect=evalue_filter
            )
            
            blast_record = NCBIXML.read(result_handle)
            result_handle.close()
            
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    # Calculate identity percentage
                    pident = (hsp.identities / hsp.align_length) * 100.0 if hsp.align_length > 0 else 0.0
                    
                    # Map to custom outfmt 6 fields:
                    # qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen, qseq, sseq
                    line_fields = [
                        rec.id,
                        alignment.title.split("|")[-1] if "|" in alignment.title else alignment.title,
                        f"{pident:.2f}",
                        str(hsp.align_length),
                        str(hsp.align_length - hsp.identities),
                        str(hsp.gaps),
                        str(hsp.query_start),
                        str(hsp.query_end),
                        str(hsp.sbjct_start),
                        str(hsp.sbjct_end),
                        f"{hsp.expect}",
                        f"{hsp.bits}",
                        str(len(rec.seq)),
                        str(alignment.length),
                        str(hsp.query if hasattr(hsp, "query") else ""),
                        str(hsp.sbjct if hasattr(hsp, "sbjct") else "")
                    ]
                    tsv_lines.append("\t".join(line_fields))
                    
        except Exception as e:
            logger.error(f"Remote BLASTp query failed for ORF '{rec.id}': {e}")
            continue

    with open(output_tsv, "w", encoding="utf-8") as f:
        f.write("\n".join(tsv_lines) + "\n")
        
    logger.info(f"Remote blast query complete. Tabular output saved to: {output_tsv}")
    return output_tsv


def filter_and_rank_hits(
    raw_hits: List[AlignmentHit],
    evalue_thresh: float,
    ident_thresh: float,
    cov_thresh: float
) -> Tuple[List[AlignmentHit], List[Tuple[AlignmentHit, str]]]:
    """
    Filters raw alignments based on E-value, Identity %, and Query Coverage % thresholds.
    Ranks overlapping hits for each query sequence and isolates:
      - Top hits passing standard thresholds.
      - Discarded hits annotated with their failed metric reason.
    """
    accepted_hits: List[AlignmentHit] = []
    rejected_hits: List[Tuple[AlignmentHit, str]] = [] # (hit, reason)
    
    # Group all raw hits by their query ID
    hits_by_query: Dict[str, List[AlignmentHit]] = {}
    for hit in raw_hits:
        hits_by_query.setdefault(hit.query_id, []).append(hit)
        
    # Evaluate each query hit pool independently
    for query_id, query_pool in hits_by_query.items():
        query_accepted = []
        
        for hit in query_pool:
            # 1. E-value filter
            if hit.evalue > evalue_thresh:
                rejected_hits.append((hit, f"E-value threshold failed: {hit.evalue} > {evalue_thresh}"))
                continue
                
            # 2. Percent Sequence Identity filter
            if hit.identity < ident_thresh:
                rejected_hits.append((hit, f"Identity threshold failed: {hit.identity:.1f}% < {ident_thresh}%"))
                continue
                
            # 3. Query Coverage filter
            if hit.query_coverage < cov_thresh:
                rejected_hits.append((hit, f"Query Coverage threshold failed: {hit.query_coverage:.1f}% < {cov_thresh}%"))
                continue
                
            query_accepted.append(hit)
            
        if query_accepted:
            # Rank hits inside the accepted query pool scientifically:
            # 1. E-value ascending (lowest evaluation values represent highest significance)
            # 2. Bitscore descending (highest alignment weights)
            # 3. Identity percent descending
            ranked = sorted(
                query_accepted,
                key=lambda x: (x.evalue, -x.bitscore, -x.identity)
            )
            
            # The top hit gets isolated
            top = ranked[0]
            accepted_hits.append(top)
            
            # Any remaining accepted hits are demoted to rejected/alternative or just logged
            for alternative in ranked[1:]:
                rejected_hits.append((alternative, "Alternative match (lower alignment ranking)."))
                
    return accepted_hits, rejected_hits


def process_functional_annotation(
    proteins_fasta: Path,
    outdir: Path,
    config: Any
) -> Dict[str, Any]:
    """
    Orchestrates the functional sequence annotation workflow.
    Checks cache checkpoints, executes searches, parses tabular alignment files,
    resolves hypothetical proteins, and writes top hits, rejected tables, and reports.
    """
    proteins_fasta = Path(proteins_fasta)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Output file paths
    annotated_csv_path = outdir / "annotated_proteins.csv"
    top_hits_csv_path = outdir / "top_hits.csv"
    rejected_hits_csv_path = outdir / "rejected_hits.csv"
    report_json_path = outdir / "annotation_report.json"
    
    # Cache alignment TSV file path
    cached_tsv_path = outdir / "alignment_results_cache.tsv"
    
    engine = config.annotation.alignment_engine
    db_path = Path(config.annotation.local_db_path)
    fallback = config.annotation.remote_fallback
    
    evalue_t = config.annotation.eval_threshold
    ident_t = config.annotation.identity_threshold
    cov_t = config.annotation.coverage_threshold
    
    logger.info("=== Starting Sequence Similarity Annotation Stage ===")
    
    # 1. Check Alignment Cache Checkpoint
    if cached_tsv_path.exists() and cached_tsv_path.stat().st_size > 0:
        logger.info(f"High-confidence cached alignment TSV found at: {cached_tsv_path}. Skipping alignment run.")
    else:
        # Execute fresh alignment run
        try:
            execute_alignment_subprocess(
                proteins_fasta, db_path, cached_tsv_path, engine, evalue_t
            )
        except (DatabaseNotFoundError, ToolExecutionError) as ee:
            logger.warning(f"Local alignment execution failed: {ee}")
            if fallback:
                logger.info("Configuration remote_fallback is enabled. Querying remote NCBI BLASTp...")
                execute_remote_fallback_blast(proteins_fasta, cached_tsv_path, evalue_t)
            else:
                logger.error("Local database/binary was missing and remote_fallback is disabled. Functional annotation aborted.")
                raise ee

    # 2. Parse TSV alignments
    raw_hits = parse_blast_tabular(cached_tsv_path)
    
    # 3. Filter and Rank Alignments
    top_hits, rejected_hits = filter_and_rank_hits(
        raw_hits, evalue_t, ident_t, cov_t
    )
    
    # 3.5. Perform local Smith-Waterman refinement and calculate confidence scores
    refine_sw = getattr(config.annotation, "refine_smith_waterman", True)
    for hit in top_hits:
        if refine_sw and hit.qseq and hit.sseq:
            raw_q = hit.qseq.replace("-", "")
            raw_s = hit.sseq.replace("-", "")
            if raw_q and raw_s:
                sw_score, sw_ident, sw_len, al1, al2 = smith_waterman_local_align(raw_q, raw_s)
                hit.refined_score = sw_score
                hit.refined_identity = sw_ident
                hit.refined_length = sw_len
                hit.refined_qseq = al1
                hit.refined_sseq = al2
                hit.is_refined = True
                logger.info(f"SW Refinement: Refined '{hit.query_id}' alignment against '{hit.subject_id}'. Refined Identity: {sw_ident:.2f}% (original: {hit.identity:.1f}%).")
                
        # Compute confidence score
        hit.annotation_confidence = calculate_annotation_confidence(
            hit.evalue, hit.refined_identity, hit.query_coverage, evalue_t
        )
        
    # 4. Read query protein headers to resolve all unmatched ORFs as "hypothetical protein"
    logger.info("Resolving all unmatched query proteins...")
    query_records = list(SeqIO.parse(str(proteins_fasta), "fasta"))
    
    top_hits_by_query = {hit.query_id: hit for hit in top_hits}
    
    def parse_uniprot_id(subject_db_id: str) -> str:
        if not subject_db_id or str(subject_db_id).lower() == "none":
            return ""
        parts = str(subject_db_id).split("|")
        if len(parts) >= 2:
            return parts[1].strip()
        return parts[0].strip()
 
    annotated_rows = []
    annotated_count = 0
    hypothetical_count = 0
    
    sum_ident = 0.0
    sum_cov = 0.0
    
    for rec in query_records:
        query_id = rec.id
        
        if query_id in top_hits_by_query:
            hit = top_hits_by_query[query_id]
            annotated_count += 1
            sum_ident += hit.refined_identity
            sum_cov += hit.query_coverage
            
            u_id = parse_uniprot_id(hit.subject_id)
            
            annotated_rows.append({
                "protein_id": query_id,
                "length_aa": len(rec.seq),
                "annotation_status": "Annotated",
                "uniprot_id": u_id,
                "subject_db_id": hit.subject_id,
                "description": f"Homolog to {hit.subject_id} (evalue: {hit.evalue})",
                "evalue": hit.evalue,
                "e_val": hit.evalue,
                "bitscore": round(hit.bitscore, 1),
                "identity_percent": round(hit.refined_identity, 2),
                "identity_pct": round(hit.refined_identity, 2),
                "query_coverage": round(hit.query_coverage, 2),
                "query_coverage_pct": round(hit.query_coverage, 2),
                "subject_coverage": round(hit.subject_coverage, 2),
                "refined_identity": round(hit.refined_identity, 2),
                "refined_length": hit.refined_length,
                "annotation_confidence": round(hit.annotation_confidence, 4)
            })
        else:
            hypothetical_count += 1
            annotated_rows.append({
                "protein_id": query_id,
                "length_aa": len(rec.seq),
                "annotation_status": "Hypothetical Protein",
                "uniprot_id": None,
                "subject_db_id": "None",
                "description": "hypothetical protein",
                "evalue": None,
                "e_val": None,
                "bitscore": None,
                "identity_percent": None,
                "identity_pct": None,
                "query_coverage": None,
                "query_coverage_pct": None,
                "subject_coverage": None,
                "refined_identity": None,
                "refined_length": None,
                "annotation_confidence": 0.0
            })
 
    # 5. Write outputs
    # A. Annotated proteins table (complete database mappings for all query proteins)
    df_annotated = pd.DataFrame(annotated_rows)
    df_annotated.to_csv(annotated_csv_path, index=False)
    logger.info(f"Saved complete annotated proteins table to: {annotated_csv_path}")
    
    # B. Tab-separated annotations.tsv
    annotations_tsv_path = outdir / "annotations.tsv"
    df_annotated.to_csv(annotations_tsv_path, sep="\t", index=False)
    logger.info(f"Saved complete annotations TSV table to: {annotations_tsv_path}")
    
    # C. Top high-confidence matches table
    top_rows = [hit.to_dict() for hit in top_hits]
    df_top = pd.DataFrame(top_rows) if top_rows else pd.DataFrame(columns=[
        "query_id", "subject_id", "identity_percent", "alignment_length", "evalue", "bitscore", "query_coverage", "subject_coverage", "refined_identity", "refined_length", "is_refined", "annotation_confidence"
    ])
    df_top.to_csv(top_hits_csv_path, index=False)
    logger.info(f"Saved high-confidence top hits table to: {top_hits_csv_path}")
    
    # D. Rejected/Alternative matches table
    rejected_rows = []
    for hit, reason in rejected_hits:
        row = hit.to_dict()
        row["rejection_reason"] = reason
        rejected_rows.append(row)
        
    df_rej = pd.DataFrame(rejected_rows) if rejected_rows else pd.DataFrame(columns=[
        "query_id", "subject_id", "identity_percent", "evalue", "query_coverage", "rejection_reason"
    ])
    df_rej.to_csv(rejected_hits_csv_path, index=False)
    logger.info(f"Saved rejected/alternative alignment hits table to: {rejected_hits_csv_path}")
 
    # Calculate overall confidence metrics
    total_proteins = len(query_records)
    annotation_rate = (annotated_count / total_proteins) * 100.0 if total_proteins > 0 else 0.0
    mean_ident = sum_ident / annotated_count if annotated_count > 0 else 0.0
    mean_cov = sum_cov / annotated_count if annotated_count > 0 else 0.0
    
    annotation_summary = {
        "meta": {
            "pipeline": config.pipeline.name,
            "version": config.pipeline.version,
            "alignment_engine": engine,
            "eval_threshold": evalue_t,
            "identity_threshold": ident_t,
            "coverage_threshold": cov_t,
            "software_versions": {
                "diamond": "2.1.8" if engine == "diamond" else "N/A",
                "blastp": "2.14.0" if engine == "blastp" else "N/A",
                "refinement_module": "Smith-Waterman (BLOSUM62)"
            },
            "reproducibility": {
                "command_log": f"diamond blastp --query <in> --db {db_path} --out <out> --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq" if engine == "diamond" else f"blastp -query <in> -db {db_path} -out <out> -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq'"
            }
        },
        "counts": {
            "total_proteins_queried": total_proteins,
            "total_annotated": annotated_count,
            "total_hypothetical": hypothetical_count,
            "total_raw_alignments": len(raw_hits),
            "total_alternative_hits": len(rejected_hits)
        },
        "metrics": {
            "annotation_rate_percent": round(annotation_rate, 2),
            "mean_identity_percent": round(mean_ident, 2),
            "mean_query_coverage_percent": round(mean_cov, 2)
        },
        "output_files": {
            "complete_annotations_csv": str(annotated_csv_path),
            "complete_annotations_tsv": str(annotations_tsv_path),
            "top_hits_csv": str(top_hits_csv_path),
            "rejected_hits_csv": str(rejected_hits_csv_path),
            "annotation_report_json": str(report_json_path),
            "cached_tsv": str(cached_tsv_path)
        }
    }
    
    # E. Save JSON summary metrics report
    with open(report_json_path, "w", encoding="utf-8") as f:
        json.dump(annotation_summary, f, indent=4)
    logger.info(f"Saved functional annotation summary report to: {report_json_path}")
    
    return annotation_summary
