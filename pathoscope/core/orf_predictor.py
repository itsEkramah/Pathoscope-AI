import json
from pathlib import Path
from typing import List, Dict, Any, Tuple
import pandas as pd
from loguru import logger
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from pathoscope.core.preprocessor import SequenceRecord

# Custom exceptions for ORF prediction
class ORFError(Exception):
    """Base exception for ORF prediction errors."""
    pass

class InvalidCoordinateError(ORFError):
    """Raised when coordinates exceed sequence length boundaries."""
    pass


class ORFRecord:
    """
    Standardized internal representation of a predicted Open Reading Frame (ORF).
    Uses 1-based, inclusive coordinates relative to the forward strand,
    consistent with the GFF3 standard.
    """
    def __init__(
        self,
        orf_id: str,
        sequence_id: str,
        start: int,       # 1-based, inclusive start coordinate
        end: int,         # 1-based, inclusive end coordinate
        strand: str,      # "+" or "-"
        frame: int,       # 1, 2, or 3 (or negative)
        nucleotide_seq: str,
        protein_seq: str,
        start_codon: str,
        confidence_score: float = 0.0
    ):
        self.id = orf_id
        self.seq_id = sequence_id
        self.start = start
        self.end = end
        self.strand = strand
        self.frame = frame
        self.nuc_seq = nucleotide_seq.upper()
        self.prot_seq = protein_seq.upper()
        self.start_codon = start_codon
        self.confidence_score = confidence_score
        
        # Metadata populated during downstream steps
        self.overlap_flag: str = "None"  # "None", "Overlap", "Nested"
        self.overlap_with: List[str] = []

    @property
    def length_bp(self) -> int:
        return len(self.nuc_seq)

    @property
    def length_aa(self) -> int:
        return len(self.prot_seq)

    def to_dict(self) -> Dict[str, Any]:
        """Convert record details to dictionary for CSV/JSON serialization."""
        return {
            "orf_id": self.id,
            "sequence_id": self.seq_id,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "frame": self.frame,
            "length_bp": self.length_bp,
            "length_aa": self.length_aa,
            "start_codon": self.start_codon,
            "confidence_score": self.confidence_score,
            "overlap_flag": self.overlap_flag,
            "overlap_with": ",".join(self.overlap_with) if self.overlap_with else ""
        }

    def to_gff_line(self) -> str:
        """Generates a standard GFF3 tab-separated format string."""
        attributes = f"ID={self.id};Name={self.id};frame={self.frame};start_codon={self.start_codon};overlap={self.overlap_flag};confidence_score={self.confidence_score}"
        if self.overlap_with:
            attributes += f";overlap_with={','.join(self.overlap_with)}"
            
        # GFF3 fields: seqid, source, type, start, end, score, strand, phase, attributes
        return "\t".join([
            self.seq_id,
            "PathoScope_ORF",
            "CDS",
            str(self.start),
            str(self.end),
            str(self.confidence_score),  # score
            self.strand,
            "0",  # phase
            attributes
        ])


def reverse_complement_coordinates(
    rc_start_idx: int,
    rc_end_idx: int,
    seq_length: int
) -> Tuple[int, int]:
    """
    Translates 0-based coordinates from the reverse complement strand
    back to 1-based, inclusive coordinates on the forward strand.
    """
    # 0-based forward coordinates:
    fwd_stop_idx = seq_length - 1 - rc_start_idx
    fwd_start_idx = seq_length - 1 - rc_end_idx
    
    # 1-based, inclusive coordinates:
    start_fwd = fwd_start_idx + 1
    end_fwd = fwd_stop_idx + 1
    
    if start_fwd > end_fwd:
        raise InvalidCoordinateError("Reverse complement coordinate mapping returned invalid start/end positions.")
        
    return start_fwd, end_fwd


def calculate_confidence(
    length_aa: int,
    start_codon: str,
    nuc_seq: str,
    bg_gc_content: float = 0.5
) -> float:
    """
    Calculates a confidence score between 0.0 and 1.0 for a predicted ORF.
    Based on:
    1. Start codon weight (ATG = 1.0, GTG = 0.7, TTG = 0.5, others = 0.1)
    2. Length scaling (sigmoidal probability; longer ORFs are highly confident)
    3. GC content bias relative to standard background (GC content within typical viral ranges)
    """
    import math
    # Start codon weight
    codon_weights = {"ATG": 1.0, "GTG": 0.7, "TTG": 0.5}
    start_weight = codon_weights.get(start_codon.upper(), 0.1)
    
    # Length scoring (sigmoidal)
    # Standard virus gene size is around 100-300 aa.
    # Let's set a midpoint at 80 aa, with transition speed 25.
    length_score = 1.0 / (1.0 + math.exp(-(length_aa - 80) / 25.0))
    
    # GC content scoring
    # Calculate GC of the ORF
    gc_count = nuc_seq.count('G') + nuc_seq.count('C')
    orf_gc = gc_count / len(nuc_seq) if len(nuc_seq) > 0 else 0.5
    
    # Deviation from background GC. High deviation in either direction is slightly penalized.
    gc_deviation = abs(orf_gc - bg_gc_content)
    gc_score = math.exp(-2.0 * (gc_deviation ** 2)) # gaussian penalty
    
    # Combine metrics with weights:
    # 50% length, 30% start codon, 20% GC score
    combined = (0.5 * length_score) + (0.3 * start_weight) + (0.2 * gc_score)
    return round(combined, 4)


def filter_nested_orfs(
    orfs: List[ORFRecord],
    resolve_nested: bool
) -> List[ORFRecord]:
    """
    Identifies ORFs that are completely nested within other larger ORFs on the same strand.
    If resolve_nested is True, discards them. Otherwise, flags them as 'Nested' and lists the container ORF.
    """
    if not orfs:
        return []
        
    # Sort by length descending, so we process larger container ORFs first
    sorted_orfs = sorted(orfs, key=lambda x: -x.length_bp)
    nested_ids = set()
    container_map = {} # maps nested_id -> container_id
    
    n = len(sorted_orfs)
    for i in range(n):
        parent = sorted_orfs[i]
        if parent.id in nested_ids:
            continue
            
        for j in range(i + 1, n):
            child = sorted_orfs[j]
            if child.id in nested_ids:
                continue
                
            # Check if child is fully contained within parent on the same strand
            if child.strand == parent.strand and child.start >= parent.start and child.end <= parent.end:
                nested_ids.add(child.id)
                container_map[child.id] = parent.id
                logger.info(f"Nested Filter: Detected ORF '{child.id}' fully nested inside '{parent.id}' on strand '{child.strand}'.")
                
    # Now update flags and rebuild the list
    filtered: List[ORFRecord] = []
    for orf in orfs:
        if orf.id in nested_ids:
            if resolve_nested:
                # Discard this ORF
                logger.info(f"Nested Filter: Discarding nested ORF '{orf.id}' nested in '{container_map[orf.id]}'.")
                continue
            else:
                # Keep but flag
                orf.overlap_flag = "Nested"
                if container_map[orf.id] not in orf.overlap_with:
                    orf.overlap_with.append(container_map[orf.id])
        filtered.append(orf)
        
    return filtered


LITERATURE_VIRAL_GENOMES = {
    "Bacteriophage MS2 (NC_001417)": {
        "length_range": (3500, 3600),
        "proteins": [
            {"name": "Maturation protein (A)", "length_aa": 393, "strand": "+", "start": 130, "end": 1309},
            {"name": "Coat protein", "length_aa": 130, "strand": "+", "start": 1335, "end": 1727},
            {"name": "Lysis protein", "length_aa": 75, "strand": "+", "start": 1678, "end": 1905},
            {"name": "Replicase (RNA polymerase)", "length_aa": 545, "strand": "+", "start": 1761, "end": 3398}
        ]
    },
    "Enterobacteria phage phiX174 (NC_001422)": {
        "length_range": (5300, 5400),
        "proteins": [
            {"name": "Protein A (DNA replication)", "length_aa": 513, "strand": "+", "start": 3981, "end": 5520},
            {"name": "Protein B", "length_aa": 120, "strand": "+", "start": 5075, "end": 5437},
            {"name": "Protein C", "length_aa": 86, "strand": "+", "start": 5434, "end": 5694},
            {"name": "Protein D", "length_aa": 152, "strand": "+", "start": 390, "end": 848},
            {"name": "Protein E (Lysis)", "length_aa": 91, "strand": "+", "start": 568, "end": 843},
            {"name": "Protein F (Capsid)", "length_aa": 427, "strand": "+", "start": 1001, "end": 2284},
            {"name": "Protein G (Spike)", "length_aa": 175, "strand": "+", "start": 2395, "end": 2922},
            {"name": "Protein H (Pilot)", "length_aa": 328, "strand": "+", "start": 2931, "end": 3917},
            {"name": "Protein J", "length_aa": 37, "strand": "+", "start": 848, "end": 961},
            {"name": "Protein K", "length_aa": 56, "strand": "+", "start": 51, "end": 221}
        ]
    }
}


def compare_predicted_to_literature(
    predicted_orfs: List[ORFRecord],
    seq_len: int,
    seq_id: str
) -> Dict[str, Any]:
    """
    Compares predicted ORFs against standard catalog templates of known viral genomes
    like Bacteriophage MS2 and phiX174, matching by length and approximate positions.
    """
    matched_genome = None
    # Match genome based on length range
    for genome_name, data in LITERATURE_VIRAL_GENOMES.items():
        min_l, max_l = data["length_range"]
        if min_l <= seq_len <= max_l:
            matched_genome = genome_name
            break
            
    if not matched_genome:
        # Fallback: check if the sequence ID contains keywords
        for genome_name, data in LITERATURE_VIRAL_GENOMES.items():
            if "ms2" in seq_id.lower() and "ms2" in genome_name.lower():
                matched_genome = genome_name
                break
            elif "phix" in seq_id.lower() and "phix" in genome_name.lower():
                matched_genome = genome_name
                break
                
    comparison_results = {
        "sequence_id": seq_id,
        "sequence_length": seq_len,
        "matched_reference_genome": matched_genome or "None (No closely matching reference virus in catalog)",
        "matches": []
    }
    
    if not matched_genome:
        return comparison_results
        
    ref_proteins = LITERATURE_VIRAL_GENOMES[matched_genome]["proteins"]
    
    # Try to match predicted ORFs
    for ref in ref_proteins:
        best_match = None
        best_score = 0.0 # combining length and coordinate overlap
        
        for orf in predicted_orfs:
            if orf.strand != ref["strand"]:
                continue
                
            # Length similarity
            len_diff = abs(orf.length_aa - ref["length_aa"])
            len_sim = max(0.0, 1.0 - (len_diff / ref["length_aa"]))
            
            # Coordinate distance (approximate start distance)
            start_diff = abs(orf.start - ref["start"])
            # Scale distance: within 100 bp is excellent, after 500 bp drops to 0
            coord_sim = max(0.0, 1.0 - (start_diff / 500.0))
            
            # Combined score (weighted: 60% length similarity, 40% coordinate similarity)
            score = (0.6 * len_sim) + (0.4 * coord_sim)
            
            # We want a high threshold for a true match, e.g., > 0.7
            if score > 0.7 and score > best_score:
                best_match = orf
                best_score = score
                
        if best_match:
            comparison_results["matches"].append({
                "reference_protein_name": ref["name"],
                "reference_length_aa": ref["length_aa"],
                "reference_start": ref["start"],
                "reference_end": ref["end"],
                "reference_strand": ref["strand"],
                "predicted_orf_id": best_match.id,
                "predicted_length_aa": best_match.length_aa,
                "predicted_start": best_match.start,
                "predicted_end": best_match.end,
                "predicted_strand": best_match.strand,
                "match_confidence": round(best_score, 4)
            })
        else:
            comparison_results["matches"].append({
                "reference_protein_name": ref["name"],
                "reference_length_aa": ref["length_aa"],
                "reference_start": ref["start"],
                "reference_end": ref["end"],
                "reference_strand": ref["strand"],
                "predicted_orf_id": "NOT_DETECTED",
                "predicted_length_aa": None,
                "predicted_start": None,
                "predicted_end": None,
                "predicted_strand": None,
                "match_confidence": 0.0
            })
            
    return comparison_results


def scan_single_strand(
    sequence_id: str,
    seq_str: str,
    strand: str,
    start_codons: List[str],
    stop_codons: List[str],
    min_len_aa: int,
    translation_table: int
) -> List[ORFRecord]:
    """
    Scans a single nucleotide sequence (can be forward or reverse complement)
    in all 3 reading frames. Identifies all possible non-nested ORFs.
    """
    seq_len = len(seq_str)
    raw_orfs: List[ORFRecord] = []
    
    # Pre-compile codons for exact matching speed
    starts = set(c.upper() for c in start_codons)
    stops = set(c.upper() for c in stop_codons)
    
    # Scan three reading frames (0, 1, 2)
    for frame_idx in range(3):
        frame_name = (frame_idx + 1) if strand == "+" else -(frame_idx + 1)
        
        # Track all active start codon indices in this frame
        active_starts = []
        
        # Step through codons in this reading frame
        for idx in range(frame_idx, seq_len - 2, 3):
            codon = seq_str[idx:idx+3]
            
            if codon in starts:
                active_starts.append(idx)
                
            elif codon in stops:
                if active_starts:
                    # Biological logic: Stop codon reached.
                    # Prevent Nested ORF Inflation:
                    # In standard scanning, if multiple start codons lead to the same stop codon,
                    # the ribosome will translate from the outermost (first) start codon.
                    # Retaining all starts would create multiple nested ORFs mapping to the same locus.
                    # Therefore, we keep only the OUTMOST (first) start codon to get the longest ORF!
                    longest_start = active_starts[0]
                    active_starts.clear()  # Clear all active starts since this stop codon consumes them
                    
                    orf_len_bp = (idx + 3) - longest_start
                    orf_len_aa = orf_len_bp // 3 - 1  # excluding stop codon
                    
                    if orf_len_aa >= min_len_aa:
                        orf_nuc = seq_str[longest_start : idx+3]
                        
                        # Handle translation using Biopython
                        biopy_seq = Seq(orf_nuc)
                        try:
                            # Translate excluding the stop codon for clean presentation
                            orf_prot = str(biopy_seq.translate(table=translation_table, to_stop=True))
                        except Exception as te:
                            logger.warning(f"Biopython translation failed for ORF at frame {frame_name}: {te}. Falling back to default table.")
                            orf_prot = str(biopy_seq.translate(to_stop=True))
                            
                        # Resolve coordinates
                        if strand == "+":
                            start_fwd = longest_start + 1
                            end_fwd = idx + 3
                        else:
                            # For reverse strand, coordinate conversion is relative to original length
                            start_fwd, end_fwd = reverse_complement_coordinates(longest_start, idx + 2, seq_len)
                            
                        # Generate temporary identifier
                        temp_id = f"ORF_{sequence_id}_{strand}{abs(frame_name)}_{start_fwd}_{end_fwd}"
                        
                        raw_orfs.append(ORFRecord(
                            orf_id=temp_id,
                            sequence_id=sequence_id,
                            start=start_fwd,
                            end=end_fwd,
                            strand=strand,
                            frame=frame_name,
                            nucleotide_seq=orf_nuc,
                            protein_seq=orf_prot,
                            start_codon=orf_nuc[:3]
                        ))
                        
    return raw_orfs


def resolve_overlaps(
    orfs: List[ORFRecord],
    overlap_threshold: int,
    policy: str
) -> List[ORFRecord]:
    """
    Interval sweep-line algorithm to identify and resolve overlaps scientifically.
    Supports keeping all ORFs and flagging overlaps (ideal for compact viral genomes)
    or keeping longest only to eliminate overlapping false positives.
    """
    if not orfs:
        return []
        
    # Step 1: Sort ORFs by start coordinate, then by length descending
    sorted_orfs = sorted(orfs, key=lambda x: (x.start, -(x.end - x.start)))
    
    resolved: List[ORFRecord] = []
    discarded_ids = set()
    
    # Step 2: Conduct overlap sweep
    n = len(sorted_orfs)
    for i in range(n):
        if sorted_orfs[i].id in discarded_ids:
            continue
            
        current = sorted_orfs[i]
        overlaps_detected: List[ORFRecord] = []
        
        # Check against all downstream ORFs that start before the current one ends
        for j in range(i + 1, n):
            candidate = sorted_orfs[j]
            if candidate.id in discarded_ids:
                continue
                
            if candidate.start >= current.end:
                # No more possible overlaps due to sorting
                break
                
            # Calculate absolute overlap length in bp
            overlap_len = min(current.end, candidate.end) - max(current.start, candidate.start)
            
            # Check if overlap exceeds the configured threshold (bp)
            if overlap_len >= overlap_threshold:
                overlaps_detected.append(candidate)
                
        if overlaps_detected:
            # Overlaps exist: apply configured scientific policy
            if policy == "longest_only":
                # Find the longest ORF in the overlapping cluster
                cluster = [current] + overlaps_detected
                longest_orf = max(cluster, key=lambda x: x.length_bp)
                
                # Keep only the longest, discard all others in the cluster
                for item in cluster:
                    if item.id != longest_orf.id:
                        discarded_ids.add(item.id)
                        logger.warning(f"Overlap Resolution: Discarding overlapping ORF '{item.id}' (length {item.length_bp} bp) in favor of longer '{longest_orf.id}' (length {longest_orf.length_bp} bp).")
                        
                resolved.append(longest_orf)
                
            elif policy == "keep_all_flag":
                # Keep all ORFs but annotate them with flags documenting the overlaps
                current.overlap_flag = "Overlap"
                for item in overlaps_detected:
                    item.overlap_flag = "Overlap"
                    current.overlap_with.append(item.id)
                    item.overlap_with.append(current.id)
                    logger.info(f"Overlap Log: Kept overlapping viral ORFs '{current.id}' & '{item.id}' (Overlap: {min(current.end, item.end) - max(current.start, item.start)} bp).")
                
                resolved.append(current)
        else:
            resolved.append(current)
            
    # Filter out any records that were marked as discarded downstream in clustering
    final_orfs = [o for o in resolved if o.id not in discarded_ids]
    return sorted(final_orfs, key=lambda x: x.start)


def predict_orfs_in_sequence(
    record: SequenceRecord,
    config: Any
) -> List[ORFRecord]:
    """
    Coordinates true 6-frame scanning and overlap resolution for a single sequence.
    """
    min_len_aa = config.orf_prediction.min_orf_length_aa
    start_codons = config.orf_prediction.start_codons
    stop_codons = config.orf_prediction.stop_codons
    table = config.orf_prediction.translation_table
    overlap_bp = config.orf_prediction.overlap_threshold_bp
    
    # Overlap resolution policy
    # We load custom configurations if defined, falling back to keep_all_flag (scientifically robust for viruses)
    policy = getattr(config.orf_prediction, "overlap_resolution_policy", "keep_all_flag")
    resolve_nested_flag = getattr(config.orf_prediction, "resolve_nested", True)
    
    seq_id = record.id
    fwd_seq = record.sequence
    
    # Compute sequence GC content for confidence scorer
    g_count = fwd_seq.upper().count('G')
    c_count = fwd_seq.upper().count('C')
    bg_gc = (g_count + c_count) / len(fwd_seq) if len(fwd_seq) > 0 else 0.5
    
    # 1. Scan Forward Strand (+1, +2, +3 frames)
    logger.info(f"Scanning forward strand of '{seq_id}' in 3 reading frames...")
    fwd_orfs = scan_single_strand(
        seq_id, fwd_seq, "+", start_codons, stop_codons, min_len_aa, table
    )
    
    # 2. Get Reverse Complement and Scan (-1, -2, -3 frames)
    logger.info(f"Scanning reverse complement strand of '{seq_id}'...")
    biopy_seq = Seq(fwd_seq)
    rev_seq = str(biopy_seq.reverse_complement())
    rev_orfs = scan_single_strand(
        seq_id, rev_seq, "-", start_codons, stop_codons, min_len_aa, table
    )
    
    all_raw_orfs = fwd_orfs + rev_orfs
    logger.info(f"Identified {len(all_raw_orfs)} total raw candidate ORFs before overlap filtering.")
    
    # Compute in silico confidence scores
    for orf in all_raw_orfs:
        orf.confidence_score = calculate_confidence(
            orf.length_aa, orf.start_codon, orf.nuc_seq, bg_gc
        )
        
    # 3. Resolve overlaps scientifically
    resolved_orfs = resolve_overlaps(all_raw_orfs, overlap_bp, policy)
    
    # 4. Filter or flag nested ORFs
    final_orfs = filter_nested_orfs(resolved_orfs, resolve_nested_flag)
    
    # 5. Standardize and rename final ORF IDs sequentially
    for idx, orf in enumerate(final_orfs):
        orf.id = f"{seq_id}_predicted_orf_{idx+1}"
        
    logger.info(f"Retained {len(final_orfs)} high-confidence ORFs for sequence '{seq_id}'.")
    return final_orfs


def process_orf_prediction(
    cleaned_fasta: Path,
    outdir: Path,
    config: Any
) -> Dict[str, Any]:
    """
    Orchestrates the ORF prediction workflow across all sequences in the preprocessed FASTA.
    Saves GFF3 files, translated FASTA files, and comprehensive statistics reports.
    """
    cleaned_fasta = Path(cleaned_fasta)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    output_gff = outdir / "coordinates.gff3"
    output_fasta = outdir / "proteins.fasta"
    report_json = outdir / "orf_statistics.json"
    coords_csv = outdir / "coordinates.csv"
    literature_json = outdir / "orf_literature_comparison.json"
    
    logger.info(f"Reading preprocessed sequence: {cleaned_fasta}")
    
    records = []
    for rec in SeqIO.parse(str(cleaned_fasta), "fasta"):
        records.append(SequenceRecord(rec.id, str(rec.seq), rec.description))
        
    all_predicted_orfs: List[ORFRecord] = []
    literature_comparisons = []
    
    for rec in records:
        predicted = predict_orfs_in_sequence(rec, config)
        all_predicted_orfs.extend(predicted)
        
        # Comparative literature check
        comp_res = compare_predicted_to_literature(predicted, len(rec.sequence), rec.id)
        literature_comparisons.append(comp_res)
        
    # Write literature comparison results
    with open(literature_json, "w", encoding="utf-8") as f:
        json.dump(literature_comparisons, f, indent=4)
    logger.info(f"Saved predicted ORF literature comparisons to: {literature_json}")
        
    # Compile output tables
    gff_lines = ["##gff-version 3"]
    csv_rows = []
    protein_records = []
    
    lengths_aa = []
    start_codon_counts = Counter()
    strand_counts = Counter()
    frame_counts = Counter()
    
    for orf in all_predicted_orfs:
        gff_lines.append(orf.to_gff_line())
        csv_rows.append(orf.to_dict())
        
        # Create protein FASTA record
        prot_desc = f"location={orf.strand}{orf.start}..{orf.end} length_aa={orf.length_aa} start_codon={orf.start_codon}"
        protein_records.append(SeqRecord(Seq(orf.prot_seq), id=orf.id, description=prot_desc))
        
        # Accumulate metrics
        lengths_aa.append(orf.length_aa)
        start_codon_counts[orf.start_codon] += 1
        strand_counts[orf.strand] += 1
        frame_counts[orf.frame] += 1

    # 1. Write GFF3 file
    with open(output_gff, "w", encoding="utf-8") as f:
        f.write("\n".join(gff_lines) + "\n")
    logger.info(f"Saved predicted ORF GFF3 coordinates to: {output_gff}")

    # 2. Write CSV coordinate file
    df_coords = pd.DataFrame(csv_rows) if csv_rows else pd.DataFrame(columns=[
        "orf_id", "sequence_id", "start", "end", "strand", "frame", "length_bp", "length_aa", "start_codon", "confidence_score", "overlap_flag", "overlap_with"
    ])
    df_coords.to_csv(coords_csv, index=False)
    logger.info(f"Saved predicted ORF CSV coordinates to: {coords_csv}")

    # 3. Write protein translated FASTA file
    if protein_records:
        SeqIO.write(protein_records, str(output_fasta), "fasta")
        logger.info(f"Saved translated proteins to: {output_fasta}")
    else:
        # Write empty FASTA
        with open(output_fasta, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("No ORFs predicted. Empty proteins.fasta created.")

    # Calculate statistics summary
    total_orfs = len(all_predicted_orfs)
    mean_len_aa = sum(lengths_aa) / total_orfs if total_orfs > 0 else 0
    overlap_count = sum(1 for orf in all_predicted_orfs if orf.overlap_flag == "Overlap")
    nested_count = sum(1 for orf in all_predicted_orfs if orf.overlap_flag == "Nested")
    mean_conf = sum(orf.confidence_score for orf in all_predicted_orfs) / total_orfs if total_orfs > 0 else 0.0
    
    stats_summary = {
        "meta": {
            "pipeline": config.pipeline.name,
            "version": config.pipeline.version,
            "min_orf_length_aa_filter": config.orf_prediction.min_orf_length_aa,
            "translation_table": config.orf_prediction.translation_table
        },
        "counts": {
            "total_orfs_predicted": total_orfs,
            "strand_distribution": dict(strand_counts),
            "frame_distribution": dict(frame_counts),
            "start_codon_frequencies": dict(start_codon_counts),
            "overlapping_orfs_flagged": overlap_count,
            "nested_orfs_flagged": nested_count
        },
        "metrics": {
            "min_length_aa": min(lengths_aa) if lengths_aa else 0,
            "max_length_aa": max(lengths_aa) if lengths_aa else 0,
            "mean_length_aa": round(mean_len_aa, 2),
            "mean_confidence_score": round(mean_conf, 4)
        },
        "output_files": {
            "gff3_coordinates": str(output_gff),
            "csv_coordinates": str(coords_csv),
            "proteins_fasta": str(output_fasta),
            "statistics_json": str(report_json),
            "literature_comparison_json": str(literature_json)
        }
    }
    
    # 4. Save statistics report json
    with open(report_json, "w", encoding="utf-8") as f:
        json.dump(stats_summary, f, indent=4)
    logger.info(f"Saved predicted ORF statistics report to: {report_json}")
    
    return stats_summary
