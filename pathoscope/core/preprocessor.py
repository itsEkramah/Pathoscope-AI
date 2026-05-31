import os
import re
import sys
import gzip
import time
import json
import hashlib
import subprocess
import multiprocessing
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional, Generator
from collections import Counter
import numpy as np
import pandas as pd
from loguru import logger
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Custom exceptions for preprocessing
class PreprocessingError(Exception):
    """Base exception for preprocessing stage errors."""
    pass

class ToolExecutionError(PreprocessingError):
    """Raised when FastQC or fastp binary fails to execute."""
    pass

class InvalidFileFormatError(PreprocessingError):
    """Raised when file extension or content layout is unsupported."""
    pass

class SequenceValidationError(PreprocessingError):
    """Raised when critical sequencing validation checks fail."""
    pass


class FastqRecord:
    """
    Standardized internal representation of a FASTQ read record.
    Provides fast, Phred-score parsing and validation calculations.
    """
    def __init__(self, header: str, sequence: str, plus: str, quality: str):
        self.header = header.strip()
        self.sequence = sequence.strip().upper()
        self.plus = plus.strip()
        self.quality = quality.strip()
        
        if not self.header:
            raise SequenceValidationError("FASTQ record missing header line.")
        if len(self.sequence) != len(self.quality):
            raise SequenceValidationError(
                f"FASTQ sequence length ({len(self.sequence)}) does not match quality length ({len(self.quality)}) for read '{self.header}'"
            )

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def phred_scores(self) -> List[int]:
        """Convert ASCII quality characters to Phred quality integers (Q = ord(char) - 33)."""
        return [ord(c) - 33 for c in self.quality]

    @property
    def mean_quality(self) -> float:
        """Calculates the average Phred quality score (Q-score) for this read."""
        scores = self.phred_scores
        return sum(scores) / len(scores) if scores else 0.0

    def to_fastq_string(self) -> str:
        """Formats the record as a standard 4-line FASTQ string."""
        return f"{self.header}\n{self.sequence}\n{self.plus}\n{self.quality}\n"


def open_sequence_stream(file_path: Path, mode: str = "rt"):
    """Gracefully opens raw or gzipped (.gz) text streams."""
    file_path = Path(file_path)
    if str(file_path).endswith(".gz"):
        return gzip.open(file_path, mode, encoding="utf-8")
    return open(file_path, mode, encoding="utf-8")


def parse_fastq_stream(stream) -> Generator[FastqRecord, None, None]:
    """Efficiently streams 4-line FASTQ records from an open file descriptor."""
    while True:
        h = stream.readline()
        if not h:
            break
        s = stream.readline()
        p = stream.readline()
        q = stream.readline()
        
        if not s or not p or not q:
            raise SequenceValidationError("Truncated FASTQ record format detected.")
            
        yield FastqRecord(h, s, p, q)


def trim_fastq_adapter(record: FastqRecord, adapter: str) -> Tuple[FastqRecord, bool]:
    """
    Scans a read for the configured adapter sequence and clips
    it along with its corresponding quality string if found.
    """
    if not adapter:
        return record, False
    
    seq = record.sequence
    idx = seq.find(adapter.upper())
    if idx != -1:
        # Trim from adapter index onwards
        trimmed_seq = seq[:idx]
        trimmed_qual = record.quality[:idx]
        
        # If read is empty after trimming, yield empty base placeholder to be caught by length filters
        if not trimmed_seq:
            trimmed_seq = "N"
            trimmed_qual = "!"
            
        return FastqRecord(record.header, trimmed_seq, record.plus, trimmed_qual), True
        
    return record, False


def sliding_window_quality_trim(
    record: FastqRecord, 
    window_size: int = 4, 
    min_qual: int = 20
) -> Tuple[FastqRecord, bool]:
    """
    purpose:
        Trims the 3' end of a FASTQ read when the average Phred quality score
        within a sliding window drops below a specified threshold.
    inputs:
        - record (FastqRecord): The raw sequence record with quality scores.
        - window_size (int): Size of the sliding window (default: 4).
        - min_qual (int): Minimum average quality score required (default: 20).
    outputs:
        - Tuple[FastqRecord, bool]: The quality-trimmed record and a boolean
          indicating whether trimming occurred.
    biological rationale:
        Sequencing quality typically deteriorates toward the 3' end of reads due to
        fluorophore degradation and phasing issues. Trimming low-quality tails reduces
        errors in downstream processes such as gene finding and annotation.
    """
    seq = record.sequence
    qual = record.quality
    
    # Vectorized computation of Phred quality scores
    scores = np.frombuffer(qual.encode('ascii'), dtype=np.uint8) - 33
    n = len(scores)
    
    if n < window_size:
        return record, False
        
    # Calculate rolling average using vectorized NumPy operations
    cumsum = np.cumsum(np.insert(scores, 0, 0))
    window_sums = cumsum[window_size:] - cumsum[:-window_size]
    window_means = window_sums / window_size
    
    # Scan backwards from the 3' end to find low-quality regions
    low_qual = window_means < min_qual
    
    if not low_qual[-1]:
        trim_idx = n
    else:
        high_indices = np.where(~low_qual)[0]
        if len(high_indices) == 0:
            trim_idx = 0
        else:
            last_high_idx = high_indices[-1]
            trim_idx = last_high_idx + 1
        
    if trim_idx < n:
        trimmed_seq = seq[:trim_idx]
        trimmed_qual = qual[:trim_idx]
        
        if not trimmed_seq:
            trimmed_seq = "N"
            trimmed_qual = "!"
            
        return FastqRecord(record.header, trimmed_seq, record.plus, trimmed_qual), True
        
    return record, False


def calculate_n50(lengths: List[int]) -> int:
    """Calculates the N50 metric of a list of sequence lengths."""
    if not lengths:
        return 0
    sorted_lengths = sorted(lengths, reverse=True)
    total_len = sum(sorted_lengths)
    half_len = total_len / 2.0
    
    cumulative_sum = 0
    for length in sorted_lengths:
        cumulative_sum += length
        if cumulative_sum >= half_len:
            return length
    return sorted_lengths[-1]


def execute_command_wrapper(cmd: List[str]) -> subprocess.CompletedProcess:
    """Safely executes subprocess execution commands, tracking errors."""
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            shell=False,
            check=True
        )
        return result
    except FileNotFoundError:
        raise ToolExecutionError(f"Bioinformatics tool binary '{cmd[0]}' not found in system environment PATH.")
    except subprocess.CalledProcessError as cpe:
        raise ToolExecutionError(f"Command execution of '{cmd[0]}' failed (Exit Code {cpe.returncode}): {cpe.stderr}")


def run_fastp_subprocess(
    in_r1: Path,
    out_r1: Path,
    in_r2: Optional[Path],
    out_r2: Optional[Path],
    out_fasta: Path,
    config: Any,
    parent_dir: Path
) -> Tuple[Dict[str, Any], str]:
    """
    Subprocess orchestrator wrapper for the industry-standard 'fastp' tool.
    Yields cleaned fastq output and automatically aggregates metrics to JSON.
    """
    report_json = parent_dir / "fastp_report.json"
    report_html = parent_dir / "fastp_report.html"
    
    cmd = [
        "fastp",
        "--in1", str(in_r1),
        "--out1", str(out_r1),
        "--json", str(report_json),
        "--html", str(report_html),
        "--length_required", str(config.preprocessing.min_length),
        "--adapter_sequence", config.preprocessing.adapter_forward,
        "--thread", "2"
    ]
    
    if in_r2 and out_r2:
        cmd.extend([
            "--in2", str(in_r2),
            "--out2", str(out_r2),
            "--adapter_sequence_r2", config.preprocessing.adapter_reverse
        ])
        
    logger.info(f"Executing fastp subprocess process: {' '.join(cmd)}")
    start_time = time.time()
    execute_command_wrapper(cmd)
    elapsed = time.time() - start_time
    
    # Parse fastp outputs
    with open(report_json, "r", encoding="utf-8") as f:
        fastp_data = json.load(f)
        
    # Convert trimmed FASTQ output to FASTA format for downstream ORF Prediction seamlessly
    logger.info("Converting fastp cleaned FASTQ output to FASTA for ORF Prediction phase...")
    records = []
    with open_sequence_stream(out_r1) as stream:
        for rec in parse_fastq_stream(stream):
            records.append(SeqRecord(Seq(rec.sequence), id=rec.header[1:].split()[0], description=""))
    SeqIO.write(records, str(out_fasta), "fasta")
    
    version = fastp_data.get("meta", {}).get("fastp_version", "unknown")
    command_str = " ".join(cmd)
    
    # Standardize output metrics
    summary = {
        "software": "fastp",
        "version": version,
        "command_history": command_str,
        "runtime_seconds": elapsed,
        "counts": {
            "total_processed": fastp_data["summary"]["before_filtering"]["total_reads"],
            "total_kept": fastp_data["summary"]["after_filtering"]["total_reads"],
            "total_discarded": fastp_data["summary"]["before_filtering"]["total_reads"] - fastp_data["summary"]["after_filtering"]["total_reads"],
            "adapter_trimmed_reads": fastp_data["filtering_result"].get("adapter_trimmed_reads", 0)
        },
        "metrics": {
            "total_bases": fastp_data["summary"]["after_filtering"]["total_bases"],
            "mean_length": fastp_data["summary"]["after_filtering"]["read1_mean_length"],
            "mean_gc_percent": round(fastp_data["summary"]["after_filtering"]["gc_content"] * 100.0, 2) if fastp_data["summary"]["after_filtering"]["gc_content"] < 1.0 else fastp_data["summary"]["after_filtering"]["gc_content"]
        }
    }
    
    return summary, command_str


def run_fastqc_subprocess(in_file: Path, parent_dir: Path):
    """Subprocess orchestrator wrapper for the Java-based 'FastQC' tool."""
    cmd = [
        "fastqc",
        "--outdir", str(parent_dir),
        "--format", "fastq",
        str(in_file)
    ]
    logger.info(f"Executing FastQC subprocess analysis: {' '.join(cmd)}")
    execute_command_wrapper(cmd)


def _process_single_read_numpy(
    record: FastqRecord,
    adapter: str,
    min_length: int,
    max_length: int,
    min_mean_qscore: float,
    max_ambig_pct: float
) -> Tuple[bool, str, Optional[FastqRecord], Dict[str, Any]]:
    """
    purpose:
        Processes a single FASTQ read by trimming adapter sequences, performing
        sliding window quality trimming, and validating length, quality, and ambiguous bases.
    inputs:
        - record (FastqRecord): The raw sequence record with quality scores.
        - adapter (str): Adapter sequence to look for and trim.
        - min_length (int): Minimum acceptable read length.
        - max_length (int): Maximum acceptable read length.
        - min_mean_qscore (float): Minimum acceptable average Phred quality score.
        - max_ambig_pct (float): Maximum acceptable percentage of ambiguous nucleotides.
    outputs:
        - Tuple[bool, str, Optional[FastqRecord], Dict[str, Any]]:
          A tuple containing validation success (bool), rejection reason (str),
          the processed FastqRecord (or None if invalid), and a dictionary of metrics/stats.
    biological rationale:
        Trimming artificial sequence adapters and filtering out low-quality or highly
        ambiguous sequences is standard bioinformatics practice to prevent false-positive annotations
        and erroneous secondary alignments.
    """
    # 1. Trimming Adapters
    trimmed_record, was_trimmed = trim_fastq_adapter(record, adapter)
    
    # 2. Sliding Window Quality Trimming
    trimmed_record, was_q_trimmed = sliding_window_quality_trim(trimmed_record, window_size=4, min_qual=20)
    
    length = trimmed_record.length
    if length == 0:
        return False, f"Read length too short (0bp < {min_length}bp cutoff)", None, {
            "read_id": record.header.split()[0],
            "length": 0,
            "mean_qscore": 0.0,
            "gc_percent": 0.0,
            "was_trimmed": was_trimmed,
            "was_q_trimmed": was_q_trimmed,
            "sequence": "",
            "quality": ""
        }
        
    # Enforce quality score evaluation via vectorized array operations (NumPy arrays)
    scores = np.frombuffer(trimmed_record.quality.encode('ascii'), dtype=np.uint8) - 33
    mean_q = float(np.mean(scores)) if len(scores) > 0 else 0.0
    
    seq_bytes = trimmed_record.sequence.encode('ascii')
    seq_arr = np.frombuffer(seq_bytes, dtype=np.uint8)
    
    # Ambiguous count: N is ASCII 78
    ambig_count = int(np.sum(seq_arr == 78))
    ambig_pct = (ambig_count / length) * 100.0 if length > 0 else 0.0
    
    # GC count: G is 71, C is 67
    gc_count = int(np.sum((seq_arr == 71) | (seq_arr == 67)))
    gc_pct = (gc_count / length) * 100.0 if length > 0 else 0.0
    
    is_valid = True
    reason = ""
    
    if length < min_length:
        is_valid = False
        reason = f"Read length too short ({length}bp < {min_length}bp cutoff)"
    elif length > max_length:
        is_valid = False
        reason = f"Read length too long ({length}bp > {max_length}bp cutoff)"
    elif mean_q < min_mean_qscore:
        is_valid = False
        reason = f"Low average quality score (Q{mean_q:.1f} < Q{min_mean_qscore})"
    elif ambig_pct > max_ambig_pct:
        is_valid = False
        reason = f"Ambiguous bases limit exceeded ({ambig_pct:.1f}% > {max_ambig_pct}%)"
        
    stats = {
        "read_id": trimmed_record.header.split()[0],
        "length": length,
        "mean_qscore": round(mean_q, 1),
        "gc_percent": round(gc_pct, 2),
        "was_trimmed": was_trimmed,
        "was_q_trimmed": was_q_trimmed,
        "sequence": trimmed_record.sequence,
        "quality": trimmed_record.quality
    }
    
    return is_valid, reason, trimmed_record if is_valid else None, stats


def _process_single_end_chunk_worker(args) -> List[Tuple[bool, str, Optional[FastqRecord], Dict[str, Any]]]:
    """
    purpose:
        Worker function to process a chunk of single-end reads.
    inputs:
        - args (tuple): Contains (records, adapter, min_length, max_length, min_mean_qscore, max_ambig_pct).
    outputs:
        - List[Tuple[bool, str, Optional[FastqRecord], Dict[str, Any]]]: List of processing results for each read.
    biological rationale:
        Parallelizes the CPU-intensive QC and trimming operations across multiple CPU cores,
        greatly reducing computational bottlenecking for large FASTQ sequencing files.
    """
    records, adapter, min_length, max_length, min_mean_qscore, max_ambig_pct = args
    results = []
    for record in records:
        res = _process_single_read_numpy(record, adapter, min_length, max_length, min_mean_qscore, max_ambig_pct)
        results.append(res)
    return results


def _process_paired_end_chunk_worker(args) -> List[Tuple[bool, str, Optional[Tuple[FastqRecord, FastqRecord]], Dict[str, Any]]]:
    """
    purpose:
        Worker function to process a chunk of paired-end reads.
    inputs:
        - args (tuple): Contains (pairs, f_adapter, r_adapter, min_length, max_length, min_mean_qscore, max_ambig_pct).
    outputs:
        - List[Tuple[bool, str, Optional[Tuple[FastqRecord, FastqRecord]], Dict[str, Any]]]: List of processing results for each pair.
    biological rationale:
        Processes paired-end reads concurrently while strictly enforcing that both forward (R1)
        and reverse (R2) reads pass filters to maintain paired-end coordinate synchronization.
    """
    pairs, f_adapter, r_adapter, min_length, max_length, min_mean_qscore, max_ambig_pct = args
    results = []
    for r1, r2 in pairs:
        v1, re1, t1, s1 = _process_single_read_numpy(r1, f_adapter, min_length, max_length, min_mean_qscore, max_ambig_pct)
        v2, re2, t2, s2 = _process_single_read_numpy(r2, r_adapter, min_length, max_length, min_mean_qscore, max_ambig_pct)
        
        is_valid = v1 and v2
        
        # Discard pair if either is invalid (synchronicity constraint)
        reason = ""
        if not is_valid:
            if not v1 and not v2:
                reason = f"R1: {re1} | R2: {re2}"
            elif not v1:
                reason = f"R1: {re1}"
            else:
                reason = f"R2: {re2}"
                
        stats = {
            "read_id": r1.header.split()[0],
            "r1_qscore": s1["mean_qscore"],
            "r2_qscore": s2["mean_qscore"],
            "r1_len": s1["length"],
            "r2_len": s2["length"],
            "r1_gc": s1["gc_percent"],
            "r2_gc": s2["gc_percent"],
            "was_trimmed": s1["was_trimmed"] or s2["was_trimmed"],
            "r1_seq": s1["sequence"],
            "r2_seq": s2["sequence"]
        }
        
        results.append((is_valid, reason, (t1, t2) if is_valid else None, stats))
    return results


def compile_pure_python_high_fidelity_preprocessing(
    in_r1: Path,
    out_r1: Path,
    in_r2: Optional[Path],
    out_r2: Optional[Path],
    out_fasta: Path,
    config: Any,
    parent_dir: Path
) -> Tuple[Dict[str, Any], str]:
    """
    purpose:
        Orchestrates high-fidelity quality control, adapter trimming, quality sliding window cut,
        and paired-end sync checks in pure Python using parallel multiprocessing chunk queues.
    inputs:
        - in_r1 (Path): Path to raw R1 FASTQ.
        - out_r1 (Path): Path to output cleaned R1 FASTQ.
        - in_r2 (Optional[Path]): Path to raw R2 FASTQ.
        - out_r2 (Optional[Path]): Path to output cleaned R2 FASTQ.
        - out_fasta (Path): Output FASTA file path.
        - config (Any): Config configuration namespace.
        - parent_dir (Path): Outdir workspace parent directory.
    outputs:
        - Tuple[Dict[str, Any], str]: Summary dict of the run metadata and execution command history.
    biological rationale:
        Ensures high-fidelity, reproducible preprocessing of sequencing reads even in environments
        lacking fastp binaries. Features high-performance NumPy parsing and multiprocessing chunks
        to prevent sequential runtime timeouts.
    """
    logger.info("Executing Pure-Python High-Fidelity preprocessing fallback engine (No local binaries detected).")
    start_time = time.time()
    
    min_len = config.preprocessing.min_length
    max_len = config.preprocessing.max_length
    max_ambig_pct = config.preprocessing.max_ambiguous_pct
    min_qscore = config.preprocessing.min_mean_qscore
    max_reads_cap = config.preprocessing.max_reads_cap
    
    f_adapter = config.preprocessing.adapter_forward
    r_adapter = config.preprocessing.adapter_reverse
    
    # 1. Parse raw streams up to the max_reads_cap (sub-sampling limits)
    records_r1 = []
    records_r2 = []
    
    if not in_r2:
        logger.info(f"Streaming single-end read file up to max_reads_cap={max_reads_cap}: {in_r1}")
        with open_sequence_stream(in_r1) as instream:
            for record in parse_fastq_stream(instream):
                records_r1.append(record)
                if len(records_r1) >= max_reads_cap:
                    break
    else:
        logger.info(f"Streaming paired-end read files up to max_reads_cap={max_reads_cap}: R1={in_r1} | R2={in_r2}")
        with open_sequence_stream(in_r1) as stream1, open_sequence_stream(in_r2) as stream2:
            iter1 = parse_fastq_stream(stream1)
            iter2 = parse_fastq_stream(stream2)
            while len(records_r1) < max_reads_cap // 2:
                try:
                    r1 = next(iter1)
                    r2 = next(iter2)
                    records_r1.append(r1)
                    records_r2.append(r2)
                except StopIteration:
                    break
                    
    # 2. Implement multithreaded chunk-based queue processing (multiprocessing pool splits)
    use_multiprocessing = True
    num_records = len(records_r1) if not in_r2 else len(records_r1) * 2
    
    # Avoid overhead if there are very few sequences
    if num_records < 100:
        use_multiprocessing = False
        
    chunk_results = []
    
    if use_multiprocessing:
        try:
            num_workers = min(multiprocessing.cpu_count(), 4)
            if not in_r2:
                chunk_size = max(1, len(records_r1) // num_workers)
                chunks = [records_r1[i:i + chunk_size] for i in range(0, len(records_r1), chunk_size)]
                chunk_args = [(chunk, f_adapter, min_len, max_len, min_qscore, max_ambig_pct) for chunk in chunks]
                with multiprocessing.Pool(processes=num_workers) as pool:
                    chunk_results = pool.map(_process_single_end_chunk_worker, chunk_args)
            else:
                pairs = list(zip(records_r1, records_r2))
                chunk_size = max(1, len(pairs) // num_workers)
                chunks = [pairs[i:i + chunk_size] for i in range(0, len(pairs), chunk_size)]
                chunk_args = [(chunk, f_adapter, r_adapter, min_len, max_len, min_qscore, max_ambig_pct) for chunk in chunks]
                with multiprocessing.Pool(processes=num_workers) as pool:
                    chunk_results = pool.map(_process_paired_end_chunk_worker, chunk_args)
        except Exception as mp_err:
            logger.warning(f"Multiprocessing execution failed: {mp_err}. Falling back to sequential processing.")
            use_multiprocessing = False
            
    if not use_multiprocessing:
        if not in_r2:
            single_chunk_res = _process_single_end_chunk_worker((records_r1, f_adapter, min_len, max_len, min_qscore, max_ambig_pct))
            chunk_results = [single_chunk_res]
        else:
            pairs = list(zip(records_r1, records_r2))
            single_chunk_res = _process_paired_end_chunk_worker((pairs, f_adapter, r_adapter, min_len, max_len, min_qscore, max_ambig_pct))
            chunk_results = [single_chunk_res]
            
    # 3. Consolidate results and save files
    total_processed = 0
    total_kept = 0
    total_discarded = 0
    
    rejection_reasons = Counter()
    lengths = []
    gc_contents = []
    base_counts = Counter()
    adapter_trimmed_count = 0
    audit_data = []
    
    if not in_r2:
        with open_sequence_stream(out_r1, "wt") as outstream_fq, \
             open(out_fasta, "w", encoding="utf-8") as outstream_fa:
             
            for chunk in chunk_results:
                for is_valid, reason, trimmed_record, stats in chunk:
                    total_processed += 1
                    if stats["was_trimmed"]:
                        adapter_trimmed_count += 1
                        
                    audit_data.append({
                        "read_id": stats["read_id"],
                        "length": stats["length"],
                        "mean_qscore": stats["mean_qscore"],
                        "gc_percent": stats["gc_percent"],
                        "status": "ACCEPTED" if is_valid else "REJECTED",
                        "reason": reason
                    })
                    
                    if not is_valid:
                        total_discarded += 1
                        rejection_reasons[reason] += 1
                        continue
                        
                    total_kept += 1
                    lengths.append(stats["length"])
                    gc_contents.append(stats["gc_percent"])
                    base_counts.update(stats["sequence"])
                    
                    outstream_fq.write(trimmed_record.to_fastq_string())
                    clean_id = trimmed_record.header[1:].split()[0]
                    outstream_fa.write(f">{clean_id}\n{trimmed_record.sequence}\n")
    else:
        with open_sequence_stream(out_r1, "wt") as out1, \
             open_sequence_stream(out_r2, "wt") as out2, \
             open(out_fasta, "w", encoding="utf-8") as outstream_fa:
             
            for chunk in chunk_results:
                for is_valid, reason, trimmed_pair, stats in chunk:
                    total_processed += 2
                    if stats["was_trimmed"]:
                        adapter_trimmed_count += 2
                        
                    audit_data.append({
                        "read_id": stats["read_id"],
                        "status": "ACCEPTED" if is_valid else "REJECTED",
                        "r1_qscore": stats["r1_qscore"],
                        "r2_qscore": stats["r2_qscore"],
                        "reason": reason
                    })
                    
                    if not is_valid:
                        total_discarded += 2
                        if "R1:" in reason:
                            r1_r = reason.split(" | ")[0].replace("R1: ", "")
                            if r1_r:
                                rejection_reasons[r1_r] += 1
                        if "R2:" in reason:
                            r2_r = reason.split(" | ")[-1].replace("R2: ", "")
                            if r2_r:
                                rejection_reasons[r2_r] += 1
                        continue
                        
                    total_kept += 2
                    t1, t2 = trimmed_pair
                    lengths.extend([stats["r1_len"], stats["r2_len"]])
                    gc_contents.extend([stats["r1_gc"], stats["r2_gc"]])
                    base_counts.update(stats["r1_seq"])
                    base_counts.update(stats["r2_seq"])
                    
                    out1.write(t1.to_fastq_string())
                    out2.write(t2.to_fastq_string())
                    
                    clean_id = t1.header[1:].split()[0]
                    outstream_fa.write(f">{clean_id}\n{t1.sequence}\n")
                    
    elapsed = time.time() - start_time
    mean_len = sum(lengths) / len(lengths) if lengths else 0
    mean_gc = sum(gc_contents) / len(gc_contents) if gc_contents else 0
    total_bases = sum(lengths)
    
    # Save the audit statistics table
    stats_csv_path = parent_dir / "audit_stats.csv"
    pd.DataFrame(audit_data).to_csv(stats_csv_path, index=False)
    
    command_str = f"python-fastq-preprocessor-fallback --min_qscore {min_qscore} --min_length {min_len} --max_reads_cap {max_reads_cap}"
    
    summary = {
        "software": "Python-fastp Fallback",
        "version": "2.0.0",
        "command_history": command_str,
        "runtime_seconds": elapsed,
        "counts": {
            "total_processed": total_processed,
            "total_kept": total_kept,
            "total_discarded": total_discarded,
            "adapter_trimmed_reads": adapter_trimmed_count,
            "detailed_rejections": dict(rejection_reasons)
        },
        "metrics": {
            "total_bases": total_bases,
            "mean_length": round(mean_len, 2),
            "n50": calculate_n50(lengths),
            "mean_gc_percent": round(mean_gc, 2)
        },
        "base_frequencies": dict(base_counts)
    }
    
    return summary, command_str


def generate_premium_qc_html_report(
    summary: Dict[str, Any],
    output_html_path: Path
):
    """
    Renders an extremely premium, responsive, dark-mode MultiQC-like aggregation dashboard
    detailing raw statistics, GC ratios, and duplication/rejection distributions.
    """
    counts = summary["counts"]
    metrics = summary["metrics"]
    
    keep_pct = (counts["total_kept"] / counts["total_processed"]) * 100.0 if counts["total_processed"] > 0 else 0.0
    disc_pct = 100.0 - keep_pct
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PathoScope AI - Quality Control Analysis Dashboard</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700;800&family=Outfit:wght@400;600;800&display=swap" rel="stylesheet">
    <style>
        :root {{
            --bg-main: #0b0f19;
            --bg-card: rgba(17, 24, 39, 0.7);
            --border-card: rgba(255, 255, 255, 0.08);
            --text-main: #f3f4f6;
            --text-muted: #9ca3af;
            --accent: #22d3ee;
            --accent-green: #34d399;
            --accent-red: #f87171;
            --gradient: linear-gradient(135deg, #1e1b4b 0%, #0b0f19 100%);
        }}
        
        body {{
            background-color: var(--bg-main);
            color: var(--text-main);
            font-family: 'Inter', sans-serif;
            margin: 0;
            padding: 40px 20px;
            background-image: var(--gradient);
            background-attachment: fixed;
        }}
        
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        
        .hero-banner {{
            background: rgba(255, 255, 255, 0.02);
            border: 1px solid var(--border-card);
            border-radius: 20px;
            padding: 40px;
            text-align: center;
            margin-bottom: 40px;
            backdrop-filter: blur(12px);
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.5);
        }}
        
        .hero-banner h1 {{
            font-family: 'Outfit', sans-serif;
            font-size: 3rem;
            margin: 0;
            font-weight: 800;
            background: linear-gradient(90deg, #c084fc, var(--accent));
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }}
        
        .hero-banner p {{
            color: var(--text-muted);
            font-size: 1.15rem;
            margin-top: 10px;
            margin-bottom: 0;
        }}
        
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 25px;
            margin-bottom: 40px;
        }}
        
        .card {{
            background: var(--bg-card);
            border: 1px solid var(--border-card);
            border-radius: 16px;
            padding: 25px;
            backdrop-filter: blur(8px);
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            transition: transform 0.2s ease, border-color 0.2s ease;
        }}
        
        .card:hover {{
            transform: translateY(-2px);
            border-color: var(--accent);
        }}
        
        .card-header {{
            font-size: 0.9rem;
            text-transform: uppercase;
            letter-spacing: 1px;
            color: var(--text-muted);
            margin-bottom: 15px;
        }}
        
        .card-value {{
            font-family: 'Outfit', sans-serif;
            font-size: 2.5rem;
            font-weight: 800;
            color: var(--accent);
        }}
        
        .card-subvalue {{
            font-size: 0.85rem;
            color: var(--text-muted);
            margin-top: 5px;
        }}
        
        .split-grid {{
            display: grid;
            grid-template-columns: 2fr 1fr;
            gap: 30px;
            margin-bottom: 40px;
        }}
        
        @media (max-width: 900px) {{
            .split-grid {{
                grid-template-columns: 1fr;
            }}
        }}
        
        .table-container {{
            background: var(--bg-card);
            border: 1px solid var(--border-card);
            border-radius: 16px;
            padding: 30px;
            backdrop-filter: blur(8px);
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
        }}
        
        .table-container h2 {{
            font-family: 'Outfit', sans-serif;
            font-size: 1.6rem;
            margin-top: 0;
            margin-bottom: 25px;
            font-weight: 800;
            border-left: 5px solid var(--accent);
            padding-left: 15px;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
        }}
        
        th, td {{
            text-align: left;
            padding: 14px 10px;
            border-bottom: 1px solid rgba(255, 255, 255, 0.05);
        }}
        
        th {{
            color: var(--text-muted);
            font-size: 0.9rem;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        
        td {{
            font-size: 1.05rem;
        }}
        
        .highlight {{
            color: var(--accent-green);
            font-weight: 600;
        }}
        
        .badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.8rem;
            font-weight: 600;
            text-transform: uppercase;
        }}
        
        .badge-success {{
            background: rgba(52, 211, 153, 0.1);
            color: var(--accent-green);
            border: 1px solid rgba(52, 211, 153, 0.2);
        }}
        
        .progress-bar-container {{
            background: rgba(255, 255, 255, 0.05);
            height: 10px;
            border-radius: 5px;
            overflow: hidden;
            margin-top: 15px;
        }}
        
        .progress-bar {{
            height: 100%;
            background: linear-gradient(90deg, var(--accent-green) 0%, var(--accent) 100%);
            border-radius: 5px;
        }}
        
        .console-box {{
            background: #020617;
            border: 1px solid rgba(255, 255, 255, 0.05);
            border-radius: 8px;
            padding: 15px;
            font-family: 'Courier New', Courier, monospace;
            font-size: 0.85rem;
            color: #38bdf8;
            overflow-x: auto;
            white-space: pre-wrap;
            word-break: break-all;
        }}
    </style>
</head>
<body>
    <div class="container">
        
        <!-- Header Banner -->
        <div class="hero-banner">
            <h1>PathoScope AI Preprocessing Dashboard</h1>
            <p>Auditable, scientifically defensible Sequence Quality Control & Preprocessing metrics</p>
        </div>
        
        <!-- Key Metrics Cards -->
        <div class="grid">
            <div class="card">
                <div class="card-header">Processed Reads</div>
                <div class="card-value">{counts["total_processed"]:,}</div>
                <div class="card-subvalue">Raw sequenced records analyzed</div>
            </div>
            <div class="card">
                <div class="card-header">Passed QC Reads</div>
                <div class="card-value" style="color: var(--accent-green);">{counts["total_kept"]:,}</div>
                <div class="card-subvalue">Accepted for gene prediction ({keep_pct:.2f}%)</div>
                <div class="progress-bar-container">
                    <div class="progress-bar" style="width: {keep_pct}%;"></div>
                </div>
            </div>
            <div class="card">
                <div class="card-header">N50 Statistics</div>
                <div class="card-value">{metrics["n50"]:,} bp</div>
                <div class="card-subvalue">Sequencing N50 length score</div>
            </div>
            <div class="card">
                <div class="card-header">Mean GC Content</div>
                <div class="card-value">{metrics["mean_gc_percent"]}%</div>
                <div class="card-subvalue">Mean nucleotide G-C bias ratio</div>
            </div>
        </div>
        
        <!-- Detailed Tables Grid -->
        <div class="split-grid">
            
            <!-- Auditing Summary Table -->
            <div class="table-container">
                <h2>Pipeline Summary</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Parameter description</th>
                            <th>Value / Status</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td>Software Utilized</td>
                            <td><span class="highlight">{summary["software"]}</span></td>
                        </tr>
                        <tr>
                            <td>Software version</td>
                            <td>v{summary["version"]}</td>
                        </tr>
                        <tr>
                            <td>Bases Generated (Cleaned)</td>
                            <td>{metrics["total_bases"]:,} bp</td>
                        </tr>
                        <tr>
                            <td>Mean Read Length</td>
                            <td>{metrics["mean_length"]} bp</td>
                        </tr>
                        <tr>
                            <td>Adapters Trimmed</td>
                            <td>{counts["adapter_trimmed_reads"]:,} reads</td>
                        </tr>
                        <tr>
                            <td>Discarded Reads</td>
                            <td style="color: var(--accent-red); font-weight: 600;">{counts["total_discarded"]:,} ({disc_pct:.2f}%)</td>
                        </tr>
                    </tbody>
                </table>
            </div>
            
            <!-- Software environment verification -->
            <div class="table-container">
                <h2>Workflow Audit</h2>
                <div class="card-header" style="margin-bottom: 10px;">Execution Command</div>
                <div class="console-box">{summary["command_history"]}</div>
                
                <div class="card-header" style="margin-top: 20px; margin-bottom: 10px;">Pipeline Metadata</div>
                <div style="font-size: 0.9rem; line-height: 1.6;">
                    <strong>Pipeline Name:</strong> {summary["meta"]["pipeline"]}<br>
                    <strong>Pipeline Version:</strong> v{summary["meta"]["version"]}<br>
                    <strong>Timestamp:</strong> {time.strftime('%Y-%m-%d %H:%M:%S')}
                </div>
            </div>
            
        </div>
        
    </div>
</body>
</html>
"""
    with open(output_html_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    logger.info(f"Successfully generated interactive HTML dashboard report at: {output_html_path}")


def validate_paired_end_synchronicity(in_r1: Path, in_r2: Path):
    """
    Validates that paired-end R1 and R2 files are fully synchronized.
    Checks that they exist, have the same sequence count, and matching read headers.
    """
    if not in_r1.exists():
        raise FileNotFoundError(f"Paired-end R1 file does not exist: {in_r1}")
    if not in_r2.exists():
        raise FileNotFoundError(f"Paired-end R2 file does not exist: {in_r2}")

    logger.info(f"Validating paired-end synchronicity between R1: {in_r1.name} and R2: {in_r2.name}...")
    
    r1_count = 0
    r2_count = 0
    
    with open_sequence_stream(in_r1) as s1, open_sequence_stream(in_r2) as s2:
        while True:
            line1 = s1.readline()
            line2 = s2.readline()
            
            if not line1 and not line2:
                break
            if (line1 and not line2) or (not line1 and line2):
                raise SequenceValidationError("Paired-end R1 and R2 files have mismatched read counts.")
                
            # FASTQ headers are lines where index % 4 == 0
            if r1_count % 4 == 0:
                h1 = line1.strip().split()[0]
                h2 = line2.strip().split()[0]
                
                # Check matching read names (allowing /1 and /2 or standard Illumina 1:N:0:1 suffix differences)
                base_h1 = h1
                if base_h1.endswith("/1") or base_h1.endswith("/2"):
                    base_h1 = base_h1[:-2]
                    
                base_h2 = h2
                if base_h2.endswith("/1") or base_h2.endswith("/2"):
                    base_h2 = base_h2[:-2]
                
                if base_h1 != base_h2:
                    raise SequenceValidationError(
                        f"Mismatched paired-end headers at record {r1_count // 4 + 1}: '{h1}' vs '{h2}'"
                    )
            
            r1_count += 1
            r2_count += 1
            
    logger.info(f"Paired-end synchronicity validated successfully: {r1_count // 4} aligned read pairs.")


def calculate_high_fidelity_qc_distributions(
    raw_r1: Path,
    raw_r2: Optional[Path],
    clean_r1: Path,
    clean_r2: Optional[Path],
    adapter_trimmed_count: int,
    rejection_reasons: Dict[str, int]
) -> Dict[str, Any]:
    """
    Ingests raw and cleaned sequence files and computes high-fidelity QC metrics:
    - Per-Base Quality Score Distribution
    - GC Content Distribution (binned 0-100)
    - Sequence Duplication Levels
    - Read Retention statistics
    """
    logger.info("Computing high-fidelity QC metrics and distributions directly from sequence files...")
    
    max_len = 1000
    qual_sums = [0.0] * max_len
    qual_counts = [0] * max_len
    
    gc_bins = {i: 0 for i in range(101)}
    seq_hashes = Counter()
    
    clean_paths = []
    if clean_r1.exists() and clean_r1.stat().st_size > 0:
        clean_paths.append(clean_r1)
    if clean_r2 and clean_r2.exists() and clean_r2.stat().st_size > 0:
        clean_paths.append(clean_r2)
        
    for p in clean_paths:
        try:
            with open_sequence_stream(p) as stream:
                for record in parse_fastq_stream(stream):
                    seq = record.sequence
                    qual_scores = record.phred_scores
                    
                    for i, q in enumerate(qual_scores):
                        if i < max_len:
                            qual_sums[i] += q
                            qual_counts[i] += 1
                            
                    length = len(seq)
                    if length > 0:
                        gc_count = sum(1 for b in seq if b in "GC")
                        gc_pct = int(round((gc_count / length) * 100.0))
                        if 0 <= gc_pct <= 100:
                            gc_bins[gc_pct] += 1
                            
                    h = hashlib.md5(seq.encode('utf-8')).hexdigest()
                    seq_hashes[h] += 1
        except Exception as e:
            logger.warning(f"Error reading file {p} for QC distributions: {e}")
            
    per_base_qualities = []
    for s, c in zip(qual_sums, qual_counts):
        if c > 0:
            per_base_qualities.append(round(s / c, 2))
        else:
            break
            
    duplication_levels = {
        "1": 0,
        "2": 0,
        "3-5": 0,
        "6-10": 0,
        "11-50": 0,
        ">50": 0
    }
    
    for h, count in seq_hashes.items():
        if count == 1:
            duplication_levels["1"] += 1
        elif count == 2:
            duplication_levels["2"] += 1
        elif 3 <= count <= 5:
            duplication_levels["3-5"] += 1
        elif 6 <= count <= 10:
            duplication_levels["6-10"] += 1
        elif 11 <= count <= 50:
            duplication_levels["11-50"] += 1
        elif count > 50:
            duplication_levels[">50"] += 1
            
    raw_reads_count = 0
    clean_reads_count = 0
    
    try:
        if raw_r1.exists() and raw_r1.stat().st_size > 0:
            with open_sequence_stream(raw_r1) as s:
                for _ in parse_fastq_stream(s):
                    raw_reads_count += 1
        if raw_r2 and raw_r2.exists() and raw_r2.stat().st_size > 0:
            with open_sequence_stream(raw_r2) as s:
                for _ in parse_fastq_stream(s):
                    raw_reads_count += 1
                    
        if clean_r1.exists() and clean_r1.stat().st_size > 0:
            with open_sequence_stream(clean_r1) as s:
                for _ in parse_fastq_stream(s):
                    clean_reads_count += 1
        if clean_r2 and clean_r2.exists() and clean_r2.stat().st_size > 0:
            with open_sequence_stream(clean_r2) as s:
                for _ in parse_fastq_stream(s):
                    clean_reads_count += 1
    except Exception as e:
        logger.warning(f"Error calculating read counts: {e}")
        
    waterfall = {
        "raw_reads": raw_reads_count,
        "adapter_trimmed": adapter_trimmed_count,
        "length_filtered": rejection_reasons.get("Read length too short...", 0) + rejection_reasons.get("Read length too long...", 0),
        "quality_filtered": rejection_reasons.get("Low average quality...", 0),
        "kept_reads": clean_reads_count
    }
    
    discarded = raw_reads_count - clean_reads_count
    if discarded < 0:
        discarded = 0
    if sum([waterfall["length_filtered"], waterfall["quality_filtered"]]) == 0:
        waterfall["quality_filtered"] = discarded
        
    return {
        "per_base_quality": per_base_qualities,
        "gc_content_distribution": {str(k): v for k, v in gc_bins.items()},
        "sequence_duplication_levels": duplication_levels,
        "read_retention_waterfall": waterfall
    }


def process_sequences(
    input_file: Path,
    output_fasta: Path,
    config: Any
) -> Dict[str, Any]:
    """
    Orchestrates the modular genomics preprocessing pipeline workflow.
    Checks available local binaries (FastQC, fastp) and executes them via subprocesses.
    On binary missing exceptions, seamlessly falls back to high-fidelity, auditable
    pure-Python adapter trimming, Phred cuts, and paired-end synchronization.
    """
    input_file = Path(input_file)
    output_fasta = Path(output_fasta)
    parent_dir = output_fasta.parent
    parent_dir.mkdir(parents=True, exist_ok=True)
    
    # Establish target output file naming paths
    out_r1 = parent_dir / "cleaned_R1.fastq.gz"
    out_r2 = parent_dir / "cleaned_R2.fastq.gz"
    report_html_path = parent_dir / "qc_report.html"
    
    # Extract file configurations
    paired = config.preprocessing.paired_end
    r2_path = config.preprocessing.fastq_r2_path
    
    in_r1 = input_file
    in_r2 = Path(r2_path) if r2_path else None
    
    # Detect file suffixes
    suffix = in_r1.name.lower()
    is_fastq = any(s in suffix for s in [".fastq", ".fq", ".fastq.gz", ".fq.gz"])
    
    # Run paired-end header and line count synchronization checks first
    if is_fastq and in_r2:
        validate_paired_end_synchronicity(in_r1, in_r2)
    
    # Standard FASTA routing - directly bypasses FASTQ quality filters and runs fasta cleaning
    if not is_fastq:
        logger.info(f"Input '{in_r1.name}' detected as assembly FASTA file. Running FASTA preprocessing...")
        
        start_time = time.time()
        raw_records = []
        for rec in SeqIO.parse(str(in_r1), "fasta"):
            raw_records.append(rec)
            
        cleaned_records = []
        rejected_records = []
        
        seen_headers = {}
        seen_sequence_hashes = {}
        
        audit_data = []
        all_lengths = []
        all_gcs = []
        base_counts = Counter()
        
        duplicate_header_count = 0
        duplicate_sequence_count = 0
        
        min_len = config.preprocessing.min_length
        max_len = config.preprocessing.max_length
        max_ambig_pct = config.preprocessing.max_ambiguous_pct
        handle_headers = config.preprocessing.handle_duplicate_headers
        remove_seqs = config.preprocessing.remove_duplicate_sequences
        
        for rec in raw_records:
            orig_id = rec.id
            modified_id = orig_id
            is_header_dup = orig_id in seen_headers
            
            if is_header_dup:
                seen_headers[orig_id] += 1
                duplicate_header_count += 1
                logger.warning(f"Duplicate header detected: '{orig_id}'. Policy set to: '{handle_headers}'.")
                
                if handle_headers == "reject":
                    reason = f"Duplicate header (ID '{orig_id}' already processed)."
                    rejected_records.append((rec, reason))
                    audit_data.append({
                        "sequence_id": orig_id,
                        "length": len(rec.seq),
                        "gc_percent": None,
                        "ambiguous_count": None,
                        "ambiguous_percent": None,
                        "status": "REJECTED",
                        "rejection_reason": reason,
                        "duplicate_type": "Duplicate Header"
                    })
                    continue
                elif handle_headers == "rename":
                    suffix_idx = seen_headers[orig_id] - 1
                    modified_id = f"{orig_id}_dup{suffix_idx}"
                    logger.info(f"Renaming duplicate header '{orig_id}' to '{modified_id}'.")
                    rec.id = modified_id
            else:
                seen_headers[orig_id] = 1
                
            seq = str(rec.seq).upper()
            length = len(seq)
            ambig = sum(1 for base in seq if base in "N")
            amb_pct = (ambig / length) * 100.0 if length > 0 else 0.0
            
            # Simple FASTA Filters
            is_valid = True
            reason = ""
            if length < min_len:
                is_valid = False
                reason = f"Sequence length ({length} bp) below minimum threshold of {min_len} bp."
            elif length > max_len:
                is_valid = False
                reason = f"Sequence length ({length} bp) exceeds maximum threshold of {max_len} bp."
            elif amb_pct > max_ambig_pct:
                is_valid = False
                reason = f"Ambiguous bases ({amb_pct:.2f}%) exceed the configuration threshold of {max_ambig_pct}%."
                
            if not is_valid:
                logger.warning(f"Discarding FASTA record '{rec.id}': {reason}")
                rejected_records.append((rec, reason))
                audit_data.append({
                    "sequence_id": rec.id,
                    "length": length,
                    "gc_percent": None,
                    "ambiguous_count": None,
                    "ambiguous_percent": None,
                    "status": "REJECTED",
                    "rejection_reason": reason,
                    "duplicate_type": "None"
                })
                continue
                
            # Sequence duplicate check
            seq_hash = hashlib.md5(seq.encode("utf-8")).hexdigest()
            is_seq_dup = seq_hash in seen_sequence_hashes
            duplicate_type = "None"
            
            if is_seq_dup:
                original_match_id = seen_sequence_hashes[seq_hash]
                duplicate_sequence_count += 1
                duplicate_type = "Duplicate Sequence Content"
                logger.warning(f"Sequence content of '{rec.id}' is identical to '{original_match_id}'.")
                
                if remove_seqs:
                    reason = f"Duplicate sequence content (identical to '{original_match_id}')."
                    rejected_records.append((rec, reason))
                    audit_data.append({
                        "sequence_id": rec.id,
                        "length": length,
                        "gc_percent": None,
                        "ambiguous_count": None,
                        "ambiguous_percent": None,
                        "status": "REJECTED",
                        "rejection_reason": reason,
                        "duplicate_type": duplicate_type
                    })
                    continue
                else:
                    logger.info(f"Duplicate sequence '{rec.id}' kept and flagged.")
                    rec.description = f"{rec.description} [Duplicate Content Match: {original_match_id}]".strip()
            else:
                seen_sequence_hashes[seq_hash] = rec.id
                
            cleaned_records.append(rec)
            
            # Stats
            all_lengths.append(length)
            gc_count = sum(1 for base in seq if base in "GC")
            gc_pct = (gc_count / length) * 100.0 if length > 0 else 0.0
            all_gcs.append(gc_pct)
            base_counts.update(seq)
            
            audit_data.append({
                "sequence_id": rec.id,
                "length": length,
                "gc_percent": round(gc_pct, 2),
                "ambiguous_count": ambig,
                "ambiguous_percent": round(amb_pct, 2),
                "status": "ACCEPTED",
                "rejection_reason": "",
                "duplicate_type": duplicate_type
            })
            
        elapsed = time.time() - start_time
        
        # Write outputs
        if cleaned_records:
            SeqIO.write(cleaned_records, str(output_fasta), "fasta")
            logger.info(f"Saved {len(cleaned_records)} cleaned sequences to: {output_fasta}")
        else:
            logger.error("Zero FASTA sequences passed QC filters.")
            
        rejected_fasta_path = parent_dir / "rejected.fasta"
        if rejected_records:
            biopy_rejected = []
            for rec_obj, reason in rejected_records:
                desc = f"[REJECTED: {reason}]"
                biopy_rejected.append(SeqRecord(rec_obj.seq, id=rec_obj.id, description=desc))
            SeqIO.write(biopy_rejected, str(rejected_fasta_path), "fasta")
        else:
            with open(rejected_fasta_path, "w", encoding="utf-8") as f:
                f.write("")
                
        stats_csv_path = parent_dir / "qc_statistics.csv"
        pd.DataFrame(audit_data).to_csv(stats_csv_path, index=False)
        
        report_json_path = parent_dir / "preprocessing_report.json"
        
        summary = {
            "software": "FASTA-Preprocessor",
            "version": "1.0.0",
            "command_history": "python-fasta-cleaner",
            "runtime_seconds": elapsed,
            "counts": {
                "total_processed": len(raw_records),
                "total_kept": len(cleaned_records),
                "total_discarded": len(rejected_records),
                "duplicate_headers_encountered": duplicate_header_count,
                "duplicate_sequences_encountered": duplicate_sequence_count,
                "adapter_trimmed_reads": 0
            },
            "metrics": {
                "total_bases": sum(all_lengths),
                "min_length": min(all_lengths) if all_lengths else 0,
                "max_length": max(all_lengths) if all_lengths else 0,
                "mean_length": round(sum(all_lengths) / len(all_lengths), 2) if all_lengths else 0,
                "n50": calculate_n50(all_lengths),
                "mean_gc_percent": round(sum(all_gcs) / len(all_gcs), 2) if all_gcs else 0
            },
            "base_frequencies": dict(base_counts)
        }
        
        fasta_gc_bins = {i: 0 for i in range(101)}
        for g in all_gcs:
            bin_gc = int(round(g))
            if 0 <= bin_gc <= 100:
                fasta_gc_bins[bin_gc] += 1
                
        summary["qc_plots"] = {
            "per_base_quality": [],
            "gc_content_distribution": {str(k): v for k, v in fasta_gc_bins.items()},
            "sequence_duplication_levels": {
                "1": len(cleaned_records) - duplicate_sequence_count,
                "2": duplicate_sequence_count,
                "3-5": 0,
                "6-10": 0,
                "11-50": 0,
                ">50": 0
            },
            "read_retention_waterfall": {
                "raw_reads": len(raw_records),
                "adapter_trimmed": 0,
                "length_filtered": len(rejected_records),
                "quality_filtered": 0,
                "kept_reads": len(cleaned_records)
            }
        }
        
        summary["meta"] = {
            "pipeline": config.pipeline.name,
            "version": config.pipeline.version
        }
        
        # Save JSON preprocessing report
        with open(report_json_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=4)
            
        generate_premium_qc_html_report(summary, report_html_path)
        return summary
        
    # --- FASTQ quality processing ---
    logger.info("Initializing FASTQ Quality Control Preprocessing workflows...")
    
    summary = None
    cmd_history = ""
    
    # 1. Attempt Subprocess execution (fastp + FastQC)
    try:
        summary, cmd_history = run_fastp_subprocess(
            in_r1, out_r1, in_r2, out_r2, output_fasta, config, parent_dir
        )
        
        # FastQC post-QC audits
        try:
            run_fastqc_subprocess(out_r1, parent_dir)
            if in_r2:
                run_fastqc_subprocess(out_r2, parent_dir)
        except Exception as fe:
            logger.warning(f"FastQC report generation failed or skipped: {fe}")
            
    except (FileNotFoundError, ToolExecutionError) as e:
        logger.warning(f"Industrial tools subprocess execution bypassed/failed: {e}")
        
        # 2. Seamless High-fidelity Pure Python fallback
        summary, cmd_history = compile_pure_python_high_fidelity_preprocessing(
            in_r1, out_r1, in_r2, out_r2, output_fasta, config, parent_dir
        )
        
    # 3. Calculate high-fidelity QC metrics and binned distributions directly from sequence files
    try:
        adapter_count = summary.get("counts", {}).get("adapter_trimmed_reads", 0)
        rejection_reasons = summary.get("counts", {}).get("detailed_rejections", {})
        
        qc_plots = calculate_high_fidelity_qc_distributions(
            raw_r1=in_r1,
            raw_r2=in_r2,
            clean_r1=out_r1,
            clean_r2=out_r2,
            adapter_trimmed_count=adapter_count,
            rejection_reasons=rejection_reasons
        )
        summary["qc_plots"] = qc_plots
    except Exception as qce:
        logger.warning(f"Failed to calculate high-fidelity QC distributions: {qce}")
        summary["qc_plots"] = {
            "per_base_quality": [],
            "gc_content_distribution": {},
            "sequence_duplication_levels": {},
            "read_retention_waterfall": {}
        }

    summary["meta"] = {
        "pipeline": config.pipeline.name,
        "version": config.pipeline.version
    }
    
    # Generate MultiQC premium html dashboard
    generate_premium_qc_html_report(summary, report_html_path)
    
    return summary


class SequenceRecord:
    """
    Standardized internal representation of a biological sequence.
    Simplifies parsing of multiple input formats (FASTA, FASTQ, CSV, TXT).
    """
    def __init__(self, seq_id: str, sequence: str, description: str = ""):
        self.id = seq_id.strip()
        self.sequence = sequence.strip().upper()
        self.description = description.strip()
        
        if not self.id:
            raise SequenceValidationError("Sequence record contains an empty or missing identifier.")
        if not self.sequence:
            raise SequenceValidationError("Sequence record contains an empty sequence.")

    @property
    def length(self) -> int:
        return len(self.sequence)

    def to_seqrecord(self, custom_description: str = None) -> SeqRecord:
        """Converts to Biopython SeqRecord for downstream file writing."""
        desc = custom_description if custom_description is not None else self.description
        return SeqRecord(Seq(self.sequence), id=self.id, description=desc)


def parse_input_file(file_path: Path) -> List[SequenceRecord]:
    """
    Ingests sequence data from various standard formats: FASTA, FASTQ, CSV, TXT.
    Returns a list of standardized SequenceRecord objects.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {file_path}")

    suffix = file_path.suffix.lower()
    records = []

    logger.info(f"Parsing input file: {file_path} (Format: {suffix})")

    try:
        # Standard FASTA or FASTQ parsing using Biopython
        if suffix in [".fasta", ".fa", ".fna", ".faa"]:
            for rec in SeqIO.parse(str(file_path), "fasta"):
                records.append(SequenceRecord(rec.id, str(rec.seq), rec.description))
                
        elif suffix in [".fastq", ".fq"]:
            for rec in SeqIO.parse(str(file_path), "fastq"):
                records.append(SequenceRecord(rec.id, str(rec.seq), rec.description))
                
        # Tabular CSV/TSV format parsing
        elif suffix in [".csv", ".tsv", ".txt"]:
            sep = "," if suffix == ".csv" else ("\t" if suffix == ".tsv" else None)
            
            # Check if TXT file is raw line-by-line sequences or simple key-value pairs
            if suffix == ".txt" and sep is None:
                with open(file_path, "r", encoding="utf-8") as f:
                    content = f.read().strip()
                
                if content.startswith(">"):
                    # It's actually a FASTA file renamed to .txt
                    from io import StringIO
                    for rec in SeqIO.parse(StringIO(content), "fasta"):
                        records.append(SequenceRecord(rec.id, str(rec.seq), rec.description))
                else:
                    # Parse each line as an independent raw sequence
                    for idx, line in enumerate(content.splitlines()):
                        line = line.strip()
                        if line:
                            records.append(SequenceRecord(f"seq_{idx+1}", line, f"Raw line {idx+1}"))
            else:
                # CSV/TSV loading using pandas
                df = pd.read_csv(file_path, sep=sep)
                
                # Check for standard column headers
                id_col = next((col for col in df.columns if col.lower() in ["id", "identifier", "header", "name"]), None)
                seq_col = next((col for col in df.columns if col.lower() in ["seq", "sequence", "nucleotides", "bases"]), None)
                desc_col = next((col for col in df.columns if col.lower() in ["desc", "description", "comment"]), None)
                
                if not seq_col:
                    raise SequenceValidationError(
                        f"Tabular file requires a sequence column (e.g., 'sequence' or 'seq'). Columns found: {list(df.columns)}"
                    )
                
                for idx, row in df.iterrows():
                    seq_id = str(row[id_col]) if id_col else f"seq_{idx+1}"
                    seq_str = str(row[seq_col]) if not pd.isna(row[seq_col]) else ""
                    desc_str = str(row[desc_col]) if desc_col and not pd.isna(row[desc_col]) else ""
                    
                    records.append(SequenceRecord(seq_id, seq_str, desc_str))
        else:
            raise InvalidFileFormatError(f"Unsupported file format extension: '{suffix}'. Supported: FASTA, FASTQ, CSV, TSV, TXT.")
            
    except Exception as e:
        if isinstance(e, (PreprocessingError, FileNotFoundError)):
            raise e
        raise PreprocessingError(f"Failed parsing file: {e}")

    logger.info(f"Successfully loaded {len(records)} raw sequence records.")
    return records


def validate_and_clean_sequence(
    record: SequenceRecord,
    min_len: int,
    max_len: int,
    max_ambig_pct: float
) -> Tuple[bool, str, Dict[str, Any]]:
    """
    Validates a single sequence record.
    Returns a tuple of:
      - is_valid (bool)
      - rejection_reason (str)
      - record_stats (dict containing length, ambiguous bases, GC content, etc.)
    """
    seq = record.sequence
    length = len(seq)
    
    # 1. Empty Identifier or Empty Sequence Check
    if not record.id:
        return False, "Sequence record is missing an identifier (empty header).", {}
    if length == 0:
        return False, "Sequence is empty (0 bp).", {}

    # 2. Length Threshold Filter
    if length < min_len:
        return False, f"Sequence length ({length} bp) below minimum threshold of {min_len} bp.", {}
    if length > max_len:
        return False, f"Sequence length ({length} bp) exceeds maximum threshold of {max_len} bp.", {}

    # 3. Invalid Character Check (Strict IUPAC nucleotide code regex)
    invalid_chars = re.findall(r"[^ACGTURYWSKMBDHVN\-]", seq)
    if invalid_chars:
        unique_invalid = set(invalid_chars)
        return False, f"Sequence contains non-IUPAC invalid characters: {unique_invalid}", {}

    # 4. Ambiguous Nucleotide Calculation
    ambig_count = sum(1 for base in seq if base in "RYSWKMBDHVN")
    ambig_pct = (ambig_count / length) * 100.0 if length > 0 else 0.0
    
    if ambig_pct > max_ambig_pct:
        return False, f"Ambiguous bases ({ambig_pct:.2f}%) exceed the configuration threshold of {max_ambig_pct}%.", {}

    # 5. GC Content Calculation
    gc_count = sum(1 for base in seq if base in "GC")
    at_count = sum(1 for base in seq if base in "ATU")
    base_total = gc_count + at_count
    gc_pct = (gc_count / base_total) * 100.0 if base_total > 0 else 0.0

    stats = {
        "length": length,
        "gc_percent": round(gc_pct, 2),
        "ambiguous_count": ambig_count,
        "ambiguous_percent": round(ambig_pct, 2),
        "md5": hashlib.md5(seq.encode("utf-8")).hexdigest()
    }
    
    return True, "", stats

