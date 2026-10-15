import hashlib
import io
import os
import re
from typing import Dict, Any, Tuple, List
from Bio import SeqIO
from pydantic import BaseModel

class GenomeQuality(BaseModel):
    total_length_bp: int
    contig_count: int
    n50_bp: int
    gc_content_pct: float
    ambiguous_bases: int
    passed_quality: bool = True
    quality_notes: List[str] = []

VALID_EXTENSIONS = {'.fasta', '.fa', '.fna'}
MAX_FILE_SIZE_BYTES = 50 * 1024 * 1024  # 50MB limit

class GenomeValidationError(Exception):
    pass

def calculate_sha256(content_bytes: bytes) -> str:
    hasher = hashlib.sha256()
    hasher.update(content_bytes)
    return hasher.hexdigest()

def calculate_n50(contig_lengths: List[int]) -> int:
    if not contig_lengths:
        return 0
    sorted_lengths = sorted(contig_lengths, reverse=True)
    total_length = sum(sorted_lengths)
    cumulative = 0
    target = total_length / 2.0
    for length in sorted_lengths:
        cumulative += length
        if cumulative >= target:
            return length
    return sorted_lengths[-1]

def validate_fasta_content(filename: str, content_bytes: bytes) -> Tuple[str, GenomeQuality, List[Dict[str, Any]]]:
    """
    Validates FASTA file content with Biopython.
    Returns:
        (file_sha256, GenomeQuality, parsed_records_summary)
    Raises:
        GenomeValidationError if malformed or invalid
    """
    # 1. File extension validation
    _, ext = os.path.splitext(filename.lower())
    if ext not in VALID_EXTENSIONS:
        raise GenomeValidationError(
            f"Unsupported file extension '{ext}'. Only .fasta, .fa, and .fna are permitted."
        )

    # 2. File size validation
    if len(content_bytes) == 0:
        raise GenomeValidationError("File is completely empty (0 bytes). Upload a valid FASTA genome file.")
    if len(content_bytes) > MAX_FILE_SIZE_BYTES:
        raise GenomeValidationError(
            f"File size ({len(content_bytes) / (1024*1024):.1f} MB) exceeds maximum allowed limit (50 MB)."
        )

    # 3. Check for binary characters / executable signatures
    if b'\x00' in content_bytes[:1024] or content_bytes.startswith(b'MZ') or content_bytes.startswith(b'\x7fELF'):
        raise GenomeValidationError("Binary or executable file detected. Only standard plaintext FASTA nucleotide files are allowed.")

    file_sha256 = calculate_sha256(content_bytes)

    # 4. Biopython parsing
    try:
        text_stream = io.StringIO(content_bytes.decode('utf-8', errors='replace'))
    except Exception as e:
        raise GenomeValidationError(f"Could not decode file content as text: {str(e)}")

    records = []
    try:
        for record in SeqIO.parse(text_stream, "fasta"):
            records.append(record)
    except Exception as e:
        raise GenomeValidationError(f"FASTA syntax error during parsing: {str(e)}")

    if not records:
        raise GenomeValidationError("No valid FASTA records found. Ensure headers start with '>' followed by valid sequence lines.")

    contig_lengths = []
    total_bases = 0
    gc_bases = 0
    ambiguous_bases = 0
    records_summary = []

    # Valid IUPAC nucleotide characters: A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N, -
    valid_nucleotides = set("ACGTURYWSKMBVDHN-acgturywskmbvdhn")

    for idx, rec in enumerate(records):
        seq_str = str(rec.seq)
        seq_len = len(seq_str)
        if seq_len == 0:
            raise GenomeValidationError(f"Contig '{rec.id}' has 0 sequence length.")
        
        # Check invalid non-nucleotide characters
        sample_chars = set(seq_str[:5000])
        invalid = sample_chars - valid_nucleotides
        if invalid:
            raise GenomeValidationError(f"Contig '{rec.id}' contains non-nucleotide characters: {list(invalid)[:5]}")

        contig_lengths.append(seq_len)
        total_bases += seq_len
        
        # Count GC and Ns
        seq_upper = seq_str.upper()
        gc_bases += seq_upper.count('G') + seq_upper.count('C')
        ambiguous_bases += seq_upper.count('N')

        records_summary.append({
            "id": rec.id,
            "description": rec.description[:100],
            "length": seq_len,
            "gc_pct": round(((seq_upper.count('G') + seq_upper.count('C')) / seq_len) * 100, 2)
        })

    n50 = calculate_n50(contig_lengths)
    gc_pct = round((gc_bases / total_bases) * 100, 2) if total_bases > 0 else 0.0

    # Quality notes
    quality_notes = []
    passed = True

    if total_bases < 500_000:
        quality_notes.append(f"Assembly size ({total_bases:,} bp) is smaller than typical bacterial genomes (>500 kbp). Likely partial assembly.")
        passed = False
    elif total_bases > 12_000_000:
        quality_notes.append(f"Assembly size ({total_bases:,} bp) is unusually large for single bacterial isolate (>12 Mbp). Potential contamination.")
        passed = False
    else:
        quality_notes.append(f"Genome length ({total_bases:,} bp) is consistent with typical bacterial assemblies.")

    if len(records) > 1000:
        quality_notes.append(f"High contig fragmentation ({len(records)} contigs). N50 is {n50:,} bp.")
    else:
        quality_notes.append(f"Good assembly contiguity ({len(records)} contigs, N50 {n50:,} bp).")

    ambiguity_pct = (ambiguous_bases / total_bases) * 100 if total_bases > 0 else 0
    if ambiguity_pct > 5.0:
        quality_notes.append(f"Elevated ambiguous 'N' base content ({ambiguity_pct:.2f}%).")
        passed = False

    quality = GenomeQuality(
        total_length_bp=total_bases,
        contig_count=len(records),
        n50_bp=n50,
        gc_content_pct=gc_pct,
        ambiguous_bases=ambiguous_bases,
        passed_quality=passed,
        quality_notes=quality_notes
    )

    return file_sha256, quality, records_summary
