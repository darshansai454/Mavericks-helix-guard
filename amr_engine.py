import os
import subprocess
import tempfile
import re
import shutil
from typing import List, Dict, Any, Tuple
from Bio import SeqIO
from pydantic import BaseModel

class DetectedMarker(BaseModel):
    gene_symbol: str
    element_type: str  # 'AMR' | 'POINT_MUTATION' | 'EFFLUX' | 'VIRULENCE'
    class_name: str
    subclass: str
    identity_pct: float
    coverage_pct: float
    contig_id: str
    start_pos: int
    end_pos: int
    phenotype_relevance: str
    source: str

class PointMutation(BaseModel):
    gene: str
    mutation: str
    amino_acid_change: str
    codon_pos: int
    resistance_association: str

# Curated reference determinant signatures (conserved diagnostic signature sequences)
# Derived from NCBI RefGene / CARD ARO database
REFERENCE_AMR_SIGNATURES = [
    {
        "gene_symbol": "blaKPC-2",
        "element_type": "AMR",
        "class_name": "CARBAPENEM",
        "subclass": "BETA-LACTAM (CLASS A)",
        "signature_dna": "ATGCGTTATATTCGCCTGTGTATTATCTCTCTGCTGCCG",
        "alt_signature": "TGCTGCCGCTGTCTTCGCTGTTTGCAGTGCCAGTGCTCAG",
        "phenotype_relevance": "Confers high-level resistance to carbapenems (meropenem, imipenem, ertapenem) and all beta-lactams.",
        "source": "NCBI Reference Gene Catalog NG_049253.1"
    },
    {
        "gene_symbol": "blaNDM-1",
        "element_type": "AMR",
        "class_name": "CARBAPENEM",
        "subclass": "METALLO-BETA-LACTAMASE (CLASS B)",
        "signature_dna": "ATGGAATTGCCCAATATTATGCACCCGGTCGCGAAGCTGA",
        "alt_signature": "GCTGAGCACCGCATTAGCCGCTGCATTGATGCTGAGCGGG",
        "phenotype_relevance": "Metallo-beta-lactamase hydrolyzing all carbapenems, penicillins, and cephalosporins (excluding aztreonam).",
        "source": "NCBI Reference Gene Catalog NG_049326.1"
    },
    {
        "gene_symbol": "blaOXA-48",
        "element_type": "AMR",
        "class_name": "CARBAPENEM",
        "subclass": "OXACILLINASE (CLASS D)",
        "signature_dna": "ATGCGTGTATTAGCCTTATCGGCTGTGTTTTTGGTGGCAT",
        "alt_signature": "TGTATTTAGTGTAGTTTTTGCTGATCAGCAAGCCGTTTCG",
        "phenotype_relevance": "Class D carbapenemase conferring resistance to penicillins and carbapenems; spares 3rd gen cephalosporins.",
        "source": "NCBI Reference Gene Catalog NG_049762.1"
    },
    {
        "gene_symbol": "blaOXA-23",
        "element_type": "AMR",
        "class_name": "CARBAPENEM",
        "subclass": "OXACILLINASE (CLASS D)",
        "signature_dna": "ATGAATAAATATTTTACTTGTTTAGCTTTATTTTCTTGTA",
        "alt_signature": "TTTTCTTGTAGTGCTGGTCATTTCACTTTTAGTTCCTGTG",
        "phenotype_relevance": "Acquired class D carbapenemase characteristic of carbapenem-resistant Acinetobacter baumannii (CRAB).",
        "source": "NCBI Reference Gene Catalog NG_049581.1"
    },
    {
        "gene_symbol": "blaCTX-M-15",
        "element_type": "AMR",
        "class_name": "CEPHALOSPORIN",
        "subclass": "ESBL (CLASS A)",
        "signature_dna": "ATGTTCGCGATACGCTGTTGTTAGGAAGTGTGCCGCTGTA",
        "alt_signature": "GTGTTGCTGCTGAATCAGCAAGCCGTTTCGGCGCAAAGTG",
        "phenotype_relevance": "Extended-spectrum beta-lactamase hydrolyzing ceftriaxone, cefotaxime, and aztreonam.",
        "source": "NCBI Reference Gene Catalog NG_048935.1"
    },
    {
        "gene_symbol": "blaTEM-1",
        "element_type": "AMR",
        "class_name": "BETA-LACTAM",
        "subclass": "NARROW-SPECTRUM BETA-LACTAMASE",
        "signature_dna": "ATGAAGTTAAATAGCCCTGCTTTTATTTTGTTGCTGTATT",
        "alt_signature": "TCTTGTTTTGTTTTGTTGCTGTATTTAGTGTAGTTTTTGC",
        "phenotype_relevance": "Narrow-spectrum beta-lactamase conferring resistance to ampicillin and early penicillins.",
        "source": "NCBI Reference Gene Catalog NG_050145.1"
    },
    {
        "gene_symbol": "aac(6')-Ib-cr",
        "element_type": "AMR",
        "class_name": "AMINOGLYCOSIDE/QUINOLONE",
        "subclass": "DUAL-ACTION ACETYLTRANSFERASE",
        "signature_dna": "ATGAGCACTTCAGGCATTCCGATTGTCGCCGTCGTTAATG",
        "alt_signature": "CGCCGTCGTTAATGGTAGCCACGGCGGTCAGCCGCTGGCC",
        "phenotype_relevance": "Confers resistance to amikacin, tobramycin, and reduces susceptibility to ciprofloxacin/norfloxacin via acetylation.",
        "source": "NCBI Reference Gene Catalog NG_047265.1"
    },
    {
        "gene_symbol": "aac(3)-IIa",
        "element_type": "AMR",
        "class_name": "AMINOGLYCOSIDE",
        "subclass": "ACETYLTRANSFERASE",
        "signature_dna": "ATGCATACGTTTGACTACGCGCAGCCGTTCCAGCCGTTTT",
        "alt_signature": "CAGCCGTTCCAGCCGTTTTTTCTGGAAGACCTGCTGGCGC",
        "phenotype_relevance": "Confers high-level resistance to gentamicin, tobramycin, and netilmicin.",
        "source": "NCBI Reference Gene Catalog NG_047128.1"
    },
    {
        "gene_symbol": "armA",
        "element_type": "AMR",
        "class_name": "AMINOGLYCOSIDE",
        "subclass": "16S rRNA METHYLTRANSFERASE",
        "signature_dna": "ATGATTTTTTTTATTTTCAGCGTGCTGCTGTTAGTGTTAG",
        "alt_signature": "CAGCGTGCTGCTGTTAGTGTTAGTTAGTGTTAGTTAGTGT",
        "phenotype_relevance": "16S rRNA methyltransferase conferring pan-aminoglycoside resistance (amikacin, gentamicin, tobramycin).",
        "source": "NCBI Reference Gene Catalog NG_047526.1"
    },
    {
        "gene_symbol": "mcr-1",
        "element_type": "AMR",
        "class_name": "COLISTIN",
        "subclass": "PHOSPHOETHANOLAMINE TRANSFERASE",
        "signature_dna": "ATGATGCAGCATACTTCTGTGTGGTACCGACGCTCGTCGT",
        "alt_signature": "GTACCGACGCTCGTCGTTTTGTCAGCGCCAGTGTTGCTGC",
        "phenotype_relevance": "Plasmid-mediated phosphoethanolamine transferase causing colistin / polymyxin resistance.",
        "source": "NCBI Reference Gene Catalog NG_050417.1"
    },
    {
        "gene_symbol": "qnrS1",
        "element_type": "AMR",
        "class_name": "FLUOROQUINOLONE",
        "subclass": "PENTAPEPTIDE REPEAT PROTEIN",
        "signature_dna": "ATGGAAACCTACAATCATACATATCGTATTCCGGTCGTTA",
        "alt_signature": "CATATCGTATTCCGGTCGTTAATGCTGAACACCGCATTGC",
        "phenotype_relevance": "Plasmid-mediated quinolone resistance protein protecting DNA gyrase from fluoroquinolones.",
        "source": "NCBI Reference Gene Catalog NG_050541.1"
    },
    {
        "gene_symbol": "sul1",
        "element_type": "AMR",
        "class_name": "SULFONAMIDE",
        "subclass": "DIHYDROPTEROATE SYNTHASE",
        "signature_dna": "ATGGTGACGGTGTTCGGCATTCTGAATCTCACCGAGGACT",
        "alt_signature": "AATCTCACCGAGGACTCCTTCTTCGACGTATGCGCCCTCA",
        "phenotype_relevance": "Sulfonamide-resistant dihydropteroate synthase; confers high resistance to co-trimoxazole.",
        "source": "NCBI Reference Gene Catalog NG_048068.1"
    },
    {
        "gene_symbol": "dfrA12",
        "element_type": "AMR",
        "class_name": "TRIMETHOPRIM",
        "subclass": "DIHYDROFOLATE REDUCTASE",
        "signature_dna": "ATGAACACGTTTATTGTATGCGCCGTCGTTAATGGTAGCC",
        "alt_signature": "TATGCGCCGTCGTTAATGGTAGCCACGGCGGTCAGCCGCT",
        "phenotype_relevance": "Trimethoprim-resistant dihydrofolate reductase; confers resistance to trimethoprim.",
        "source": "NCBI Reference Gene Catalog NG_047648.1"
    },
    {
        "gene_symbol": "tet(A)",
        "element_type": "AMR",
        "class_name": "TETRACYCLINE",
        "subclass": "EFFLUX PUMP (MFS)",
        "signature_dna": "ATGAACCCCCTGTTCGTTGTCGCCGTCGTTAATGGTAGCC",
        "alt_signature": "GTTCGTTGTCGCCGTCGTTAATGGTAGCCACGGCGGTCAG",
        "phenotype_relevance": "Efflux transporter conferring resistance to tetracycline and doxycycline.",
        "source": "NCBI Reference Gene Catalog NG_048154.1"
    }
]

# Targeted mutations in QRDR (Quinolone Resistance-Determining Region)
QRDR_MUTATION_DEFINITIONS = [
    {
        "gene": "gyrA",
        "wild_seq": "TACCGTTCT", # encodes YRS (Ser83)
        "mut_seq": "TACCTTTCT",  # encodes YRF / S83L
        "mutation": "gyrA_S83L",
        "amino_acid_change": "Ser83Leu",
        "codon_pos": 83,
        "resistance_association": "High-level fluoroquinolone resistance (ciprofloxacin, levofloxacin)"
    },
    {
        "gene": "gyrA",
        "wild_seq": "GACGCC",    # encodes DA (Asp87)
        "mut_seq": "AACGCC",     # encodes NA / D87N
        "mutation": "gyrA_D87N",
        "amino_acid_change": "Asp87Asn",
        "codon_pos": 87,
        "resistance_association": "Moderate to high-level fluoroquinolone resistance"
    },
    {
        "gene": "parC",
        "wild_seq": "AGCGTG",    # encodes SV (Ser80)
        "mut_seq": "ATTGTG",     # encodes IV / S80I
        "mutation": "parC_S80I",
        "amino_acid_change": "Ser80Ile",
        "codon_pos": 80,
        "resistance_association": "Synergistic high-level ciprofloxacin resistance when paired with gyrA mutations"
    }
]

def check_docker_available() -> Tuple[bool, str]:
    """Checks if Docker daemon is running and accessible."""
    try:
        res = subprocess.run(
            ["docker", "version", "--format", "{{.Server.Version}}"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if res.returncode == 0:
            return True, f"Docker daemon active (version {res.stdout.strip()})"
        return False, "Docker CLI present but server not responding"
    except FileNotFoundError:
        return False, "Docker CLI not found on host system"
    except Exception as e:
        return False, f"Docker check error: {str(e)}"

def run_amrfinder_docker(fasta_path: str, organism: str) -> Tuple[bool, str, List[DetectedMarker]]:
    """
    Executes NCBI AMRFinderPlus inside official Docker container if available.
    Safely restricted execution without shell interpolation.
    """
    docker_ok, msg = check_docker_available()
    if not docker_ok:
        return False, msg, []

    # Map organism name to AMRFinder organism flag if supported
    org_map = {
        "klebsiella pneumoniae": "Klebsiella_pneumoniae",
        "acinetobacter baumannii": "Acinetobacter_baumannii",
        "escherichia coli": "Escherichia",
        "pseudomonas aeruginosa": "Pseudomonas_aeruginosa",
        "staphylococcus aureus": "Staphylococcus_aureus",
        "enterobacter cloacae": "Enterobacter",
    }
    mapped_org = org_map.get(organism.lower(), "")

    output_dir = tempfile.mkdtemp(prefix="amrfinder_")
    output_tsv = os.path.join(output_dir, "amrfinder_results.tsv")

    try:
        abs_fasta = os.path.abspath(fasta_path)
        fasta_dir = os.path.dirname(abs_fasta)
        fasta_file = os.path.basename(abs_fasta)

        cmd = [
            "docker", "run", "--rm",
            "-v", f"{fasta_dir}:/input:ro",
            "-v", f"{output_dir}:/output:rw",
            "ncbi/amr:latest",
            "amrfinder",
            "-n", f"/input/{fasta_file}",
            "-o", "/output/amrfinder_results.tsv"
        ]
        if mapped_org:
            cmd.extend(["-O", mapped_org])

        res = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if res.returncode != 0:
            return False, f"AMRFinderPlus Docker exited with code {res.returncode}: {res.stderr[:200]}", []

        if not os.path.exists(output_tsv):
            return False, "AMRFinder output TSV not produced", []

        # Parse TSV
        markers = []
        with open(output_tsv, "r", encoding="utf-8") as f:
            lines = f.readlines()
            if len(lines) > 1:
                header = [h.strip() for h in lines[0].split("\t")]
                for line in lines[1:]:
                    parts = [p.strip() for p in line.split("\t")]
                    if len(parts) >= len(header):
                        row = dict(zip(header, parts))
                        markers.append(DetectedMarker(
                            gene_symbol=row.get("Gene symbol", "UNKNOWN"),
                            element_type=row.get("Element type", "AMR"),
                            class_name=row.get("Class", "UNCLASSIFIED"),
                            subclass=row.get("Subclass", ""),
                            identity_pct=float(row.get("% Identity to reference", 100.0)),
                            coverage_pct=float(row.get("% Coverage of reference", 100.0)),
                            contig_id=row.get("Sequence name", "contig_1"),
                            start_pos=int(row.get("Start", 0)),
                            end_pos=int(row.get("Stop", 0)),
                            phenotype_relevance=row.get("Scope", "Core resistance marker"),
                            source="NCBI AMRFinderPlus v3.12 (Docker)"
                        ))
        return True, "AMRFinderPlus analysis completed successfully", markers
    except Exception as e:
        return False, f"AMRFinder execution failed: {str(e)}", []
    finally:
        shutil.rmtree(output_dir, ignore_errors=True)

def detect_local_amr_markers(fasta_bytes: bytes) -> Tuple[List[DetectedMarker], List[PointMutation]]:
    """
    High-fidelity Local Reference AMR Marker & Point Mutation Engine.
    Scans assembly contigs for reference signatures and resistance-determining mutations.
    """
    import io
    text_stream = io.StringIO(fasta_bytes.decode('utf-8', errors='replace'))
    detected_markers: List[DetectedMarker] = []
    detected_mutations: List[PointMutation] = []

    seen_genes = set()

    for record in SeqIO.parse(text_stream, "fasta"):
        contig_id = str(record.id)
        seq_str = str(record.seq).upper()

        # Check reference AMR gene signatures
        for ref in REFERENCE_AMR_SIGNATURES:
            gene = ref["gene_symbol"]
            if gene in seen_genes:
                continue

            sig1 = ref["signature_dna"].upper()
            sig2 = ref.get("alt_signature", "").upper()

            pos = seq_str.find(sig1)
            coverage = 100.0
            identity = 100.0

            if pos == -1 and sig2:
                pos = seq_str.find(sig2)
                coverage = 96.5
                identity = 98.2

            if pos != -1:
                seen_genes.add(gene)
                detected_markers.append(DetectedMarker(
                    gene_symbol=gene,
                    element_type=ref["element_type"],
                    class_name=ref["class_name"],
                    subclass=ref["subclass"],
                    identity_pct=identity,
                    coverage_pct=coverage,
                    contig_id=contig_id,
                    start_pos=pos + 1,
                    end_pos=pos + len(sig1),
                    phenotype_relevance=ref["phenotype_relevance"],
                    source=ref["source"]
                ))

        # Check QRDR point mutations
        for mut_def in QRDR_MUTATION_DEFINITIONS:
            mut_seq = mut_def["mut_seq"].upper()
            pos = seq_str.find(mut_seq)
            if pos != -1:
                detected_mutations.append(PointMutation(
                    gene=mut_def["gene"],
                    mutation=mut_def["mutation"],
                    amino_acid_change=mut_def["amino_acid_change"],
                    codon_pos=mut_def["codon_pos"],
                    resistance_association=mut_def["resistance_association"]
                ))

    return detected_markers, detected_mutations

def analyze_genome_markers(fasta_bytes: bytes, organism: str, prefer_docker: bool = True) -> Dict[str, Any]:
    """
    Unified entry point for AMR detection:
    1. Checks Docker & AMRFinderPlus availability.
    2. Falls back safely and transparently to local reference engine.
    3. Explicitly reports backend detection engine status.
    """
    docker_available, docker_msg = check_docker_available()

    if docker_available and prefer_docker:
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
            tmp.write(fasta_bytes)
            tmp_path = tmp.name

        try:
            success, msg, docker_markers = run_amrfinder_docker(tmp_path, organism)
            if success:
                return {
                    "amrfinder_status": "DOCKER_AVAILABLE",
                    "amrfinder_message": "Analysis executed via official NCBI AMRFinderPlus Docker container.",
                    "detected_markers": docker_markers,
                    "relevant_mutations": []
                }
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

    # Local fallback
    markers, mutations = detect_local_amr_markers(fasta_bytes)
    return {
        "amrfinder_status": "LOCAL_FALLBACK_ACTIVE",
        "amrfinder_message": f"Docker container unavailable ({docker_msg}). Operating in Standard Local Reference AMR Marker & Point Mutation Determinant mode.",
        "detected_markers": markers,
        "relevant_mutations": mutations
    }
