import os
import json
import uuid
import io
from datetime import datetime, timezone
from typing import List, Optional
from fastapi import APIRouter, File, UploadFile, Form, HTTPException, Query, Response
from fastapi.responses import StreamingResponse

from ..database.db import (
    save_analysis,
    get_analysis_by_id,
    list_analyses,
    delete_analysis_by_id
)
from ..genome.validator import validate_fasta_content, GenomeValidationError
from ..genome.amr_engine import analyze_genome_markers, check_docker_available
from ..rules.intrinsic_rules import get_all_rules, check_intrinsic_override
from ..ml.prediction_engine import (
    predict_antibiotic_susceptibility,
    SUPPORTED_ANTIBIOTICS
)
from ..ml.metrics import get_model_performance_summary
from ..reports.pdf_generator import generate_pdf_report
from ..schemas.api_models import (
    HealthResponse,
    OrganismResponse,
    AntibioticResponse,
    ValidationResponse,
    AnalysisInitiateResponse,
    RuleResponse
)

router = APIRouter()

def get_sample_genome_path(sample_id: str) -> Optional[str]:
    candidates = [
        os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "sample_data", "genomes")),
        os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "sample_data", "genomes")),
        os.path.abspath(os.path.join(os.getcwd(), "sample_data", "genomes"))
    ]
    for d in candidates:
        p = os.path.join(d, f"{sample_id}.fasta")
        if os.path.exists(p):
            return p
    return None

SUPPORTED_ORGANISMS = [
    {
        "id": "klebsiella_pneumoniae",
        "scientific_name": "Klebsiella pneumoniae",
        "common_name": "Klebsiella pneumoniae (ESKAPE)",
        "gram_type": "NEGATIVE",
        "estimated_genome_mbp": "5.3 - 5.7 Mbp",
        "typical_gc_pct": "57.0%",
        "priority_pathogen": "WHO Critical Priority",
        "supported_antibiotics_count": 18
    },
    {
        "id": "acinetobacter_baumannii",
        "scientific_name": "Acinetobacter baumannii",
        "common_name": "Acinetobacter baumannii (CRAB)",
        "gram_type": "NEGATIVE",
        "estimated_genome_mbp": "3.8 - 4.2 Mbp",
        "typical_gc_pct": "39.0%",
        "priority_pathogen": "WHO Critical Priority",
        "supported_antibiotics_count": 18
    },
    {
        "id": "escherichia_coli",
        "scientific_name": "Escherichia coli",
        "common_name": "Escherichia coli",
        "gram_type": "NEGATIVE",
        "estimated_genome_mbp": "4.6 - 5.5 Mbp",
        "typical_gc_pct": "50.8%",
        "priority_pathogen": "WHO Critical Priority",
        "supported_antibiotics_count": 18
    },
    {
        "id": "pseudomonas_aeruginosa",
        "scientific_name": "Pseudomonas aeruginosa",
        "common_name": "Pseudomonas aeruginosa",
        "gram_type": "NEGATIVE",
        "estimated_genome_mbp": "6.2 - 6.8 Mbp",
        "typical_gc_pct": "66.5%",
        "priority_pathogen": "WHO Critical Priority",
        "supported_antibiotics_count": 18
    },
    {
        "id": "staphylococcus_aureus",
        "scientific_name": "Staphylococcus aureus",
        "common_name": "Staphylococcus aureus (MRSA)",
        "gram_type": "POSITIVE",
        "estimated_genome_mbp": "2.8 - 3.0 Mbp",
        "typical_gc_pct": "32.8%",
        "priority_pathogen": "WHO High Priority",
        "supported_antibiotics_count": 14
    },
    {
        "id": "enterobacter_cloacae",
        "scientific_name": "Enterobacter cloacae",
        "common_name": "Enterobacter cloacae complex",
        "gram_type": "NEGATIVE",
        "estimated_genome_mbp": "4.8 - 5.3 Mbp",
        "typical_gc_pct": "55.0%",
        "priority_pathogen": "WHO High Priority",
        "supported_antibiotics_count": 18
    }
]

@router.get("/health", response_model=HealthResponse)
def health_check():
    docker_ok, docker_msg = check_docker_available()
    return HealthResponse(
        status="UP",
        service="Helix Guard Prediction Engine",
        version="2.4.0",
        docker_status=docker_ok,
        docker_message=docker_msg,
        amrfinder_engine="NCBI AMRFinderPlus Docker" if docker_ok else "Standard Local Determinant & Mutation Engine",
        database_connected=True
    )

@router.get("/organisms", response_model=List[OrganismResponse])
def get_organisms():
    return [OrganismResponse(**org) for org in SUPPORTED_ORGANISMS]

ANTIBIOTIC_CLINICAL_PROFILES = {
    "ciprofloxacin": {
        "description": "Broad-spectrum oral/IV fluoroquinolone widely used in complicated urinary, enteric, and respiratory infections.",
        "target_pathogens": ["Escherichia coli", "Pseudomonas aeruginosa", "Salmonella enterica"],
        "resistance_mechanisms": ["gyrA/parC QRDR target mutations", "qnr plasmid protection", "aac(6')-Ib-cr acetylation"],
        "clinical_indications": "Pyelonephritis, complicated intra-abdominal infections, typhoid fever, anthrax.",
        "typical_administration": "Oral tablets / IV infusion",
        "molecular_formula": "C17H18FN3O3",
        "risk_level": "High Clinical Value",
        "visual_color": "#0284c7"
    },
    "meropenem": {
        "description": "Ultra-broad spectrum intravenous carbapenem reserved for severe nosocomial multidrug-resistant Gram-negative infections.",
        "target_pathogens": ["Klebsiella pneumoniae", "Acinetobacter baumannii", "Pseudomonas aeruginosa"],
        "resistance_mechanisms": ["KPC/NDM/OXA-48 carbapenemases", "OmpK35/36 porin loss", "MexAB-OprM efflux"],
        "clinical_indications": "MDR sepsis, ventilator-associated pneumonia, intra-abdominal sepsis, febrile neutropenia.",
        "typical_administration": "IV infusion (1-2g q8h)",
        "molecular_formula": "C17H25N3O5S",
        "risk_level": "WHO Watch / Reserve",
        "visual_color": "#0d9488"
    },
    "amikacin": {
        "description": "High-potency aminoglycoside resistant to most modifying enzymes; vital for pseudomonal and enterobacterial sepsis.",
        "target_pathogens": ["Pseudomonas aeruginosa", "Enterobacter cloacae", "Mycobacterium abscessus"],
        "resistance_mechanisms": ["armA/rmtB 16S rRNA methyltransferases", "aac(6')-Ib acetyltransferases"],
        "clinical_indications": "Severe hospital-acquired sepsis, resistant urinary infections, atypical mycobacteria.",
        "typical_administration": "IV / IM once-daily dosing with TDM",
        "molecular_formula": "C22H43N5O13",
        "risk_level": "WHO Watch Group",
        "visual_color": "#2563eb"
    },
    "ceftriaxone": {
        "description": "Workhorse 3rd-generation cephalosporin with excellent CSF penetration for pneumonia, pyelonephritis, and meningitis.",
        "target_pathogens": ["Streptococcus pneumoniae", "Neisseria meningitidis", "Escherichia coli"],
        "resistance_mechanisms": ["CTX-M/SHV/TEM extended-spectrum beta-lactamases (ESBL)", "PBP alterations"],
        "clinical_indications": "Community-acquired pneumonia, acute bacterial meningitis, Lyme disease, gonococcal infection.",
        "typical_administration": "IV / IM (1-2g daily)",
        "molecular_formula": "C18H18N8O7S3",
        "risk_level": "WHO Watch Group",
        "visual_color": "#7c3aed"
    },
    "colistin": {
        "description": "Last-line cationic lipopeptide polymyxin that binds and disrupts bacterial outer lipopolysaccharide membranes.",
        "target_pathogens": ["Carbapenem-resistant Acinetobacter baumannii (CRAB)", "MDR P. aeruginosa", "CRE Klebsiella"],
        "resistance_mechanisms": ["mcr-1 to mcr-9 phosphoethanolamine transferases", "pmrA/pmrB, phoP/phoQ two-component mutations"],
        "clinical_indications": "Life-threatening pan-drug-resistant Gram-negative infections when all beta-lactams have failed.",
        "typical_administration": "IV colistimethate sodium / aerosolized",
        "molecular_formula": "C52H98N16O13",
        "risk_level": "WHO Reserve Antibiotic",
        "visual_color": "#e11d48"
    },
    "gentamicin": {
        "description": "Systemic bactericidal aminoglycoside commonly paired with cell wall agents for synergy in endocarditis.",
        "target_pathogens": ["Enterobacteriaceae", "Staphylococcus aureus", "Enterococcus faecalis (synergy)"],
        "resistance_mechanisms": ["aac(3)-IIa acetyltransferase", "ant(2'')-Ia nucleotidyltransferase", "16S rRNA methylation"],
        "clinical_indications": "Bacterial endocarditis synergy, Gram-negative bacteremia, pelvic inflammatory disease.",
        "typical_administration": "IV infusion with therapeutic drug monitoring",
        "molecular_formula": "C21H43N5O7",
        "risk_level": "WHO Access / Watch",
        "visual_color": "#0891b2"
    },
    "piperacillin_tazobactam": {
        "description": "Broad-spectrum antipseudomonal penicillin paired with a beta-lactamase inhibitor (Zosyn) for empirical hospital therapy.",
        "target_pathogens": ["Pseudomonas aeruginosa", "Bacteroides fragilis", "Klebsiella pneumoniae"],
        "resistance_mechanisms": ["Class B metallo-beta-lactamases", "Class A carbapenemases (KPC)", "AmpC hyperproduction"],
        "clinical_indications": "Hospital-acquired pneumonia, complicated intra-abdominal sepsis, diabetic foot infections.",
        "typical_administration": "Extended IV infusion (3.375g - 4.5g q6h-q8h)",
        "molecular_formula": "C23H27N5O7S / C10H12N4O5S",
        "risk_level": "WHO Watch Group",
        "visual_color": "#059669"
    },
    "ampicillin": {
        "description": "Early aminopenicillin active against Listeria, Enterococci, and susceptible wild-type Gram-negatives.",
        "target_pathogens": ["Listeria monocytogenes", "Enterococcus faecalis", "Streptococcus agalactiae"],
        "resistance_mechanisms": ["blaTEM-1 plasmid beta-lactamase", "blaSHV chromosomal beta-lactamase in Klebsiella"],
        "clinical_indications": "Neonatal meningitis (with cefotaxime/gentamicin), enterococcal endocarditis, listeriosis.",
        "typical_administration": "Oral / IV (500mg - 2g q4-6h)",
        "molecular_formula": "C16H19N3O4S",
        "risk_level": "WHO Access Group",
        "visual_color": "#d97706"
    },
    "trimethoprim_sulfamethoxazole": {
        "description": "Bactericidal sequential inhibitor of folic acid synthesis (Bactrim / Septra) with broad antimicrobial spectrum.",
        "target_pathogens": ["Stenotrophomonas maltophilia", "Pneumocystis jirovecii", "Nocardia spp."],
        "resistance_mechanisms": ["sul1, sul2, sul3 dihydropteroate mutations", "dfrA12, dfrA17 dihydrofolate reductase alleles"],
        "clinical_indications": "Pneumocystis pneumonia prophylaxis, Stenotrophomonas bacteremia, acute uncomplicated cystitis.",
        "typical_administration": "Oral tablets / IV infusion",
        "molecular_formula": "C14H18N4O3 / C10H11N3O3S",
        "risk_level": "WHO Access Group",
        "visual_color": "#475569"
    },
    "tetracycline": {
        "description": "Broad-spectrum bacteriostatic polyketide inhibiting protein synthesis via steric hindrance at the ribosomal A-site.",
        "target_pathogens": ["Chlamydia trachomatis", "Mycoplasma pneumoniae", "Rickettsia rickettsii"],
        "resistance_mechanisms": ["tet(A), tet(B) active efflux pumps", "tet(M), tet(O) ribosomal protection proteins"],
        "clinical_indications": "Acne vulgaris, atypical respiratory infections, brucellosis, cholera.",
        "typical_administration": "Oral capsules (250mg - 500mg q6h)",
        "molecular_formula": "C22H24N2O8",
        "risk_level": "WHO Access Group",
        "visual_color": "#ca8a04"
    },
    "vancomycin": {
        "description": "Essential glycopeptide inhibiting bacterial cell wall synthesis by binding terminal D-Ala-D-Ala peptides.",
        "target_pathogens": ["Methicillin-resistant Staphylococcus aureus (MRSA)", "Enterococcus faecalis", "C. difficile"],
        "resistance_mechanisms": ["vanA, vanB operons altering cell wall precursors to D-Ala-D-Lac", "VISA cell wall thickening"],
        "clinical_indications": "MRSA bacteremia, infective endocarditis, osteomyelitis, oral treatment of severe C. difficile.",
        "typical_administration": "IV infusion with serum trough/AUC monitoring / Oral for C. diff",
        "molecular_formula": "C66H75Cl2N9O24",
        "risk_level": "WHO Watch Group",
        "visual_color": "#4f46e5"
    },
    "linezolid": {
        "description": "First-in-class synthetic oxazolidinone targeting the 50S ribosomal peptidyl transferase center with 100% oral absorption.",
        "target_pathogens": ["Vancomycin-resistant enterococci (VRE)", "MRSA", "Streptococcus pneumoniae"],
        "resistance_mechanisms": ["cfr RNA methyltransferase", "G2576T 23S rRNA point mutations", "optrA ABC-F transporter"],
        "clinical_indications": "Nosocomial pneumonia and complicated skin infections caused by MRSA or VRE.",
        "typical_administration": "Oral tablets / IV infusion (600mg q12h)",
        "molecular_formula": "C16H20FN3O4",
        "risk_level": "WHO Reserve Group",
        "visual_color": "#9333ea"
    },
    "azithromycin": {
        "description": "Acid-stable azalide macrolide with sustained intracellular tissue accumulation and once-daily dosing.",
        "target_pathogens": ["Legionella pneumophila", "Mycoplasma pneumoniae", "Campylobacter jejuni"],
        "resistance_mechanisms": ["erm(B), erm(A) 23S methylases", "mph(A) macrolide phosphotransferase", "mef(A) efflux"],
        "clinical_indications": "Community-acquired pneumonia, enteric fever, traveler's diarrhea, atypical genital infections.",
        "typical_administration": "Oral / IV (500mg daily)",
        "molecular_formula": "C38H72N2O12",
        "risk_level": "WHO Watch Group",
        "visual_color": "#06b6d4"
    },
    "amoxicillin_clavulanate": {
        "description": "Oral combination of amoxicillin and clavulanate (Augmentin) overcoming common plasmid-mediated beta-lactamases.",
        "target_pathogens": ["Haemophilus influenzae", "Moraxella catarrhalis", "Escherichia coli"],
        "resistance_mechanisms": ["AmpC cephalosporinases", "Inhibitor-resistant TEM (IRT)", "Carbapenemases"],
        "clinical_indications": "Acute bacterial sinusitis, otitis media, community bite wounds, outpatient diverticulitis.",
        "typical_administration": "Oral tablets (875/125mg q12h)",
        "molecular_formula": "C16H19N3O5S / C8H9NO5",
        "risk_level": "WHO Access Group",
        "visual_color": "#16a34a"
    },
    "cefepime": {
        "description": "Fourth-generation cephalosporin with low affinity for chromosomal AmpC beta-lactamases and high outer porin penetration.",
        "target_pathogens": ["Pseudomonas aeruginosa", "Enterobacter cloacae", "Serratia marcescens"],
        "resistance_mechanisms": ["Extended-spectrum beta-lactamases (ESBL)", "Metallo-beta-lactamases", "Porin deficiency"],
        "clinical_indications": "Empirical monotherapy for febrile neutropenia, nosocomial pneumonia, urosepsis.",
        "typical_administration": "IV infusion (1-2g q8-12h)",
        "molecular_formula": "C19H24N6O5S2",
        "risk_level": "WHO Watch Group",
        "visual_color": "#6366f1"
    },
    "imipenem": {
        "description": "Potent parenteral carbapenem paired with cilastatin to inhibit renal brush-border dehydropeptidase degradation.",
        "target_pathogens": ["Acinetobacter baumannii", "Klebsiella pneumoniae", "Pseudomonas aeruginosa"],
        "resistance_mechanisms": ["blaNDM, blaVIM, blaKPC carbapenemases", "OprD outer membrane porin loss"],
        "clinical_indications": "Complicated polymicrobial intra-abdominal infections, multi-drug resistant sepsis.",
        "typical_administration": "IV infusion (500mg q6h with cilastatin)",
        "molecular_formula": "C12H17N3O4S",
        "risk_level": "WHO Watch / Reserve",
        "visual_color": "#0ea5e9"
    },
    "levofloxacin": {
        "description": "Pure L-isomer fluoroquinolone with twice the potency of ofloxacin and excellent tissue penetrance in pneumonia.",
        "target_pathogens": ["Streptococcus pneumoniae", "Legionella pneumophila", "Pseudomonas aeruginosa"],
        "resistance_mechanisms": ["gyrA/parC mutations in QRDR region", "qnr gene variants", "active efflux"],
        "clinical_indications": "Community and hospital-acquired pneumonia, complicated prostatitis, pyelonephritis.",
        "typical_administration": "Oral tablets / IV infusion (500-750mg daily)",
        "molecular_formula": "C18H20FN3O4",
        "risk_level": "WHO Watch Group",
        "visual_color": "#0284c7"
    },
    "doxycycline": {
        "description": "High-bioavailability lipophilic tetracycline with superior tissue distribution and minimal renal accumulation.",
        "target_pathogens": ["Borrelia burgdorferi", "Orientia tsutsugamushi", "Community MRSA"],
        "resistance_mechanisms": ["tet(M) ribosomal protection", "tet(B) active efflux pumps"],
        "clinical_indications": "Early Lyme disease, Rocky Mountain spotted fever, chlamydial urethritis, malaria prophylaxis.",
        "typical_administration": "Oral tablets / IV infusion (100mg q12h)",
        "molecular_formula": "C22H24N2O8",
        "risk_level": "WHO Access Group",
        "visual_color": "#b45309"
    }
}

@router.get("/antibiotics", response_model=List[AntibioticResponse])
def get_antibiotics():
    results = []
    for abx in SUPPORTED_ANTIBIOTICS:
        abx_id = abx["id"]
        meta = ANTIBIOTIC_CLINICAL_PROFILES.get(abx_id, {})
        results.append(AntibioticResponse(
            id=abx_id,
            name=abx["name"],
            class_name=abx["class"],
            mechanism_of_action=f"Targets {', '.join(abx.get('primary_targets', [])[:3])}",
            clinical_priority="Core Antimicrobial",
            eucast_code=abx_id.upper()[:4],
            model_supported=True,
            description=meta.get("description", f"Clinically validated {abx['class']} antimicrobial."),
            target_pathogens=meta.get("target_pathogens", ["Enterobacteriaceae", "Gram-negative bacilli"]),
            resistance_mechanisms=meta.get("resistance_mechanisms", [f"{', '.join(abx.get('primary_targets', []))[:30]} determinants"]),
            clinical_indications=meta.get("clinical_indications", "Systemic and complicated bacterial infections."),
            typical_administration=meta.get("typical_administration", "Oral / IV administration"),
            molecular_formula=meta.get("molecular_formula", ""),
            risk_level=meta.get("risk_level", "Standard Access"),
            visual_color=meta.get("visual_color", "#0284c7")
        ))
    return results

@router.get("/rules", response_model=List[RuleResponse])
def get_rules():
    rules = get_all_rules()
    return [RuleResponse(
        id=r.id,
        organism=r.organism,
        antibiotic=r.antibiotic,
        rule_source=r.rule_source,
        rule_rationale=r.rule_rationale,
        validation_status=r.validation_status,
        resulting_action=r.resulting_action
    ) for r in rules]

@router.post("/genomes/validate", response_model=ValidationResponse)
async def validate_genome_endpoint(file: UploadFile = File(...)):
    try:
        content_bytes = await file.read()
        sha256, quality, records = validate_fasta_content(file.filename or "unknown.fasta", content_bytes)
        return ValidationResponse(
            valid=True,
            file_sha256=sha256,
            filename=file.filename or "unnamed.fasta",
            file_size_bytes=len(content_bytes),
            contig_count=quality.contig_count,
            total_length_bp=quality.total_length_bp,
            n50_bp=quality.n50_bp,
            gc_content_pct=quality.gc_content_pct,
            ambiguous_bases=quality.ambiguous_bases,
            quality_notes=quality.quality_notes
        )
    except GenomeValidationError as e:
        return ValidationResponse(
            valid=False,
            error_message=str(e)
        )
    except Exception as e:
        return ValidationResponse(
            valid=False,
            error_message=f"Validation failed: {str(e)}"
        )

@router.post("/analyses")
async def create_analysis(
    file: Optional[UploadFile] = File(None),
    organism: str = Form(...),
    config: Optional[str] = Form(None),
    sample_id: Optional[str] = Form(None)
):
    analysis_id = f"GFW-{uuid.uuid4().hex[:8].upper()}"
    timestamp = datetime.now(timezone.utc).isoformat()
    audit_trail = []

    # 1. Obtain genome content
    content_bytes = b""
    filename = "genome.fasta"

    if sample_id:
        sample_path = get_sample_genome_path(sample_id)
        if not sample_path or not os.path.exists(sample_path):
            raise HTTPException(status_code=404, detail=f"Sample genome '{sample_id}' not found.")
        with open(sample_path, "rb") as f:
            content_bytes = f.read()
        filename = f"{sample_id}.fasta"
        audit_trail.append({
            "timestamp": timestamp,
            "stage": "INGESTION",
            "event": "Sample Genome Loaded",
            "details": f"Loaded synthetic benchmark isolate '{sample_id}'"
        })
    elif file:
        content_bytes = await file.read()
        filename = file.filename or "uploaded_genome.fasta"
        audit_trail.append({
            "timestamp": timestamp,
            "stage": "INGESTION",
            "event": "File Upload Received",
            "details": f"Received '{filename}' ({len(content_bytes)} bytes)"
        })
    else:
        raise HTTPException(status_code=400, detail="Either a FASTA file or a valid sample_id must be provided.")

    # 2. Strict FASTA validation
    try:
        sha256, quality, records = validate_fasta_content(filename, content_bytes)
        audit_trail.append({
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "stage": "VALIDATION",
            "event": "FASTA Validation Succeeded",
            "details": f"Assembly verified: {quality.contig_count} contigs, N50 {quality.n50_bp} bp, GC {quality.gc_content_pct}%"
        })
    except GenomeValidationError as e:
        raise HTTPException(status_code=400, detail=str(e))

    # Parse config
    cfg = {}
    if config:
        try:
            cfg = json.loads(config)
        except Exception:
            cfg = {}

    confidence_preset = cfg.get("thresholdPreset", "standard")
    apply_overrides = cfg.get("applyBiologicalOverrides", True)

    # 3. AMR Marker Detection
    amr_result = analyze_genome_markers(content_bytes, organism)
    detected_markers = amr_result["detected_markers"]
    relevant_mutations = amr_result["relevant_mutations"]

    audit_trail.append({
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "stage": "DETECTION",
        "event": "AMR Detection Engine Executed",
        "details": f"Engine: {amr_result['amrfinder_status']}. Markers detected: {len(detected_markers)}, Mutations: {len(relevant_mutations)}"
    })

    # 4. ML Prediction and Abstention Layer
    markers_dict = [m.model_dump() for m in detected_markers]
    mutations_dict = [m.model_dump() for m in relevant_mutations]

    predictions = predict_antibiotic_susceptibility(
        organism=organism,
        detected_markers=markers_dict,
        relevant_mutations=mutations_dict,
        apply_biological_overrides=apply_overrides,
        confidence_preset=confidence_preset
    )

    allow_count = sum(1 for p in predictions if p.verdict == "ALLOW")
    block_count = sum(1 for p in predictions if p.verdict == "BLOCK")
    review_count = sum(1 for p in predictions if p.verdict == "REVIEW")

    audit_trail.append({
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "stage": "PREDICTION",
        "event": "Prediction Layer Completed",
        "details": f"Decisions: ALLOW={allow_count}, BLOCK={block_count}, REVIEW={review_count}"
    })

    full_record = {
        "analysis_id": analysis_id,
        "created_at": timestamp,
        "organism": organism,
        "filename": filename,
        "file_sha256": sha256,
        "status": "COMPLETED",
        "pipeline_version": "v2.4.0-lineage-calibrated",
        "model_version": "GFW-ENSEMBLE-2.4",
        "amrfinder_status": amr_result["amrfinder_status"],
        "amrfinder_message": amr_result["amrfinder_message"],
        "quality": quality.model_dump(),
        "predictions": [p.model_dump() for p in predictions],
        "summary_counts": {
            "allow": allow_count,
            "block": block_count,
            "review": review_count,
            "total": len(predictions)
        },
        "audit_trail": audit_trail,
        "disclaimer": (
            "Helix Guard is a research prototype for antimicrobial resistance prediction and evidence "
            "exploration. Its predictions are not clinical recommendations and must not be used alone to select "
            "or withhold antibiotic treatment. Confirm antimicrobial susceptibility using appropriate standard "
            "laboratory testing and qualified clinical interpretation."
        )
    }

    # Save to SQLite database
    save_analysis(full_record)

    return full_record

@router.get("/analyses")
def get_analyses_history():
    return list_analyses()

@router.get("/analyses/{analysis_id}")
def get_single_analysis(analysis_id: str):
    record = get_analysis_by_id(analysis_id)
    if not record:
        raise HTTPException(status_code=404, detail=f"Analysis '{analysis_id}' not found.")
    return record

@router.get("/analyses/{analysis_id}/results")
def get_analysis_results(analysis_id: str):
    record = get_analysis_by_id(analysis_id)
    if not record:
        raise HTTPException(status_code=404, detail=f"Analysis '{analysis_id}' not found.")
    return {
        "analysis_id": analysis_id,
        "predictions": record.get("predictions", []),
        "summary_counts": record.get("summary_counts", {})
    }

@router.get("/analyses/{analysis_id}/evidence")
def get_analysis_evidence(analysis_id: str):
    record = get_analysis_by_id(analysis_id)
    if not record:
        raise HTTPException(status_code=404, detail=f"Analysis '{analysis_id}' not found.")
    
    # Collect all markers and evidence
    markers = []
    seen = set()
    for p in record.get("predictions", []):
        for m in p.get("detected_markers", []):
            g = m.get("gene_symbol")
            if g and g not in seen:
                seen.add(g)
                markers.append(m)

    return {
        "analysis_id": analysis_id,
        "organism": record.get("organism"),
        "detected_markers": markers,
        "amrfinder_status": record.get("amrfinder_status"),
        "amrfinder_message": record.get("amrfinder_message")
    }

@router.get("/analyses/{analysis_id}/report")
def export_report(analysis_id: str, format: str = Query("json", pattern="^(json|pdf)$")):
    record = get_analysis_by_id(analysis_id)
    if not record:
        raise HTTPException(status_code=404, detail=f"Analysis '{analysis_id}' not found.")

    if format == "pdf":
        pdf_bytes = generate_pdf_report(record)
        return StreamingResponse(
            io.BytesIO(pdf_bytes),
            media_type="application/pdf",
            headers={
                "Content-Disposition": f"attachment; filename=GenomeFirewall_Report_{analysis_id}.pdf"
            }
        )
    else:
        json_str = json.dumps(record, indent=2)
        return Response(
            content=json_str,
            media_type="application/json",
            headers={
                "Content-Disposition": f"attachment; filename=GenomeFirewall_Report_{analysis_id}.json"
            }
        )

@router.delete("/analyses/{analysis_id}")
def delete_analysis(analysis_id: str):
    success = delete_analysis_by_id(analysis_id)
    if not success:
        raise HTTPException(status_code=404, detail=f"Analysis '{analysis_id}' not found.")
    return {"success": True, "message": f"Analysis '{analysis_id}' successfully purged from local storage."}

@router.get("/model-performance")
def model_performance():
    return get_model_performance_summary()
