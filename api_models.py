from typing import List, Dict, Any, Optional
from pydantic import BaseModel

class HealthResponse(BaseModel):
    status: str
    service: str
    version: str
    docker_status: bool
    docker_message: str
    amrfinder_engine: str
    database_connected: bool

class OrganismResponse(BaseModel):
    id: str
    scientific_name: str
    common_name: str
    gram_type: str
    estimated_genome_mbp: str
    typical_gc_pct: str
    priority_pathogen: str
    supported_antibiotics_count: int

class AntibioticResponse(BaseModel):
    id: str
    name: str
    class_name: str
    mechanism_of_action: str
    clinical_priority: str
    eucast_code: str
    model_supported: bool
    description: str = ""
    target_pathogens: List[str] = []
    resistance_mechanisms: List[str] = []
    clinical_indications: str = ""
    typical_administration: str = "Oral / IV"
    molecular_formula: str = ""
    risk_level: str = "Standard Access"
    visual_color: str = "#0284c7"

class ValidationResponse(BaseModel):
    valid: bool
    file_sha256: str = ""
    filename: str = ""
    file_size_bytes: int = 0
    contig_count: int = 0
    total_length_bp: int = 0
    n50_bp: int = 0
    gc_content_pct: float = 0.0
    ambiguous_bases: int = 0
    quality_notes: List[str] = []
    error_message: Optional[str] = None

class AnalysisInitiateResponse(BaseModel):
    analysis_id: str
    status: str
    message: str

class RuleResponse(BaseModel):
    id: str
    organism: str
    antibiotic: str
    rule_source: str
    rule_rationale: str
    validation_status: str
    resulting_action: str
