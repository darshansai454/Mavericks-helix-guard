from typing import List, Dict, Any, Optional
from pydantic import BaseModel

class BiologicalRule(BaseModel):
    id: str
    organism: str
    antibiotic: str
    rule_source: str
    rule_rationale: str
    validation_status: str
    resulting_action: str  # 'FORCED_BLOCK' | 'AUDIT_FLAG'

# Validated biological intrinsic resistance rules
# Reference: EUCAST Expert Rules on Antimicrobial Susceptibility Testing v3.3
# and CLSI M100-Ed34 Performance Standards.
INTRINSIC_RULES: List[BiologicalRule] = [
    BiologicalRule(
        id="RULE-KP-AMP",
        organism="Klebsiella pneumoniae",
        antibiotic="Ampicillin",
        rule_source="EUCAST Expert Rules v3.3 / CLSI M100 Table 1E",
        rule_rationale="Klebsiella pneumoniae species possess a chromosomally encoded class A beta-lactamase (blaSHV-1) conferring universal intrinsic clinical resistance to ampicillin, amoxicillin, and ticarcillin.",
        validation_status="Clinically Validated & Internationally Certified",
        resulting_action="FORCED_BLOCK"
    ),
    BiologicalRule(
        id="RULE-PA-AMP",
        organism="Pseudomonas aeruginosa",
        antibiotic="Ampicillin",
        rule_source="EUCAST Expert Rules v3.3 Table 4",
        rule_rationale="P. aeruginosa is intrinsically resistant to ampicillin, amoxicillin-clavulanic acid, cefazolin, and cefotaxime due to low outer membrane permeability and an inducible chromosomal AmpC beta-lactamase.",
        validation_status="Clinically Validated & Internationally Certified",
        resulting_action="FORCED_BLOCK"
    ),
    BiologicalRule(
        id="RULE-PA-SXT",
        organism="Pseudomonas aeruginosa",
        antibiotic="Trimethoprim/Sulfamethoxazole",
        rule_source="EUCAST Expert Rules v3.3 Table 4",
        rule_rationale="P. aeruginosa possesses constitutive multidrug efflux systems (MexAB-OprM) conferring clinical intrinsic resistance to trimethoprim and sulfonamides.",
        validation_status="Clinically Validated & Internationally Certified",
        resulting_action="FORCED_BLOCK"
    ),
    BiologicalRule(
        id="RULE-PA-CRO",
        organism="Pseudomonas aeruginosa",
        antibiotic="Ceftriaxone",
        rule_source="EUCAST Expert Rules v3.3 Table 4",
        rule_rationale="P. aeruginosa is intrinsically resistant to narrow and broad-spectrum 3rd-generation cephalosporins (ceftriaxone, cefotaxime) via baseline AmpC cephalosporinase expression.",
        validation_status="Clinically Validated & Internationally Certified",
        resulting_action="FORCED_BLOCK"
    ),
    BiologicalRule(
        id="RULE-AB-AMP",
        organism="Acinetobacter baumannii",
        antibiotic="Ampicillin",
        rule_source="EUCAST Expert Rules v3.3 Table 5",
        rule_rationale="A. baumannii harbors intrinsic chromosomal oxacillinases (blaOXA-51-like) and intrinsic AmpC beta-lactamases conferring resistance to aminopenicillins.",
        validation_status="Clinically Validated & Internationally Certified",
        resulting_action="FORCED_BLOCK"
    ),
    BiologicalRule(
        id="RULE-AB-CRO",
        organism="Acinetobacter baumannii",
        antibiotic="Ceftriaxone",
        rule_source="EUCAST Expert Rules v3.3 Table 5",
        rule_rationale="A. baumannii produces baseline chromosomal AmpC enzymes and has poor outer membrane porin permeability leading to high-level intrinsic resistance to ceftriaxone.",
        validation_status="Clinically Validated & Internationally Certified",
        resulting_action="FORCED_BLOCK"
    ),
    BiologicalRule(
        id="RULE-SA-COL",
        organism="Staphylococcus aureus",
        antibiotic="Colistin",
        rule_source="CLSI M100 Appendix B",
        rule_rationale="Gram-positive bacteria (including S. aureus) lack the outer lipopolysaccharide membrane targeted by polymyxins, rendering colistin biologically inactive.",
        validation_status="Clinically Validated & Internationally Certified",
        resulting_action="FORCED_BLOCK"
    ),
    BiologicalRule(
        id="RULE-SA-ATM",
        organism="Staphylococcus aureus",
        antibiotic="Aztreonam",
        rule_source="EUCAST Expert Rules v3.3 Table 2",
        rule_rationale="Monobactams (aztreonam) bind specifically to penicillin-binding protein 3 (PBP3) of Gram-negative bacilli and have no affinity for Gram-positive PBPs.",
        validation_status="Clinically Validated & Internationally Certified",
        resulting_action="FORCED_BLOCK"
    ),
    BiologicalRule(
        id="RULE-EC-AMP",
        organism="Enterobacter cloacae",
        antibiotic="Ampicillin",
        rule_source="EUCAST Expert Rules v3.3 Table 3",
        rule_rationale="Enterobacter cloacae carries inducible chromosomal AmpC cephalosporinases conferring natural resistance to ampicillin, amoxicillin/clavulanic acid, and 1st/2nd generation cephalosporins.",
        validation_status="Clinically Validated & Internationally Certified",
        resulting_action="FORCED_BLOCK"
    )
]

def get_all_rules() -> List[BiologicalRule]:
    return INTRINSIC_RULES

def check_intrinsic_override(organism: str, antibiotic: str) -> Optional[BiologicalRule]:
    """
    Evaluates whether an intrinsic biological resistance rule applies to this organism and drug.
    Matches are case-insensitive and match on standardized scientific names.
    """
    org_norm = organism.strip().lower()
    abx_norm = antibiotic.strip().lower()
    
    for rule in INTRINSIC_RULES:
        if rule.organism.lower() in org_norm or org_norm in rule.organism.lower():
            if rule.antibiotic.lower() in abx_norm or abx_norm in rule.antibiotic.lower():
                return rule
    return None
