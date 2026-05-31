import yaml
from pathlib import Path
from typing import List, Optional
from pydantic import BaseModel, Field, field_validator

class PipelineConfig(BaseModel):
    name: str = "PathoScope AI Default Run"
    version: str = "1.0.0"
    random_seed: int = 42

class PreprocessingConfig(BaseModel):
    min_length: int = Field(100, ge=1)
    max_length: int = Field(500000, ge=1)
    max_ambiguous_pct: float = Field(5.0, ge=0.0, le=100.0)
    handle_duplicate_headers: str = "reject"
    remove_duplicate_sequences: bool = False
    
    # Advanced FASTQ Preprocessing additions
    adapter_forward: str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    adapter_reverse: str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    paired_end: bool = False
    fastq_r2_path: Optional[str] = None
    min_mean_qscore: int = Field(30, ge=0, le=60)
    low_complexity_filter: bool = True
    max_reads_cap: int = Field(20000, ge=100)


    @field_validator('max_length')
    def validate_lengths(cls, v, info):
        return v

class OrfPredictionConfig(BaseModel):
    translation_table: int = Field(1, ge=1, le=33)
    min_orf_length_aa: int = Field(30, ge=1)
    start_codons: List[str] = ["ATG", "GTG", "TTG"]
    stop_codons: List[str] = ["TAA", "TAG", "TGA"]
    overlap_threshold_bp: int = Field(50, ge=0)
    resolve_nested: bool = True
    overlap_resolution_policy: str = "keep_all_flag"

class AnnotationConfig(BaseModel):
    alignment_engine: str = "diamond"
    local_db_path: str = "data/reference/viral_proteins.dmnd"
    remote_fallback: bool = False
    eval_threshold: float = Field(1e-5, gt=0.0)
    identity_threshold: float = Field(30.0, ge=0.0, le=100.0)
    coverage_threshold: float = Field(50.0, ge=0.0, le=100.0)

class DomainSearchConfig(BaseModel):
    hmmer_db_path: str = "data/reference/Pfam-A.hmm"
    eval_threshold: float = Field(1e-4, gt=0.0)

class PathwayMappingConfig(BaseModel):
    db_cache_path: str = "data/reference/pathways_cache.db"
    use_kegg_api: bool = True
    kegg_organism: str = "ko"

class StatisticsConfig(BaseModel):
    fdr_threshold: float = Field(0.05, gt=0.0, le=1.0)
    bg_universe_size: int = Field(10000, ge=1)
    
    # Advanced Functional Genomics additions
    log2_fc_threshold: float = Field(1.5, ge=0.0)
    control_replicates: Optional[List[str]] = None
    treated_replicates: Optional[List[str]] = None


class AiInterpretationConfig(BaseModel):
    provider: str = "ollama"
    model_name: str = "llama3"
    api_key_env_var: str = "OPENAI_API_KEY"
    temperature: float = Field(0.2, ge=0.0, le=2.0)

class ReportingConfig(BaseModel):
    theme: str = "dark"
    pdf_generation: bool = True

class AppConfig(BaseModel):
    pipeline: PipelineConfig = Field(default_factory=PipelineConfig)
    preprocessing: PreprocessingConfig = Field(default_factory=PreprocessingConfig)
    orf_prediction: OrfPredictionConfig = Field(default_factory=OrfPredictionConfig)
    annotation: AnnotationConfig = Field(default_factory=AnnotationConfig)
    domain_search: DomainSearchConfig = Field(default_factory=DomainSearchConfig)
    pathway_mapping: PathwayMappingConfig = Field(default_factory=PathwayMappingConfig)
    statistics: StatisticsConfig = Field(default_factory=StatisticsConfig)
    ai_interpretation: AiInterpretationConfig = Field(default_factory=AiInterpretationConfig)
    reporting: ReportingConfig = Field(default_factory=ReportingConfig)


def load_config(config_path: Path) -> AppConfig:
    """
    Reads a YAML file and parses it into a fully validated AppConfig Pydantic model.
    """
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
    with open(config_path, "r", encoding="utf-8") as f:
        try:
            raw_data = yaml.safe_load(f) or {}
            # Standardize empty values to empty dicts so default values are populated
            if not isinstance(raw_data, dict):
                raw_data = {}
            return AppConfig(**raw_data)
        except yaml.YAMLError as ye:
            raise ValueError(f"Invalid YAML syntax in configuration file: {ye}")
        except Exception as e:
            raise ValueError(f"Configuration validation failed: {e}")
