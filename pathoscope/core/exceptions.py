"""
exceptions.py

Consolidated pipeline custom exceptions system to isolate module failures 
and prevent silent crashes.
"""

class PathoScopeError(Exception):
    """Base exception class for all PathoScope AI errors."""
    pass

class PreprocessingError(PathoScopeError):
    """Raised when sequence validation or QC preprocessing fails."""
    pass

class ORFError(PathoScopeError):
    """Raised when open reading frame coordinate prediction fails."""
    pass

class AnnotationError(PathoScopeError):
    """Raised when sequence similarity annotation fails."""
    pass

class PathwayError(PathoScopeError):
    """Raised when pathway/conserved domain mapping fails."""
    pass

class StatisticsError(PathoScopeError):
    """Raised when statistical overrepresentation or GSEA fails."""
    pass

class BiologicalInconsistencyError(StatisticsError):
    """Raised when statistical parameters violate biological/mathematical boundaries."""
    pass

class NormalizationError(PathoScopeError):
    """Raised when Gene ID or library count normalization fails."""
    pass

class DifferentialExpressionError(PathoScopeError):
    """Raised when differential expression analysis fails."""
    pass

class ToolExecutionError(PathoScopeError):
    """Raised when external bioinformatics tools (e.g. fastp, HMMER, DIAMOND) fail."""
    pass
