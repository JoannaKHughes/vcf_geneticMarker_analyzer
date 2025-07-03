"""
Genetic Marker Analysis System

A modular system for analyzing genetic markers against VCF genotype data
with optimizations for large files and comprehensive reporting capabilities.
"""

__version__ = "1.0.0"
__author__ = "Genetic Analysis Team"

# Import main components for easy access
from .config import (
    get_config_manager, get_config, get_logger, 
    AnalysisConfig, ConfigManager, Constants
)

from .vcf_processor import (
    VCFProcessor, VCFHeader, VCFRecord, GenotypeStats, VCFCache
)

from .excel_processor import (
    ExcelProcessor, MarkerInfo, ExcelValidationResult
)

from .marker_analyzer import (
    MarkerAnalyzer, MarkerMatch, AnalysisResult, BatchAnalysisResults,
    ProximityMatcher
)

from .report_generator import ReportGenerator

from .main_controller import (
    GeneticMarkerAnalysisController, ProgressTracker
)

# Define what's available when importing with "from package import *"
__all__ = [
    # Configuration
    'get_config_manager', 'get_config', 'get_logger',
    'AnalysisConfig', 'ConfigManager', 'Constants',
    
    # VCF Processing
    'VCFProcessor', 'VCFHeader', 'VCFRecord', 'GenotypeStats', 'VCFCache',
    
    # Excel Processing
    'ExcelProcessor', 'MarkerInfo', 'ExcelValidationResult',
    
    # Marker Analysis
    'MarkerAnalyzer', 'MarkerMatch', 'AnalysisResult', 'BatchAnalysisResults',
    'ProximityMatcher',
    
    # Report Generation
    'ReportGenerator',
    
    # Main Controller
    'GeneticMarkerAnalysisController', 'ProgressTracker',
    
    # Package metadata
    '__version__', '__author__'
]

# Quick access function for common use cases
def analyze_markers(vcf_file: str, excel_file: str, output_dir: str, 
                   config_file: str = None, **kwargs):
    """
    Convenience function to run complete genetic marker analysis.
    
    Args:
        vcf_file: Path to VCF file
        excel_file: Path to Excel marker file
        output_dir: Output directory for results
        config_file: Optional configuration file path
        **kwargs: Additional parameters for the analysis
        
    Returns:
        Dictionary with analysis results
    """
    controller = GeneticMarkerAnalysisController(config_file=config_file)
    return controller.run_analysis(vcf_file, excel_file, output_dir, **kwargs)

def quick_analysis(vcf_file: str, excel_file: str, output_file: str,
                  config_file: str = None):
    """
    Convenience function to run quick analysis with minimal output.
    
    Args:
        vcf_file: Path to VCF file
        excel_file: Path to Excel marker file
        output_file: Output file path
        config_file: Optional configuration file path
        
    Returns:
        Dictionary with basic analysis results
    """
    controller = GeneticMarkerAnalysisController(config_file=config_file)
    return controller.run_quick_analysis(vcf_file, excel_file, output_file)

def get_file_info(vcf_file: str, excel_file: str, config_file: str = None):
    """
    Convenience function to get information about input files.
    
    Args:
        vcf_file: Path to VCF file
        excel_file: Path to Excel marker file
        config_file: Optional configuration file path
        
    Returns:
        Dictionary with file information
    """
    controller = GeneticMarkerAnalysisController(config_file=config_file)
    return controller.get_file_info(vcf_file, excel_file)