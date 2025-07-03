"""
Configuration management and shared constants for the genetic marker analysis system.

This module provides centralized configuration management, logging setup,
and shared constants used across all components of the analysis system.
"""

import logging
import os
from typing import Dict, Any, Optional
from dataclasses import dataclass
from pathlib import Path


@dataclass
class AnalysisConfig:
    """Configuration class for genetic marker analysis parameters."""
    
    # File processing settings
    max_file_size_mb: int = 500
    chunk_size: int = 10000
    proximity_bp: int = 100
    
    # Performance settings
    enable_caching: bool = True
    cache_size_mb: int = 100
    max_workers: int = 4
    
    # Analysis settings
    min_call_rate: float = 0.8
    min_maf: float = 0.01
    max_missing_rate: float = 0.2
    
    # Output settings
    output_format: str = "xlsx"
    include_detailed_stats: bool = True
    generate_plots: bool = False
    
    # Logging settings
    log_level: str = "INFO"
    log_file: Optional[str] = None
    
    def validate(self) -> None:
        """Validate configuration parameters."""
        if self.max_file_size_mb <= 0:
            raise ValueError("max_file_size_mb must be positive")
        if self.chunk_size <= 0:
            raise ValueError("chunk_size must be positive")
        if self.proximity_bp < 0:
            raise ValueError("proximity_bp must be non-negative")
        if not 0 <= self.min_call_rate <= 1:
            raise ValueError("min_call_rate must be between 0 and 1")
        if not 0 <= self.min_maf <= 0.5:
            raise ValueError("min_maf must be between 0 and 0.5")
        if not 0 <= self.max_missing_rate <= 1:
            raise ValueError("max_missing_rate must be between 0 and 1")


class ConfigManager:
    """Manages configuration settings for the genetic marker analysis system."""
    
    def __init__(self, config_file: Optional[str] = None):
        """
        Initialize configuration manager.
        
        Args:
            config_file: Optional path to configuration file
        """
        self.config = AnalysisConfig()
        self._config_file = config_file
        self._logger = None
        
        if config_file and os.path.exists(config_file):
            self.load_config(config_file)
        
        self.setup_logging()
    
    def load_config(self, config_file: str) -> None:
        """
        Load configuration from file.
        
        Args:
            config_file: Path to configuration file
        """
        try:
            import json
            with open(config_file, 'r') as f:
                config_data = json.load(f)
            
            # Update config with loaded values
            for key, value in config_data.items():
                if hasattr(self.config, key):
                    setattr(self.config, key, value)
            
            self.config.validate()
            
        except Exception as e:
            raise RuntimeError(f"Failed to load configuration from {config_file}: {e}")
    
    def save_config(self, config_file: str) -> None:
        """
        Save current configuration to file.
        
        Args:
            config_file: Path to save configuration file
        """
        try:
            import json
            config_dict = {
                field.name: getattr(self.config, field.name)
                for field in self.config.__dataclass_fields__.values()
            }
            
            with open(config_file, 'w') as f:
                json.dump(config_dict, f, indent=2)
                
        except Exception as e:
            raise RuntimeError(f"Failed to save configuration to {config_file}: {e}")
    
    def setup_logging(self) -> None:
        """Setup logging configuration."""
        log_level = getattr(logging, self.config.log_level.upper())
        
        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # Setup root logger
        root_logger = logging.getLogger()
        root_logger.setLevel(log_level)
        
        # Clear existing handlers
        root_logger.handlers.clear()
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(log_level)
        console_handler.setFormatter(formatter)
        root_logger.addHandler(console_handler)
        
        # File handler if specified
        if self.config.log_file:
            try:
                file_handler = logging.FileHandler(self.config.log_file)
                file_handler.setLevel(log_level)
                file_handler.setFormatter(formatter)
                root_logger.addHandler(file_handler)
            except Exception as e:
                print(f"Warning: Could not setup file logging: {e}")
        
        self._logger = logging.getLogger(__name__)
        self._logger.info("Logging configured successfully")
    
    def get_logger(self, name: str) -> logging.Logger:
        """
        Get a logger instance.
        
        Args:
            name: Logger name
            
        Returns:
            Logger instance
        """
        return logging.getLogger(name)
    
    def update_config(self, **kwargs) -> None:
        """
        Update configuration parameters.
        
        Args:
            **kwargs: Configuration parameters to update
        """
        for key, value in kwargs.items():
            if hasattr(self.config, key):
                setattr(self.config, key, value)
            else:
                raise ValueError(f"Unknown configuration parameter: {key}")
        
        self.config.validate()
        
        # Reconfigure logging if log settings changed
        if any(key.startswith('log_') for key in kwargs):
            self.setup_logging()


# Shared constants
class Constants:
    """Shared constants for the genetic marker analysis system."""
    
    # VCF format constants
    VCF_HEADER_PREFIX = "##"
    VCF_COLUMN_HEADER = "#CHROM"
    VCF_REQUIRED_COLUMNS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    
    # Genotype encoding
    MISSING_GENOTYPE = "./."
    HOMOZYGOUS_REF = "0/0"
    HETEROZYGOUS = "0/1"
    HOMOZYGOUS_ALT = "1/1"
    
    # Excel column headers
    EXCEL_MARKER_COLUMNS = {
        'marker_name': ['Marker', 'MarkerName', 'ID', 'Name'],
        'chromosome': ['Chr', 'Chromosome', 'CHROM', 'CHR'],
        'position': ['Pos', 'Position', 'POS', 'BP'],
        'ref_allele': ['Ref', 'Reference', 'REF', 'Allele1'],
        'alt_allele': ['Alt', 'Alternative', 'ALT', 'Allele2']
    }
    
    # File extensions
    SUPPORTED_VCF_EXTENSIONS = ['.vcf', '.vcf.gz']
    SUPPORTED_EXCEL_EXTENSIONS = ['.xlsx', '.xls']
    
    # Error messages
    ERROR_MESSAGES = {
        'file_not_found': "File not found: {filepath}",
        'invalid_file_format': "Invalid file format: {filepath}",
        'file_too_large': "File too large (>{max_size}MB): {filepath}",
        'invalid_vcf_format': "Invalid VCF format in file: {filepath}",
        'no_genotype_data': "No genotype data found in VCF file: {filepath}",
        'invalid_excel_format': "Invalid Excel format or missing required columns: {filepath}",
        'memory_error': "Insufficient memory to process file: {filepath}",
        'processing_error': "Error processing file {filepath}: {error}"
    }
    
    # Performance thresholds
    PERFORMANCE_THRESHOLDS = {
        'large_file_mb': 50,
        'very_large_file_mb': 200,
        'max_cache_entries': 1000,
        'batch_size_small': 1000,
        'batch_size_large': 10000
    }


# Global configuration instance
_config_manager = None


def get_config_manager(config_file: Optional[str] = None) -> ConfigManager:
    """
    Get or create the global configuration manager instance.
    
    Args:
        config_file: Optional path to configuration file
        
    Returns:
        ConfigManager instance
    """
    global _config_manager
    
    if _config_manager is None:
        _config_manager = ConfigManager(config_file)
    
    return _config_manager


def get_config() -> AnalysisConfig:
    """
    Get the current analysis configuration.
    
    Returns:
        AnalysisConfig instance
    """
    return get_config_manager().config


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance.
    
    Args:
        name: Logger name
        
    Returns:
        Logger instance
    """
    return get_config_manager().get_logger(name)


# Utility functions
def validate_file_path(filepath: str, check_exists: bool = True) -> Path:
    """
    Validate and convert file path to Path object.
    
    Args:
        filepath: File path to validate
        check_exists: Whether to check if file exists
        
    Returns:
        Path object
        
    Raises:
        FileNotFoundError: If file doesn't exist and check_exists is True
        ValueError: If filepath is invalid
    """
    if not filepath:
        raise ValueError("File path cannot be empty")
    
    path = Path(filepath)
    
    if check_exists and not path.exists():
        raise FileNotFoundError(Constants.ERROR_MESSAGES['file_not_found'].format(filepath=filepath))
    
    return path


def get_file_size_mb(filepath: str) -> float:
    """
    Get file size in megabytes.
    
    Args:
        filepath: Path to file
        
    Returns:
        File size in MB
    """
    path = validate_file_path(filepath)
    return path.stat().st_size / (1024 * 1024)


def check_file_size(filepath: str, max_size_mb: Optional[float] = None) -> None:
    """
    Check if file size is within limits.
    
    Args:
        filepath: Path to file
        max_size_mb: Maximum file size in MB (uses config default if None)
        
    Raises:
        ValueError: If file is too large
    """
    if max_size_mb is None:
        max_size_mb = get_config().max_file_size_mb
    
    file_size = get_file_size_mb(filepath)
    
    if file_size > max_size_mb:
        raise ValueError(Constants.ERROR_MESSAGES['file_too_large'].format(
            filepath=filepath, max_size=max_size_mb
        ))