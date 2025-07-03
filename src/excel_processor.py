"""
Excel marker data processing and validation module.

This module handles reading Excel files containing genetic marker information,
validates the data format, and provides utilities for marker data manipulation.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Set, Tuple, Any
from dataclasses import dataclass
from pathlib import Path
import re

from .config import get_config, get_logger, Constants, validate_file_path, check_file_size


@dataclass
class MarkerInfo:
    """Information about a genetic marker."""
    name: str
    chromosome: str
    position: int
    ref_allele: str
    alt_allele: str
    marker_type: Optional[str] = None
    gene: Optional[str] = None
    annotation: Optional[str] = None


@dataclass
class ExcelValidationResult:
    """Result of Excel file validation."""
    is_valid: bool
    errors: List[str]
    warnings: List[str]
    detected_columns: Dict[str, str]
    row_count: int
    marker_count: int


class ExcelProcessor:
    """Processor for Excel files containing genetic marker data."""
    
    def __init__(self):
        """Initialize Excel processor."""
        self.config = get_config()
        self.logger = get_logger(__name__)
        self._column_mapping_cache = {}
    
    def _detect_column_mapping(self, df: pd.DataFrame) -> Dict[str, str]:
        """
        Automatically detect column mapping for marker data.
        
        Args:
            df: DataFrame with Excel data
            
        Returns:
            Dictionary mapping standard names to actual column names
        """
        columns = df.columns.tolist()
        mapping = {}
        
        # Iterate through each expected column type
        for standard_name, possible_names in Constants.EXCEL_MARKER_COLUMNS.items():
            found_column = None
            
            # Check for exact matches first (case insensitive)
            for col in columns:
                if col.lower() in [name.lower() for name in possible_names]:
                    found_column = col
                    break
            
            # If no exact match, check for partial matches
            if not found_column:
                for col in columns:
                    for possible_name in possible_names:
                        if possible_name.lower() in col.lower() or col.lower() in possible_name.lower():
                            found_column = col
                            break
                    if found_column:
                        break
            
            if found_column:
                mapping[standard_name] = found_column
                self.logger.debug(f"Mapped {standard_name} -> {found_column}")
        
        return mapping
    
    def _validate_chromosome(self, chrom_value: Any) -> Tuple[bool, str]:
        """
        Validate and standardize chromosome value.
        
        Args:
            chrom_value: Raw chromosome value
            
        Returns:
            Tuple of (is_valid, standardized_value)
        """
        if pd.isna(chrom_value):
            return False, "Missing chromosome"
        
        chrom_str = str(chrom_value).strip().upper()
        
        # Remove 'CHR' prefix if present
        if chrom_str.startswith('CHR'):
            chrom_str = chrom_str[3:]
        
        # Check for valid chromosome names
        valid_autosomes = set(map(str, range(1, 100)))  # Support up to 99 autosomes
        valid_sex_chroms = {'X', 'Y', 'Z', 'W', 'MT', 'M'}
        valid_unplaced = {'UN', 'UNPLACED', 'UNKNOWN'}
        
        if chrom_str in valid_autosomes or chrom_str in valid_sex_chroms or chrom_str in valid_unplaced:
            return True, chrom_str
        
        # Try to extract number from mixed strings
        import re
        match = re.search(r'\d+', chrom_str)
        if match:
            num = match.group()
            if num in valid_autosomes:
                return True, num
        
        return False, f"Invalid chromosome: {chrom_value}"
    
    def _validate_position(self, pos_value: Any) -> Tuple[bool, int]:
        """
        Validate genomic position.
        
        Args:
            pos_value: Raw position value
            
        Returns:
            Tuple of (is_valid, position)
        """
        if pd.isna(pos_value):
            return False, 0
        
        try:
            # Handle string values that might contain commas
            if isinstance(pos_value, str):
                pos_value = pos_value.replace(',', '').strip()
            
            position = int(float(pos_value))
            
            if position <= 0:
                return False, position
            
            # Check for unreasonably large positions (> 1 billion)
            if position > 1_000_000_000:
                return False, position
            
            return True, position
            
        except (ValueError, TypeError):
            return False, 0
    
    def _validate_allele(self, allele_value: Any) -> Tuple[bool, str]:
        """
        Validate allele sequence.
        
        Args:
            allele_value: Raw allele value
            
        Returns:
            Tuple of (is_valid, standardized_allele)
        """
        if pd.isna(allele_value):
            return False, ""
        
        allele_str = str(allele_value).strip().upper()
        
        # Empty allele
        if not allele_str:
            return False, ""
        
        # Check for valid DNA bases
        valid_bases = set('ACGTN-')
        if not all(base in valid_bases for base in allele_str):
            return False, allele_str
        
        # Check for reasonable length (< 1000 bp for SNPs/indels)
        if len(allele_str) > 1000:
            return False, allele_str
        
        return True, allele_str
    
    def validate_excel_file(self, filepath: str) -> ExcelValidationResult:
        """
        Validate Excel file format and content.
        
        Args:
            filepath: Path to Excel file
            
        Returns:
            ExcelValidationResult object
        """
        errors = []
        warnings = []
        detected_columns = {}
        row_count = 0
        marker_count = 0
        
        try:
            # Validate file existence and size
            validate_file_path(filepath)
            check_file_size(filepath)
            
            # Check file extension
            if not any(filepath.lower().endswith(ext) for ext in Constants.SUPPORTED_EXCEL_EXTENSIONS):
                errors.append(f"Unsupported file format. Expected: {Constants.SUPPORTED_EXCEL_EXTENSIONS}")
                return ExcelValidationResult(False, errors, warnings, detected_columns, 0, 0)
            
            # Read Excel file
            try:
                df = pd.read_excel(filepath)
            except Exception as e:
                errors.append(f"Failed to read Excel file: {e}")
                return ExcelValidationResult(False, errors, warnings, detected_columns, 0, 0)
            
            row_count = len(df)
            
            if row_count == 0:
                errors.append("Excel file is empty")
                return ExcelValidationResult(False, errors, warnings, detected_columns, 0, 0)
            
            # Detect column mapping
            detected_columns = self._detect_column_mapping(df)
            
            # Check for required columns
            required_fields = ['marker_name', 'chromosome', 'position']
            missing_required = []
            
            for field in required_fields:
                if field not in detected_columns:
                    missing_required.append(field)
            
            if missing_required:
                errors.append(f"Missing required columns: {missing_required}")
                return ExcelValidationResult(False, errors, warnings, detected_columns, row_count, 0)
            
            # Validate data content
            valid_markers = 0
            
            for idx, row in df.iterrows():
                row_errors = []
                
                # Validate marker name
                marker_name = row.get(detected_columns['marker_name'])
                if pd.isna(marker_name) or str(marker_name).strip() == '':
                    row_errors.append(f"Row {idx + 2}: Missing marker name")
                
                # Validate chromosome
                chrom_value = row.get(detected_columns['chromosome'])
                chrom_valid, chrom_msg = self._validate_chromosome(chrom_value)
                if not chrom_valid:
                    row_errors.append(f"Row {idx + 2}: {chrom_msg}")
                
                # Validate position
                pos_value = row.get(detected_columns['position'])
                pos_valid, pos_result = self._validate_position(pos_value)
                if not pos_valid:
                    row_errors.append(f"Row {idx + 2}: Invalid position: {pos_value}")
                
                # Validate alleles if present
                if 'ref_allele' in detected_columns:
                    ref_value = row.get(detected_columns['ref_allele'])
                    ref_valid, ref_result = self._validate_allele(ref_value)
                    if not ref_valid and not pd.isna(ref_value):
                        row_errors.append(f"Row {idx + 2}: Invalid reference allele: {ref_value}")
                
                if 'alt_allele' in detected_columns:
                    alt_value = row.get(detected_columns['alt_allele'])
                    alt_valid, alt_result = self._validate_allele(alt_value)
                    if not alt_valid and not pd.isna(alt_value):
                        row_errors.append(f"Row {idx + 2}: Invalid alternative allele: {alt_value}")
                
                # If no row errors, count as valid marker
                if not row_errors:
                    valid_markers += 1
                else:
                    # Add first few errors as warnings (don't fail entire file)
                    if len(warnings) < 50:  # Limit warnings
                        warnings.extend(row_errors[:3])
            
            marker_count = valid_markers
            
            # Check for reasonable number of valid markers
            if valid_markers == 0:
                errors.append("No valid markers found in file")
            elif valid_markers < row_count * 0.5:
                warnings.append(f"Only {valid_markers}/{row_count} markers are valid (< 50%)")
            
            # Check for duplicates
            if 'marker_name' in detected_columns:
                marker_names = df[detected_columns['marker_name']].dropna()
                duplicates = marker_names[marker_names.duplicated()].unique()
                if len(duplicates) > 0:
                    warnings.append(f"Found {len(duplicates)} duplicate marker names")
            
            self.logger.info(f"Excel validation complete: {valid_markers}/{row_count} valid markers")
            
            is_valid = len(errors) == 0
            return ExcelValidationResult(is_valid, errors, warnings, detected_columns, row_count, marker_count)
            
        except Exception as e:
            error_msg = f"Excel validation failed: {e}"
            self.logger.error(error_msg)
            errors.append(error_msg)
            return ExcelValidationResult(False, errors, warnings, detected_columns, row_count, marker_count)
    
    def read_markers(self, filepath: str, validate: bool = True) -> List[MarkerInfo]:
        """
        Read genetic markers from Excel file.
        
        Args:
            filepath: Path to Excel file
            validate: Whether to validate file first
            
        Returns:
            List of MarkerInfo objects
        """
        if validate:
            validation_result = self.validate_excel_file(filepath)
            if not validation_result.is_valid:
                raise ValueError(f"Excel file validation failed: {validation_result.errors}")
            
            if validation_result.warnings:
                for warning in validation_result.warnings[:10]:  # Show first 10 warnings
                    self.logger.warning(warning)
        
        try:
            df = pd.read_excel(filepath)
            column_mapping = self._detect_column_mapping(df)
            
            markers = []
            
            for idx, row in df.iterrows():
                try:
                    # Extract marker name
                    marker_name = str(row[column_mapping['marker_name']]).strip()
                    if not marker_name or marker_name.lower() == 'nan':
                        continue
                    
                    # Extract and validate chromosome
                    chrom_valid, chromosome = self._validate_chromosome(row[column_mapping['chromosome']])
                    if not chrom_valid:
                        continue
                    
                    # Extract and validate position
                    pos_valid, position = self._validate_position(row[column_mapping['position']])
                    if not pos_valid:
                        continue
                    
                    # Extract alleles if available
                    ref_allele = ""
                    alt_allele = ""
                    
                    if 'ref_allele' in column_mapping:
                        ref_valid, ref_allele = self._validate_allele(row[column_mapping['ref_allele']])
                        if not ref_valid:
                            ref_allele = ""
                    
                    if 'alt_allele' in column_mapping:
                        alt_valid, alt_allele = self._validate_allele(row[column_mapping['alt_allele']])
                        if not alt_valid:
                            alt_allele = ""
                    
                    # Extract optional fields
                    marker_type = None
                    gene = None
                    annotation = None
                    
                    # Look for additional columns
                    for col in df.columns:
                        col_lower = col.lower()
                        if 'type' in col_lower and marker_type is None:
                            marker_type = str(row[col]).strip() if not pd.isna(row[col]) else None
                        elif 'gene' in col_lower and gene is None:
                            gene = str(row[col]).strip() if not pd.isna(row[col]) else None
                        elif any(term in col_lower for term in ['annotation', 'note', 'comment']) and annotation is None:
                            annotation = str(row[col]).strip() if not pd.isna(row[col]) else None
                    
                    marker = MarkerInfo(
                        name=marker_name,
                        chromosome=chromosome,
                        position=position,
                        ref_allele=ref_allele,
                        alt_allele=alt_allele,
                        marker_type=marker_type,
                        gene=gene,
                        annotation=annotation
                    )
                    
                    markers.append(marker)
                    
                except Exception as e:
                    self.logger.warning(f"Skipping row {idx + 2}: {e}")
                    continue
            
            self.logger.info(f"Successfully read {len(markers)} markers from {filepath}")
            return markers
            
        except Exception as e:
            error_msg = Constants.ERROR_MESSAGES['processing_error'].format(
                filepath=filepath, error=str(e)
            )
            self.logger.error(error_msg)
            raise RuntimeError(error_msg) from e
    
    def filter_markers_by_chromosome(self, markers: List[MarkerInfo], 
                                   chromosomes: Set[str]) -> List[MarkerInfo]:
        """
        Filter markers by chromosome.
        
        Args:
            markers: List of MarkerInfo objects
            chromosomes: Set of chromosome names to keep
            
        Returns:
            Filtered list of markers
        """
        filtered_markers = [
            marker for marker in markers 
            if marker.chromosome in chromosomes
        ]
        
        self.logger.info(f"Filtered {len(markers)} -> {len(filtered_markers)} markers by chromosome")
        return filtered_markers
    
    def filter_markers_by_position_range(self, markers: List[MarkerInfo], 
                                       chrom: str, start: int, end: int) -> List[MarkerInfo]:
        """
        Filter markers by genomic position range.
        
        Args:
            markers: List of MarkerInfo objects
            chrom: Chromosome name
            start: Start position
            end: End position
            
        Returns:
            Filtered list of markers
        """
        filtered_markers = [
            marker for marker in markers 
            if marker.chromosome == chrom and start <= marker.position <= end
        ]
        
        self.logger.info(f"Filtered markers in {chrom}:{start}-{end}: {len(filtered_markers)} found")
        return filtered_markers
    
    def get_marker_chromosomes(self, markers: List[MarkerInfo]) -> Set[str]:
        """
        Get set of unique chromosomes from markers.
        
        Args:
            markers: List of MarkerInfo objects
            
        Returns:
            Set of chromosome names
        """
        chromosomes = {marker.chromosome for marker in markers}
        self.logger.info(f"Found markers on {len(chromosomes)} chromosomes: {sorted(chromosomes)}")
        return chromosomes
    
    def get_markers_by_name(self, markers: List[MarkerInfo], names: Set[str]) -> Dict[str, MarkerInfo]:
        """
        Create mapping of marker names to MarkerInfo objects.
        
        Args:
            markers: List of MarkerInfo objects
            names: Set of marker names to include
            
        Returns:
            Dictionary mapping names to MarkerInfo objects
        """
        marker_dict = {}
        
        for marker in markers:
            if marker.name in names:
                marker_dict[marker.name] = marker
        
        self.logger.info(f"Created mapping for {len(marker_dict)}/{len(names)} requested markers")
        return marker_dict
    
    def export_markers_to_excel(self, markers: List[MarkerInfo], output_path: str, 
                               include_stats: bool = False) -> None:
        """
        Export markers to Excel file.
        
        Args:
            markers: List of MarkerInfo objects
            output_path: Output file path
            include_stats: Whether to include additional statistics
        """
        try:
            # Prepare data for export
            data = []
            for marker in markers:
                row = {
                    'Marker_Name': marker.name,
                    'Chromosome': marker.chromosome,
                    'Position': marker.position,
                    'Reference_Allele': marker.ref_allele,
                    'Alternative_Allele': marker.alt_allele,
                    'Marker_Type': marker.marker_type,
                    'Gene': marker.gene,
                    'Annotation': marker.annotation
                }
                data.append(row)
            
            df = pd.DataFrame(data)
            
            # Add summary statistics if requested
            if include_stats:
                summary_data = []
                chromosomes = self.get_marker_chromosomes(markers)
                
                for chrom in sorted(chromosomes):
                    chrom_markers = [m for m in markers if m.chromosome == chrom]
                    summary_data.append({
                        'Chromosome': chrom,
                        'Marker_Count': len(chrom_markers),
                        'Min_Position': min(m.position for m in chrom_markers),
                        'Max_Position': max(m.position for m in chrom_markers)
                    })
                
                summary_df = pd.DataFrame(summary_data)
                
                # Write to Excel with multiple sheets
                with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
                    df.to_excel(writer, sheet_name='Markers', index=False)
                    summary_df.to_excel(writer, sheet_name='Summary', index=False)
            else:
                df.to_excel(output_path, index=False)
            
            self.logger.info(f"Exported {len(markers)} markers to {output_path}")
            
        except Exception as e:
            error_msg = f"Failed to export markers to {output_path}: {e}"
            self.logger.error(error_msg)
            raise RuntimeError(error_msg) from e
    
    def validate_marker_names(self, markers: List[MarkerInfo]) -> Tuple[Set[str], Set[str]]:
        """
        Validate marker names for duplicates and format issues.
        
        Args:
            markers: List of MarkerInfo objects
            
        Returns:
            Tuple of (valid_names, problematic_names)
        """
        all_names = [marker.name for marker in markers]
        name_counts = {}
        
        for name in all_names:
            name_counts[name] = name_counts.get(name, 0) + 1
        
        valid_names = set()
        problematic_names = set()
        
        for name, count in name_counts.items():
            if count > 1:
                problematic_names.add(name)
                self.logger.warning(f"Duplicate marker name: {name} (appears {count} times)")
            elif not name or name.isspace():
                problematic_names.add(name)
                self.logger.warning(f"Invalid marker name: '{name}'")
            else:
                valid_names.add(name)
        
        self.logger.info(f"Marker name validation: {len(valid_names)} valid, {len(problematic_names)} problematic")
        return valid_names, problematic_names