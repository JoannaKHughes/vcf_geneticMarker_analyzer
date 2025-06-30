#!/usr/bin/env python3
"""
Genetic Marker Analysis Script

This script analyzes genetic markers from an Excel workbook against a large collection of VCF files.
It performs exact matching, proximity checking, and frequency calculations for genetic variants.
Supports both uncompressed (.vcf) and gzip-compressed (.vcf.gz) VCF files.

Author: Genetic Analysis Team
Version: 1.0
"""

import argparse
import gzip
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional, Any
from collections import defaultdict, Counter
import time

import pandas as pd
from tqdm import tqdm
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('genetic_marker_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class GeneticMarkerAnalyzer:
    """
    Main class for genetic marker analysis against VCF files.
    """
    
    def __init__(self, excel_path: str, vcf_folder: str, output_path: str):
        """
        Initialize the analyzer with input and output paths.
        
        Args:
            excel_path: Path to Excel workbook with marker data
            vcf_folder: Path to folder containing VCF files
            output_path: Path for output Excel file
        """
        self.excel_path = Path(excel_path)
        self.vcf_folder = Path(vcf_folder)
        self.output_path = Path(output_path)
        
        # Data structures
        self.excel_data = {}  # Store data from each Excel sheet
        self.vcf_files = []   # List of VCF files to process
        self.vcf_identifiers = {}  # Map filename to 10-char identifier
        self.marker_results = {}  # Store analysis results
        self.vcf_data_cache = {}  # Cache VCF data for efficiency
        
        # Statistics
        self.stats = {
            'total_markers': 0,
            'total_vcf_files': 0,
            'exact_matches': 0,
            'proximity_matches': 0,
            'missing_markers': 0
        }
        
        logger.info(f"Initialized analyzer: Excel={excel_path}, VCF_folder={vcf_folder}, Output={output_path}")
    
    def validate_inputs(self) -> bool:
        """
        Validate that input files and folders exist.
        
        Returns:
            bool: True if all inputs are valid
        """
        if not self.excel_path.exists():
            logger.error(f"Excel file not found: {self.excel_path}")
            return False
            
        if not self.vcf_folder.exists():
            logger.error(f"VCF folder not found: {self.vcf_folder}")
            return False
            
        # Check for VCF files
        vcf_files = list(self.vcf_folder.glob("*.vcf")) + list(self.vcf_folder.glob("*.vcf.gz"))
        if not vcf_files:
            logger.error(f"No VCF files found in {self.vcf_folder}")
            return False
            
        logger.info(f"Found {len(vcf_files)} VCF files")
        return True
    
    def read_excel_data(self) -> Dict[str, pd.DataFrame]:
        """
        Read genetic marker data from Excel workbook.
        
        Returns:
            Dict mapping sheet names to DataFrames
        """
        logger.info("Reading Excel workbook...")
        
        try:
            # Read all three required sheets
            required_sheets = ["ISAG Core", "ISAG Back up", "X_Y"]
            excel_data = {}
            
            for sheet_name in required_sheets:
                try:
                    df = pd.read_excel(self.excel_path, sheet_name=sheet_name)
                    
                    # Validate required columns
                    required_columns = ['CHROM', 'POS', 'IDs', 'List', 'Reference Allele', 'Variant Allele', 'Strand', 'Sequence']
                    missing_columns = [col for col in required_columns if col not in df.columns]
                    
                    if missing_columns:
                        logger.warning(f"Missing columns in {sheet_name}: {missing_columns}")
                    
                    # Clean and standardize data
                    df['CHROM'] = df['CHROM'].astype(str)
                    df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
                    
                    # Remove rows with invalid positions
                    df = df.dropna(subset=['POS'])
                    df['POS'] = df['POS'].astype(int)
                    
                    excel_data[sheet_name] = df
                    logger.info(f"Loaded {len(df)} markers from sheet '{sheet_name}'")
                    
                except Exception as e:
                    logger.error(f"Error reading sheet '{sheet_name}': {e}")
                    continue
            
            if not excel_data:
                raise ValueError("No valid sheets found in Excel file")
                
            self.excel_data = excel_data
            self.stats['total_markers'] = sum(len(df) for df in excel_data.values())
            
            return excel_data
            
        except Exception as e:
            logger.error(f"Error reading Excel file: {e}")
            raise
    
    def discover_vcf_files(self) -> List[Path]:
        """
        Discover and catalog VCF files in the input folder.
        
        Returns:
            List of VCF file paths
        """
        logger.info("Discovering VCF files...")
        
        # Find all VCF files (both .vcf and .vcf.gz)
        vcf_files = []
        vcf_files.extend(self.vcf_folder.glob("*.vcf"))
        vcf_files.extend(self.vcf_folder.glob("*.vcf.gz"))
        
        # Sort for consistent processing order
        vcf_files = sorted(vcf_files)
        
        # Generate 10-character identifiers from filenames
        identifiers = {}
        for vcf_file in vcf_files:
            # Handle both .vcf and .vcf.gz files for identifier extraction
            filename = vcf_file.name
            if filename.endswith('.vcf.gz'):
                filename = filename[:-7]  # Remove .vcf.gz extension
            elif filename.endswith('.vcf'):
                filename = filename[:-4]  # Remove .vcf extension
            
            # Take first 10 characters as identifier
            identifier = filename[:10]
            identifiers[vcf_file] = identifier
            
            logger.debug(f"VCF file: {vcf_file.name} -> Identifier: {identifier}")
        
        self.vcf_files = vcf_files
        self.vcf_identifiers = identifiers
        self.stats['total_vcf_files'] = len(vcf_files)
        
        logger.info(f"Found {len(vcf_files)} VCF files")
        return vcf_files
    
    def parse_vcf_file(self, vcf_path: Path) -> Dict[Tuple[str, int], Dict]:
        """
        Parse a single VCF file and extract relevant variant information.
        
        Args:
            vcf_path: Path to VCF file
            
        Returns:
            Dict mapping (CHROM, POS) tuples to variant information
        """
        variants = {}
        
        try:
            # Handle both compressed and uncompressed VCF files
            if vcf_path.name.endswith('.vcf.gz'):
                file_opener = gzip.open
                mode = 'rt'
            else:
                file_opener = open
                mode = 'r'
            
            
            with file_opener(vcf_path, mode) as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    # Skip header lines
                    if line.startswith('#'):
                        continue
                    
                    # Parse variant line
                    fields = line.split('\t')
                    if len(fields) < 9:  # Minimum VCF format requirements
                        logger.warning(f"Invalid VCF line in {vcf_path.name}:{line_num}")
                        continue
                    
                    try:
                        chrom = fields[0]
                        pos = int(fields[1])
                        ref = fields[3]
                        alt = fields[4]
                        
                        # Store variant information
                        variant_key = (chrom, pos)
                        variants[variant_key] = {
                            'CHROM': chrom,
                            'POS': pos,
                            'REF': ref,
                            'ALT': alt,
                            'samples': fields[9:] if len(fields) > 9 else [],
                            'line_fields': fields
                        }
                        
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Error parsing line {line_num} in {vcf_path.name}: {e}")
                        continue
                        
        except Exception as e:
            logger.error(f"Error reading VCF file {vcf_path}: {e}")
            return {}
        
        logger.debug(f"Parsed {len(variants)} variants from {vcf_path.name}")
        return variants
    
    def calculate_allele_frequencies(self, variant_data: Dict, vcf_identifier: str) -> Dict[str, float]:
        """
        Calculate allele frequencies using the exact logic specified.
        
        Args:
            variant_data: Variant information from VCF
            vcf_identifier: 10-character VCF identifier
            
        Returns:
            Dict with frequency calculations
        """
        fields = variant_data['line_fields']
        
        # Initialize counters
        homo_ref = 0
        homo_alt = 0
        het = 0
        missing_genotypes = 0
        total_allele = 0
        
        # Process sample data starting at column 9 (0-indexed)
        for i in range(9, len(fields)):
            total_allele += 1
            genotype_field = fields[i]
            
            # Extract genotype (first field before colon)
            g = genotype_field.split(':')[0]
            
            if g == '0/0':
                homo_ref += 1
                total_allele += 2
            elif g == '0/1' or g == '1/0':
                het += 1
                total_allele += 2
            elif g == '1/1':
                homo_alt += 1
                total_allele += 2
            else:
                missing_genotypes += 1
        
        # Calculate frequencies
        if total_allele > 0:
            afreq = (homo_alt * 2 + het) / total_allele  # Alt allele frequency
            rfreq = (homo_ref * 2 + het) / total_allele  # Ref allele frequency
        else:
            afreq = 0.0
            rfreq = 0.0
        
        return {
            f'freq_ref_of_{vcf_identifier}': rfreq,
            f'freq_alt_of_{vcf_identifier}': afreq,
            'homo_ref': homo_ref,
            'homo_alt': homo_alt,
            'het': het,
            'missing_genotypes': missing_genotypes,
            'total_samples': len(fields) - 9,
            'total_alleles': total_allele
        }
    
    def find_exact_match(self, marker_chrom: str, marker_pos: int, 
                        marker_ref: str, marker_alt: str, 
                        vcf_variants: Dict) -> Optional[Dict]:
        """
        Find exact match for a marker in VCF variants.
        
        Args:
            marker_chrom: Chromosome from Excel
            marker_pos: Position from Excel
            marker_ref: Reference allele from Excel
            marker_alt: Variant allele from Excel
            vcf_variants: Dictionary of VCF variants
            
        Returns:
            Variant data if exact match found, None otherwise
        """
        variant_key = (str(marker_chrom), int(marker_pos))
        
        if variant_key in vcf_variants:
            variant = vcf_variants[variant_key]
            # Check if alleles also match (optional strict matching)
            return variant
            
        return None
    
    def find_proximity_matches(self, marker_chrom: str, marker_pos: int, 
                             vcf_variants: Dict, window: int = 100) -> List[Dict]:
        """
        Find variants within proximity window (±100 bp by default).
        
        Args:
            marker_chrom: Chromosome from Excel
            marker_pos: Position from Excel
            vcf_variants: Dictionary of VCF variants
            window: Base pair window for proximity search
            
        Returns:
            List of nearby variants
        """
        proximity_matches = []
        
        for (v_chrom, v_pos), variant in vcf_variants.items():
            if (str(v_chrom) == str(marker_chrom) and 
                abs(v_pos - marker_pos) <= window and
                v_pos != marker_pos):  # Exclude exact position matches
                
                proximity_matches.append({
                    'variant': variant,
                    'distance': abs(v_pos - marker_pos)
                })
        
        # Sort by distance
        proximity_matches.sort(key=lambda x: x['distance'])
        return proximity_matches
    
    def process_markers_against_vcf(self, sheet_name: str, markers_df: pd.DataFrame) -> Dict:
        """
        Process all markers from a sheet against all VCF files.
        
        Args:
            sheet_name: Name of the Excel sheet
            markers_df: DataFrame containing marker data
            
        Returns:
            Dictionary with analysis results
        """
        logger.info(f"Processing {len(markers_df)} markers from sheet '{sheet_name}'...")
        
        results = {}
        
        # Process each marker
        for idx, marker in tqdm(markers_df.iterrows(), total=len(markers_df), 
                              desc=f"Processing {sheet_name} markers"):
            
            marker_id = marker.get('IDs', f"marker_{idx}")
            marker_chrom = str(marker['CHROM'])
            marker_pos = int(marker['POS'])
            marker_ref = marker.get('Reference Allele', '')
            marker_alt = marker.get('Variant Allele', '')
            
            # Initialize result structure for this marker
            marker_result = {
                'Marker_ID': marker_id,
                'Excel_CHROM': marker_chrom,
                'Excel_POS': marker_pos,
                'Excel_Ref': marker_ref,
                'Excel_Var': marker_alt,
                'Strand': marker.get('Strand', ''),
                'List': marker.get('List', ''),
                'Match_All': True,  # Will be set to False if missing from any VCF
                'missing_from_vcfs': [],
                'present_in_vcfs': [],
                'proximity_matches': []
            }
            
            # Combined frequency data across all VCFs
            combined_freq_data = {
                'total_ref_alleles': 0,
                'total_alt_alleles': 0,
                'total_alleles': 0
            }
            
            # Process against each VCF file
            for vcf_file in self.vcf_files:
                vcf_identifier = self.vcf_identifiers[vcf_file]
                
                # Load VCF data if not cached
                if vcf_file not in self.vcf_data_cache:
                    logger.debug(f"Loading VCF file: {vcf_file.name}")
                    self.vcf_data_cache[vcf_file] = self.parse_vcf_file(vcf_file)
                
                vcf_variants = self.vcf_data_cache[vcf_file]
                
                # Look for exact match
                exact_match = self.find_exact_match(marker_chrom, marker_pos, 
                                                  marker_ref, marker_alt, vcf_variants)
                
                if exact_match:
                    # Exact match found
                    marker_result[f'{vcf_identifier}_CHROM_exact'] = exact_match['CHROM']
                    marker_result[f'{vcf_identifier}_POS_exact'] = exact_match['POS']
                    marker_result[f'{vcf_identifier}_Match'] = 'EXACT'
                    marker_result['present_in_vcfs'].append(vcf_identifier)
                    
                    # Calculate frequencies
                    freq_data = self.calculate_allele_frequencies(exact_match, vcf_identifier)
                    marker_result.update(freq_data)
                    
                    # Add to combined frequency calculation
                    combined_freq_data['total_ref_alleles'] += freq_data['homo_ref'] * 2 + freq_data['het']
                    combined_freq_data['total_alt_alleles'] += freq_data['homo_alt'] * 2 + freq_data['het']
                    combined_freq_data['total_alleles'] += freq_data['total_alleles']
                    
                else:
                    # No exact match, check proximity
                    proximity_matches = self.find_proximity_matches(marker_chrom, marker_pos, vcf_variants)
                    
                    if proximity_matches:
                        # Found proximity matches
                        closest_match = proximity_matches[0]
                        marker_result[f'{vcf_identifier}_CHROM_exact'] = closest_match['variant']['CHROM']
                        marker_result[f'{vcf_identifier}_POS_exact'] = closest_match['variant']['POS']
                        marker_result[f'{vcf_identifier}_Match'] = f"PROXIMITY_{closest_match['distance']}bp"
                        
                        # Store proximity match info
                        marker_result['proximity_matches'].append({
                            'vcf': vcf_identifier,
                            'distance': closest_match['distance'],
                            'position': closest_match['variant']['POS']
                        })
                    else:
                        # No match found
                        marker_result[f'{vcf_identifier}_CHROM_exact'] = ''
                        marker_result[f'{vcf_identifier}_POS_exact'] = ''
                        marker_result[f'{vcf_identifier}_Match'] = 'NO_MATCH'
                        marker_result['Match_All'] = False
                        marker_result['missing_from_vcfs'].append(vcf_identifier)
                    
                    # Set frequency values to 0 for missing markers
                    marker_result[f'freq_ref_of_{vcf_identifier}'] = 0.0
                    marker_result[f'freq_alt_of_{vcf_identifier}'] = 0.0
            
            # Calculate combined frequencies across all VCFs
            if combined_freq_data['total_alleles'] > 0:
                marker_result['freq_ref_of_all_vcf'] = combined_freq_data['total_ref_alleles'] / combined_freq_data['total_alleles']
                marker_result['freq_alt_of_all_vcf'] = combined_freq_data['total_alt_alleles'] / combined_freq_data['total_alleles']
            else:
                marker_result['freq_ref_of_all_vcf'] = 0.0
                marker_result['freq_alt_of_all_vcf'] = 0.0
            
            results[marker_id] = marker_result
        
        return results
    
    def generate_summary_statistics(self) -> Dict:
        """
        Generate summary statistics for the analysis.
        
        Returns:
            Dictionary with summary statistics
        """
        summary = {
            'total_markers_analyzed': 0,
            'total_vcf_files': len(self.vcf_files),
            'markers_by_sheet': {},
            'match_statistics': {
                'exact_matches': 0,
                'proximity_matches': 0,
                'no_matches': 0,
                'match_all_vcfs': 0
            },
            'vcf_file_info': []
        }
        
        # VCF file information
        for vcf_file in self.vcf_files:
            identifier = self.vcf_identifiers[vcf_file]
            summary['vcf_file_info'].append({
                'filename': vcf_file.name,
                'identifier': identifier,
                'file_size_mb': vcf_file.stat().st_size / (1024 * 1024)
            })
        
        # Analyze results from each sheet
        for sheet_name, results in self.marker_results.items():
            sheet_stats = {
                'total_markers': len(results),
                'exact_matches': 0,
                'proximity_matches': 0,
                'no_matches': 0,
                'match_all_vcfs': 0
            }
            
            for marker_id, marker_data in results.items():
                summary['total_markers_analyzed'] += 1
                
                # Count match types
                match_all = marker_data.get('Match_All', False)
                if match_all:
                    sheet_stats['match_all_vcfs'] += 1
                    sheet_stats['exact_matches'] += 1
                else:
                    # Check if any proximity matches
                    has_proximity = any('PROXIMITY' in str(marker_data.get(f'{vcf_id}_Match', '')) 
                                      for vcf_id in self.vcf_identifiers.values())
                    if has_proximity:
                        sheet_stats['proximity_matches'] += 1
                    else:
                        sheet_stats['no_matches'] += 1
            
            summary['markers_by_sheet'][sheet_name] = sheet_stats
            
            # Add to overall statistics
            summary['match_statistics']['exact_matches'] += sheet_stats['exact_matches']
            summary['match_statistics']['proximity_matches'] += sheet_stats['proximity_matches']
            summary['match_statistics']['no_matches'] += sheet_stats['no_matches']
            summary['match_statistics']['match_all_vcfs'] += sheet_stats['match_all_vcfs']
        
        return summary
    
    def create_output_excel(self, summary_stats: Dict) -> None:
        """
        Create comprehensive Excel output with multiple sheets.
        
        Args:
            summary_stats: Summary statistics dictionary
        """
        logger.info("Creating output Excel file...")
        
        with pd.ExcelWriter(self.output_path, engine='openpyxl') as writer:
            
            # Create comparison sheets for each Excel input sheet
            for sheet_name, results in self.marker_results.items():
                output_sheet_name = f"{sheet_name.replace(' ', '_')}_Results"
                
                # Convert results to DataFrame
                results_list = []
                for marker_id, marker_data in results.items():
                    row = marker_data.copy()
                    
                    # Convert lists to strings for Excel compatibility
                    if 'missing_from_vcfs' in row:
                        row['missing_from_vcfs'] = ','.join(row['missing_from_vcfs'])
                    if 'present_in_vcfs' in row:
                        row['present_in_vcfs'] = ','.join(row['present_in_vcfs'])
                    if 'proximity_matches' in row:
                        row['proximity_matches'] = str(row['proximity_matches'])
                    
                    results_list.append(row)
                
                results_df = pd.DataFrame(results_list)
                results_df.to_excel(writer, sheet_name=output_sheet_name, index=False)
                logger.info(f"Created sheet: {output_sheet_name} with {len(results_df)} markers")
            
            # Create summary sheet
            summary_data = []
            summary_data.append(['Analysis Summary', ''])
            summary_data.append(['Total Markers Analyzed', summary_stats['total_markers_analyzed']])
            summary_data.append(['Total VCF Files', summary_stats['total_vcf_files']])
            summary_data.append(['', ''])
            summary_data.append(['Match Statistics', ''])
            summary_data.append(['Exact Matches', summary_stats['match_statistics']['exact_matches']])
            summary_data.append(['Proximity Matches', summary_stats['match_statistics']['proximity_matches']])
            summary_data.append(['No Matches', summary_stats['match_statistics']['no_matches']])
            summary_data.append(['Match All VCFs', summary_stats['match_statistics']['match_all_vcfs']])
            summary_data.append(['', ''])
            summary_data.append(['Markers by Sheet', ''])
            
            for sheet_name, sheet_stats in summary_stats['markers_by_sheet'].items():
                summary_data.append([f'{sheet_name} - Total', sheet_stats['total_markers']])
                summary_data.append([f'{sheet_name} - Exact Matches', sheet_stats['exact_matches']])
                summary_data.append([f'{sheet_name} - Proximity Matches', sheet_stats['proximity_matches']])
                summary_data.append([f'{sheet_name} - No Matches', sheet_stats['no_matches']])
            
            summary_data.append(['', ''])
            summary_data.append(['VCF File Information', ''])
            summary_data.append(['Filename', 'Identifier', 'Size (MB)'])
            
            for vcf_info in summary_stats['vcf_file_info']:
                summary_data.append([vcf_info['filename'], vcf_info['identifier'], 
                                   round(vcf_info['file_size_mb'], 2)])
            
            summary_df = pd.DataFrame(summary_data, columns=['Parameter', 'Value', 'Extra'])
            summary_df.to_excel(writer, sheet_name='Summary', index=False)
            
            # Create missing markers sheet
            missing_markers = []
            for sheet_name, results in self.marker_results.items():
                for marker_id, marker_data in results.items():
                    if not marker_data.get('Match_All', True):
                        missing_markers.append({
                            'Sheet': sheet_name,
                            'Marker_ID': marker_id,
                            'Excel_CHROM': marker_data['Excel_CHROM'],
                            'Excel_POS': marker_data['Excel_POS'],
                            'Missing_from_VCFs': ','.join(marker_data.get('missing_from_vcfs', [])),
                            'Present_in_VCFs': ','.join(marker_data.get('present_in_vcfs', [])),
                            'Proximity_matches': str(marker_data.get('proximity_matches', []))
                        })
            
            if missing_markers:
                missing_df = pd.DataFrame(missing_markers)
                missing_df.to_excel(writer, sheet_name='Missing_markers', index=False)
                logger.info(f"Created Missing_markers sheet with {len(missing_markers)} entries")
            
            # Create allele representation sheet
            allele_analysis = []
            for sheet_name, results in self.marker_results.items():
                for marker_id, marker_data in results.items():
                    ref_freq = marker_data.get('freq_ref_of_all_vcf', 0)
                    alt_freq = marker_data.get('freq_alt_of_all_vcf', 0)
                    
                    # Flag unusual frequency patterns
                    unusual_pattern = ''
                    if ref_freq == 0 and alt_freq == 0:
                        unusual_pattern = 'NO_DATA'
                    elif ref_freq < 0.01 or alt_freq < 0.01:
                        unusual_pattern = 'RARE_ALLELE'
                    elif abs(ref_freq - alt_freq) < 0.1:
                        unusual_pattern = 'BALANCED'
                    
                    allele_analysis.append({
                        'Sheet': sheet_name,
                        'Marker_ID': marker_id,
                        'Excel_CHROM': marker_data['Excel_CHROM'],
                        'Excel_POS': marker_data['Excel_POS'],
                        'REF_Frequency': ref_freq,
                        'ALT_Frequency': alt_freq,
                        'Pattern_Flag': unusual_pattern,
                        'Total_VCFs_Present': len(marker_data.get('present_in_vcfs', []))
                    })
            
            allele_df = pd.DataFrame(allele_analysis)
            allele_df.to_excel(writer, sheet_name='Allele_representation', index=False)
            logger.info(f"Created Allele_representation sheet with {len(allele_analysis)} entries")
        
        logger.info(f"Output Excel file created: {self.output_path}")
    
    def run_analysis(self) -> None:
        """
        Execute the complete genetic marker analysis pipeline.
        """
        start_time = time.time()
        logger.info("Starting genetic marker analysis...")
        
        try:
            # Step 1: Validate inputs
            if not self.validate_inputs():
                raise ValueError("Input validation failed")
            
            # Step 2: Read Excel data
            self.read_excel_data()
            
            # Step 3: Discover VCF files
            self.discover_vcf_files()
            
            # Step 4: Process markers against VCF files
            for sheet_name, markers_df in self.excel_data.items():
                logger.info(f"Processing sheet: {sheet_name}")
                results = self.process_markers_against_vcf(sheet_name, markers_df)
                self.marker_results[sheet_name] = results
            
            # Step 5: Generate summary statistics
            summary_stats = self.generate_summary_statistics()
            
            # Step 6: Create output Excel file
            self.create_output_excel(summary_stats)
            
            # Final statistics
            end_time = time.time()
            duration = end_time - start_time
            
            logger.info("Analysis completed successfully!")
            logger.info(f"Total runtime: {duration:.2f} seconds")
            logger.info(f"Processed {summary_stats['total_markers_analyzed']} markers")
            logger.info(f"Analyzed {summary_stats['total_vcf_files']} VCF files")
            logger.info(f"Exact matches: {summary_stats['match_statistics']['exact_matches']}")
            logger.info(f"Proximity matches: {summary_stats['match_statistics']['proximity_matches']}")
            logger.info(f"No matches: {summary_stats['match_statistics']['no_matches']}")
            
        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            raise


def main():
    """
    Main function to handle command-line interface.
    """
    parser = argparse.ArgumentParser(
        description="Genetic Marker Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python genetic_marker_analyzer.py --excel_file input.xlsx --vcf_folder /path/to/vcf/files --output_file output.xlsx
  python genetic_marker_analyzer.py --excel_file markers.xlsx --vcf_folder ./vcf_data --output_file results.xlsx --verbose
        """
    )
    
    parser.add_argument(
        '--excel_file',
        required=True,
        help='Path to Excel workbook with genetic marker data'
    )
    
    parser.add_argument(
        '--vcf_folder',
        required=True,
        help='Path to folder containing VCF files'
    )
    
    parser.add_argument(
        '--output_file',
        required=True,
        help='Path for output Excel file'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    parser.add_argument(
        '--log_file',
        default='genetic_marker_analysis.log',
        help='Path to log file (default: genetic_marker_analysis.log)'
    )
    
    args = parser.parse_args()
    
    # Configure logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Update log file path
    for handler in logging.getLogger().handlers:
        if isinstance(handler, logging.FileHandler):
            handler.close()
            logging.getLogger().removeHandler(handler)
    
    # Add new file handler
    file_handler = logging.FileHandler(args.log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logging.getLogger().addHandler(file_handler)
    
    try:
        # Create analyzer and run analysis
        analyzer = GeneticMarkerAnalyzer(
            excel_path=args.excel_file,
            vcf_folder=args.vcf_folder,
            output_path=args.output_file
        )
        
        analyzer.run_analysis()
        
        print(f"\nAnalysis completed successfully!")
        print(f"Results saved to: {args.output_file}")
        print(f"Log file: {args.log_file}")
        
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()