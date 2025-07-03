"""
Output generation and Excel report creation module.

This module handles the creation of comprehensive reports from analysis results,
including Excel files with multiple sheets, visualizations, and summary statistics.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
import json
from datetime import datetime
import warnings

from .config import get_config, get_logger, Constants
from .marker_analyzer import AnalysisResult, BatchAnalysisResults, MarkerMatch
from .excel_processor import MarkerInfo


class ReportGenerator:
    """Generates comprehensive reports from genetic marker analysis results."""
    
    def __init__(self):
        """Initialize report generator."""
        self.config = get_config()
        self.logger = get_logger(__name__)
    
    def create_detailed_report(self, batch_results: BatchAnalysisResults, 
                             output_path: str, vcf_filepath: str,
                             excel_filepath: str) -> None:
        """
        Create a detailed Excel report with multiple sheets.
        
        Args:
            batch_results: Results from batch analysis
            output_path: Output file path for the report
            vcf_filepath: Path to the VCF file used
            excel_filepath: Path to the Excel file used
        """
        self.logger.info(f"Creating detailed report: {output_path}")
        
        try:
            with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
                # Sheet 1: Summary
                self._create_summary_sheet(batch_results, writer, vcf_filepath, excel_filepath)
                
                # Sheet 2: Detailed Results
                self._create_detailed_results_sheet(batch_results, writer)
                
                # Sheet 3: Match Statistics
                self._create_match_statistics_sheet(batch_results, writer)
                
                # Sheet 4: Allele Frequencies
                self._create_allele_frequencies_sheet(batch_results, writer)
                
                # Sheet 5: Quality Metrics
                self._create_quality_metrics_sheet(batch_results, writer)
                
                # Sheet 6: Genotype Matrix (if requested and manageable size)
                if self.config.include_detailed_stats:
                    self._create_genotype_matrix_sheet(batch_results, writer)
                
                # Sheet 7: Failed Analyses
                self._create_failed_analyses_sheet(batch_results, writer)
            
            self.logger.info(f"Detailed report created successfully: {output_path}")
            
        except Exception as e:
            error_msg = f"Failed to create detailed report: {e}"
            self.logger.error(error_msg)
            raise RuntimeError(error_msg) from e
    
    def _create_summary_sheet(self, batch_results: BatchAnalysisResults, 
                            writer: pd.ExcelWriter, vcf_filepath: str, 
                            excel_filepath: str) -> None:
        """Create summary sheet with overall statistics."""
        summary_data = []
        
        # Analysis overview
        summary_data.extend([
            ['Analysis Overview', ''],
            ['Total Markers Analyzed', batch_results.total_markers],
            ['Successful Matches', batch_results.successful_matches],
            ['Failed Analyses', batch_results.failed_analyses],
            ['Processing Time (seconds)', f"{batch_results.processing_time:.2f}"],
            ['Analysis Date', datetime.now().strftime('%Y-%m-%d %H:%M:%S')],
            ['', ''],
        ])
        
        # File information
        summary_data.extend([
            ['Input Files', ''],
            ['VCF File', Path(vcf_filepath).name],
            ['Excel File', Path(excel_filepath).name],
            ['', ''],
        ])
        
        # Match type summary
        summary_data.extend([
            ['Match Types', ''],
            ['Exact Matches', batch_results.summary_stats.get('exact_matches', 0)],
            ['Allele Matches', batch_results.summary_stats.get('allele_matches', 0)],
            ['Proximity Matches', batch_results.summary_stats.get('proximity_matches', 0)],
            ['No Matches', batch_results.summary_stats.get('no_matches', 0)],
            ['Errors', batch_results.summary_stats.get('errors', 0)],
            ['', ''],
        ])
        
        # Quality statistics
        summary_data.extend([
            ['Quality Statistics', ''],
            ['Average Distance (bp)', f"{batch_results.summary_stats.get('avg_distance', 0):.1f}"],
            ['Average Call Rate', f"{batch_results.summary_stats.get('avg_call_rate', 0):.3f}"],
            ['Average MAF', f"{batch_results.summary_stats.get('avg_maf', 0):.3f}"],
            ['', ''],
        ])
        
        # Chromosome distribution
        chromosomes = batch_results.summary_stats.get('chromosomes', [])
        summary_data.extend([
            ['Chromosome Distribution', ''],
            ['Number of Chromosomes', len(chromosomes)],
            ['Chromosomes', ', '.join(map(str, chromosomes))],
        ])
        
        # Create DataFrame and write to Excel
        summary_df = pd.DataFrame(summary_data, columns=['Metric', 'Value'])
        summary_df.to_excel(writer, sheet_name='Summary', index=False)
    
    def _create_detailed_results_sheet(self, batch_results: BatchAnalysisResults, 
                                     writer: pd.ExcelWriter) -> None:
        """Create detailed results sheet with all marker matches."""
        detailed_data = []
        
        for result in batch_results.results:
            marker = result.marker
            
            if result.best_match:
                match = result.best_match
                vcf_record = match.vcf_record
                stats = match.genotype_stats
                
                row = {
                    'Marker_Name': marker.name,
                    'Marker_Chromosome': marker.chromosome,
                    'Marker_Position': marker.position,
                    'Marker_Ref_Allele': marker.ref_allele,
                    'Marker_Alt_Allele': marker.alt_allele,
                    'VCF_Chromosome': vcf_record.chrom,
                    'VCF_Position': vcf_record.pos,
                    'VCF_ID': vcf_record.id,
                    'VCF_Ref_Allele': vcf_record.ref,
                    'VCF_Alt_Allele': vcf_record.alt,
                    'Distance_bp': match.distance,
                    'Match_Type': match.match_type,
                    'Allele_Match': match.allele_match,
                    'Call_Rate': stats.call_rate,
                    'MAF': stats.maf,
                    'Total_Samples': stats.total_samples,
                    'Missing_Count': stats.missing_count,
                    'Homozygous_Ref': stats.homozygous_ref_count,
                    'Heterozygous': stats.heterozygous_count,
                    'Homozygous_Alt': stats.homozygous_alt_count,
                    'Analysis_Status': result.analysis_status,
                    'Total_Variants_Nearby': result.total_variants_nearby
                }
            else:
                row = {
                    'Marker_Name': marker.name,
                    'Marker_Chromosome': marker.chromosome,
                    'Marker_Position': marker.position,
                    'Marker_Ref_Allele': marker.ref_allele,
                    'Marker_Alt_Allele': marker.alt_allele,
                    'VCF_Chromosome': '',
                    'VCF_Position': '',
                    'VCF_ID': '',
                    'VCF_Ref_Allele': '',
                    'VCF_Alt_Allele': '',
                    'Distance_bp': '',
                    'Match_Type': '',
                    'Allele_Match': '',
                    'Call_Rate': '',
                    'MAF': '',
                    'Total_Samples': '',
                    'Missing_Count': '',
                    'Homozygous_Ref': '',
                    'Heterozygous': '',
                    'Homozygous_Alt': '',
                    'Analysis_Status': result.analysis_status,
                    'Total_Variants_Nearby': result.total_variants_nearby
                }
            
            detailed_data.append(row)
        
        detailed_df = pd.DataFrame(detailed_data)
        detailed_df.to_excel(writer, sheet_name='Detailed_Results', index=False)
    
    def _create_match_statistics_sheet(self, batch_results: BatchAnalysisResults,
                                     writer: pd.ExcelWriter) -> None:
        """Create match statistics sheet."""
        # Distance distribution
        distance_dist = batch_results.summary_stats.get('distance_distribution', {})
        distance_data = [
            {'Distance_Range_bp': f"{k}-{k+9}" if k < 100 else f"{k}+", 'Count': v}
            for k, v in sorted(distance_dist.items())
        ]
        
        if distance_data:
            distance_df = pd.DataFrame(distance_data)
            distance_df.to_excel(writer, sheet_name='Match_Statistics', 
                                startrow=0, index=False)
        
        # Match type distribution
        match_types = batch_results.summary_stats.get('match_types', {})
        match_data = [
            {'Match_Type': k, 'Count': v}
            for k, v in match_types.items()
        ]
        
        if match_data:
            match_df = pd.DataFrame(match_data)
            start_row = len(distance_data) + 3 if distance_data else 0
            match_df.to_excel(writer, sheet_name='Match_Statistics', 
                             startrow=start_row, index=False)
        
        # Chromosome distribution
        chrom_counts = {}
        for result in batch_results.results:
            chrom = result.marker.chromosome
            chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        
        chrom_data = [
            {'Chromosome': k, 'Marker_Count': v}
            for k, v in sorted(chrom_counts.items())
        ]
        
        if chrom_data:
            chrom_df = pd.DataFrame(chrom_data)
            start_row = len(distance_data) + len(match_data) + 6 if distance_data or match_data else 0
            chrom_df.to_excel(writer, sheet_name='Match_Statistics', 
                             startrow=start_row, index=False)
    
    def _create_allele_frequencies_sheet(self, batch_results: BatchAnalysisResults,
                                       writer: pd.ExcelWriter) -> None:
        """Create allele frequencies sheet."""
        frequency_data = []
        
        for result in batch_results.results:
            if result.best_match:
                marker = result.marker
                stats = result.best_match.genotype_stats
                
                # Calculate allele frequencies
                total_called = stats.total_samples - stats.missing_count
                total_alleles = total_called * 2
                
                if total_alleles > 0:
                    ref_freq = (stats.homozygous_ref_count * 2 + stats.heterozygous_count) / total_alleles
                    alt_freq = (stats.homozygous_alt_count * 2 + stats.heterozygous_count) / total_alleles
                    
                    row = {
                        'Marker_Name': marker.name,
                        'Chromosome': marker.chromosome,
                        'Position': marker.position,
                        'Ref_Allele': result.best_match.vcf_record.ref,
                        'Alt_Allele': result.best_match.vcf_record.alt,
                        'Ref_Frequency': ref_freq,
                        'Alt_Frequency': alt_freq,
                        'MAF': min(ref_freq, alt_freq),
                        'Call_Rate': stats.call_rate,
                        'Total_Samples': stats.total_samples,
                        'Called_Samples': total_called,
                        'HWE_Expected_Het': 2 * ref_freq * alt_freq if ref_freq > 0 and alt_freq > 0 else 0,
                        'Observed_Het_Rate': stats.heterozygous_count / total_called if total_called > 0 else 0
                    }
                    
                    frequency_data.append(row)
        
        if frequency_data:
            frequency_df = pd.DataFrame(frequency_data)
            frequency_df.to_excel(writer, sheet_name='Allele_Frequencies', index=False)
    
    def _create_quality_metrics_sheet(self, batch_results: BatchAnalysisResults,
                                    writer: pd.ExcelWriter) -> None:
        """Create quality metrics sheet."""
        quality_data = []
        
        for result in batch_results.results:
            if result.best_match:
                marker = result.marker
                match = result.best_match
                stats = match.genotype_stats
                
                # Quality flags
                high_call_rate = stats.call_rate >= self.config.min_call_rate
                sufficient_maf = stats.maf >= self.config.min_maf
                close_proximity = match.distance <= self.config.proximity_bp
                
                row = {
                    'Marker_Name': marker.name,
                    'Distance_bp': match.distance,
                    'Call_Rate': stats.call_rate,
                    'MAF': stats.maf,
                    'Missing_Rate': stats.missing_count / stats.total_samples if stats.total_samples > 0 else 1,
                    'High_Call_Rate': high_call_rate,
                    'Sufficient_MAF': sufficient_maf,
                    'Close_Proximity': close_proximity,
                    'All_QC_Pass': high_call_rate and sufficient_maf and close_proximity,
                    'Match_Type': match.match_type,
                    'Allele_Match': match.allele_match,
                    'Analysis_Status': result.analysis_status
                }
                
                quality_data.append(row)
        
        if quality_data:
            quality_df = pd.DataFrame(quality_data)
            quality_df.to_excel(writer, sheet_name='Quality_Metrics', index=False)
    
    def _create_genotype_matrix_sheet(self, batch_results: BatchAnalysisResults,
                                    writer: pd.ExcelWriter) -> None:
        """Create genotype matrix sheet (if data is manageable size)."""
        try:
            # Get all samples
            all_samples = set()
            valid_results = []
            
            for result in batch_results.results:
                if result.best_match:
                    all_samples.update(result.best_match.vcf_record.genotypes.keys())
                    valid_results.append(result)
            
            samples = sorted(list(all_samples))
            
            # Check if matrix size is reasonable for Excel
            max_markers = 1000  # Limit for Excel export
            max_samples = 100   # Limit for Excel export
            
            if len(valid_results) > max_markers:
                self.logger.warning(f"Too many markers ({len(valid_results)}) for genotype matrix export. "
                                  f"Limiting to first {max_markers}.")
                valid_results = valid_results[:max_markers]
            
            if len(samples) > max_samples:
                self.logger.warning(f"Too many samples ({len(samples)}) for genotype matrix export. "
                                  f"Limiting to first {max_samples}.")
                samples = samples[:max_samples]
            
            if not valid_results or not samples:
                return
            
            # Create genotype matrix
            matrix_data = []
            
            for result in valid_results:
                row = {'Marker_Name': result.marker.name}
                genotypes = result.best_match.vcf_record.genotypes
                
                for sample in samples:
                    gt = genotypes.get(sample, Constants.MISSING_GENOTYPE)
                    # Convert to readable format
                    if gt == Constants.HOMOZYGOUS_REF:
                        row[sample] = '0/0'
                    elif gt == Constants.HETEROZYGOUS:
                        row[sample] = '0/1'
                    elif gt == Constants.HOMOZYGOUS_ALT:
                        row[sample] = '1/1'
                    else:
                        row[sample] = './.'
                
                matrix_data.append(row)
            
            if matrix_data:
                matrix_df = pd.DataFrame(matrix_data)
                matrix_df.to_excel(writer, sheet_name='Genotype_Matrix', index=False)
                
        except Exception as e:
            self.logger.warning(f"Failed to create genotype matrix sheet: {e}")
    
    def _create_failed_analyses_sheet(self, batch_results: BatchAnalysisResults,
                                    writer: pd.ExcelWriter) -> None:
        """Create failed analyses sheet."""
        failed_data = []
        
        for result in batch_results.results:
            if result.analysis_status in ['error', 'no_variants', 'no_matches']:
                row = {
                    'Marker_Name': result.marker.name,
                    'Chromosome': result.marker.chromosome,
                    'Position': result.marker.position,
                    'Analysis_Status': result.analysis_status,
                    'Error_Message': result.error_message or '',
                    'Variants_Nearby': result.total_variants_nearby
                }
                failed_data.append(row)
        
        if failed_data:
            failed_df = pd.DataFrame(failed_data)
            failed_df.to_excel(writer, sheet_name='Failed_Analyses', index=False)
    
    def create_summary_report(self, batch_results: BatchAnalysisResults,
                            output_path: str) -> None:
        """
        Create a simple summary report.
        
        Args:
            batch_results: Results from batch analysis
            output_path: Output file path for the summary
        """
        self.logger.info(f"Creating summary report: {output_path}")
        
        try:
            # Prepare summary data
            summary_data = []
            
            for result in batch_results.results:
                if result.best_match:
                    match = result.best_match
                    row = {
                        'Marker_Name': result.marker.name,
                        'Chromosome': result.marker.chromosome,
                        'Position': result.marker.position,
                        'VCF_Position': match.vcf_record.pos,
                        'Distance_bp': match.distance,
                        'Match_Type': match.match_type,
                        'Allele_Match': match.allele_match,
                        'Call_Rate': match.genotype_stats.call_rate,
                        'MAF': match.genotype_stats.maf,
                        'Status': result.analysis_status
                    }
                else:
                    row = {
                        'Marker_Name': result.marker.name,
                        'Chromosome': result.marker.chromosome,
                        'Position': result.marker.position,
                        'VCF_Position': '',
                        'Distance_bp': '',
                        'Match_Type': '',
                        'Allele_Match': '',
                        'Call_Rate': '',
                        'MAF': '',
                        'Status': result.analysis_status
                    }
                
                summary_data.append(row)
            
            # Create DataFrame and save
            summary_df = pd.DataFrame(summary_data)
            
            if output_path.endswith('.xlsx'):
                summary_df.to_excel(output_path, index=False)
            elif output_path.endswith('.csv'):
                summary_df.to_csv(output_path, index=False)
            else:
                # Default to CSV
                summary_df.to_csv(output_path + '.csv', index=False)
            
            self.logger.info(f"Summary report created: {output_path}")
            
        except Exception as e:
            error_msg = f"Failed to create summary report: {e}"
            self.logger.error(error_msg)
            raise RuntimeError(error_msg) from e
    
    def export_to_json(self, batch_results: BatchAnalysisResults,
                      output_path: str) -> None:
        """
        Export analysis results to JSON format.
        
        Args:
            batch_results: Results from batch analysis
            output_path: Output file path for JSON
        """
        self.logger.info(f"Exporting results to JSON: {output_path}")
        
        try:
            # Convert results to JSON-serializable format
            json_data = {
                'analysis_metadata': {
                    'total_markers': batch_results.total_markers,
                    'successful_matches': batch_results.successful_matches,
                    'failed_analyses': batch_results.failed_analyses,
                    'processing_time': batch_results.processing_time,
                    'analysis_date': datetime.now().isoformat()
                },
                'summary_statistics': batch_results.summary_stats,
                'results': []
            }
            
            for result in batch_results.results:
                result_data = {
                    'marker': {
                        'name': result.marker.name,
                        'chromosome': result.marker.chromosome,
                        'position': result.marker.position,
                        'ref_allele': result.marker.ref_allele,
                        'alt_allele': result.marker.alt_allele,
                        'marker_type': result.marker.marker_type,
                        'gene': result.marker.gene,
                        'annotation': result.marker.annotation
                    },
                    'analysis_status': result.analysis_status,
                    'total_variants_nearby': result.total_variants_nearby,
                    'error_message': result.error_message
                }
                
                if result.best_match:
                    match = result.best_match
                    result_data['best_match'] = {
                        'vcf_record': {
                            'chrom': match.vcf_record.chrom,
                            'pos': match.vcf_record.pos,
                            'id': match.vcf_record.id,
                            'ref': match.vcf_record.ref,
                            'alt': match.vcf_record.alt,
                            'qual': match.vcf_record.qual,
                            'filter': match.vcf_record.filter
                        },
                        'distance': match.distance,
                        'match_type': match.match_type,
                        'allele_match': match.allele_match,
                        'genotype_stats': {
                            'total_samples': match.genotype_stats.total_samples,
                            'missing_count': match.genotype_stats.missing_count,
                            'homozygous_ref_count': match.genotype_stats.homozygous_ref_count,
                            'heterozygous_count': match.genotype_stats.heterozygous_count,
                            'homozygous_alt_count': match.genotype_stats.homozygous_alt_count,
                            'call_rate': match.genotype_stats.call_rate,
                            'maf': match.genotype_stats.maf,
                            'allele_counts': match.genotype_stats.allele_counts
                        }
                    }
                
                json_data['results'].append(result_data)
            
            # Write to JSON file
            with open(output_path, 'w') as f:
                json.dump(json_data, f, indent=2)
            
            self.logger.info(f"JSON export completed: {output_path}")
            
        except Exception as e:
            error_msg = f"Failed to export to JSON: {e}"
            self.logger.error(error_msg)
            raise RuntimeError(error_msg) from e
    
    def create_qc_report(self, batch_results: BatchAnalysisResults,
                        output_path: str) -> None:
        """
        Create a quality control report highlighting potential issues.
        
        Args:
            batch_results: Results from batch analysis
            output_path: Output file path for QC report
        """
        self.logger.info(f"Creating QC report: {output_path}")
        
        try:
            qc_issues = []
            
            for result in batch_results.results:
                issues = []
                
                if result.analysis_status == 'error':
                    issues.append(f"Analysis failed: {result.error_message}")
                elif result.analysis_status in ['no_variants', 'no_matches']:
                    issues.append(f"No suitable variants found")
                elif result.best_match:
                    match = result.best_match
                    stats = match.genotype_stats
                    
                    # Check various QC criteria
                    if stats.call_rate < self.config.min_call_rate:
                        issues.append(f"Low call rate: {stats.call_rate:.3f}")
                    
                    if stats.maf < self.config.min_maf:
                        issues.append(f"Low MAF: {stats.maf:.3f}")
                    
                    if match.distance > self.config.proximity_bp:
                        issues.append(f"Large distance: {match.distance} bp")
                    
                    if not match.allele_match and result.marker.ref_allele and result.marker.alt_allele:
                        issues.append("Allele mismatch")
                    
                    if stats.missing_count / stats.total_samples > self.config.max_missing_rate:
                        issues.append(f"High missing rate: {stats.missing_count / stats.total_samples:.3f}")
                
                if issues:
                    qc_issues.append({
                        'Marker_Name': result.marker.name,
                        'Chromosome': result.marker.chromosome,
                        'Position': result.marker.position,
                        'Issues': '; '.join(issues),
                        'Analysis_Status': result.analysis_status
                    })
            
            if qc_issues:
                qc_df = pd.DataFrame(qc_issues)
                qc_df.to_excel(output_path, index=False)
                self.logger.info(f"QC report created with {len(qc_issues)} flagged markers")
            else:
                # Create empty report indicating no issues
                empty_df = pd.DataFrame(columns=['Message'])
                empty_df.loc[0] = ['No quality control issues detected']
                empty_df.to_excel(output_path, index=False)
                self.logger.info("QC report created: No issues detected")
            
        except Exception as e:
            error_msg = f"Failed to create QC report: {e}"
            self.logger.error(error_msg)
            raise RuntimeError(error_msg) from e