"""
Core analysis logic for combining genetic markers with VCF genotype data.

This module provides the main analysis functionality including proximity matching,
allele frequency calculations, and comprehensive genetic marker analysis.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Set, Any
from dataclasses import dataclass, field
from collections import defaultdict
import concurrent.futures
from threading import Lock

from .config import get_config, get_logger, Constants
from .vcf_processor import VCFProcessor, VCFRecord, GenotypeStats
from .excel_processor import ExcelProcessor, MarkerInfo


@dataclass
class MarkerMatch:
    """Represents a match between a marker and VCF variant."""
    marker: MarkerInfo
    vcf_record: VCFRecord
    distance: int
    match_type: str  # 'exact', 'position', 'proximity'
    allele_match: bool
    genotype_stats: GenotypeStats


@dataclass
class AnalysisResult:
    """Complete analysis result for a marker."""
    marker: MarkerInfo
    matches: List[MarkerMatch]
    best_match: Optional[MarkerMatch]
    total_variants_nearby: int
    analysis_status: str
    error_message: Optional[str] = None


@dataclass
class BatchAnalysisResults:
    """Results from batch analysis of multiple markers."""
    results: List[AnalysisResult] = field(default_factory=list)
    summary_stats: Dict[str, Any] = field(default_factory=dict)
    processing_time: float = 0.0
    total_markers: int = 0
    successful_matches: int = 0
    failed_analyses: int = 0


class ProximityMatcher:
    """Handles proximity-based matching between markers and VCF variants."""
    
    def __init__(self, proximity_bp: int = 100):
        """
        Initialize proximity matcher.
        
        Args:
            proximity_bp: Maximum distance in base pairs for proximity matching
        """
        self.proximity_bp = proximity_bp
        self.logger = get_logger(__name__)
    
    def find_nearby_variants(self, marker: MarkerInfo, vcf_records: List[VCFRecord]) -> List[Tuple[VCFRecord, int]]:
        """
        Find VCF variants near a marker position.
        
        Args:
            marker: MarkerInfo object
            vcf_records: List of VCF records to search
            
        Returns:
            List of tuples (VCFRecord, distance)
        """
        nearby_variants = []
        
        for record in vcf_records:
            if record.chrom == marker.chromosome:
                distance = abs(record.pos - marker.position)
                if distance <= self.proximity_bp:
                    nearby_variants.append((record, distance))
        
        # Sort by distance
        nearby_variants.sort(key=lambda x: x[1])
        return nearby_variants
    
    def determine_match_type(self, marker: MarkerInfo, vcf_record: VCFRecord, distance: int) -> str:
        """
        Determine the type of match between marker and VCF record.
        
        Args:
            marker: MarkerInfo object
            vcf_record: VCF record
            distance: Distance in base pairs
            
        Returns:
            Match type string
        """
        if distance == 0:
            return 'exact'
        elif distance <= 10:  # Very close
            return 'position'
        else:
            return 'proximity'
    
    def check_allele_match(self, marker: MarkerInfo, vcf_record: VCFRecord) -> bool:
        """
        Check if marker alleles match VCF record alleles.
        
        Args:
            marker: MarkerInfo object
            vcf_record: VCF record
            
        Returns:
            True if alleles match
        """
        if not marker.ref_allele or not marker.alt_allele:
            return False
        
        # Check direct match
        if (marker.ref_allele == vcf_record.ref and marker.alt_allele == vcf_record.alt):
            return True
        
        # Check reverse match
        if (marker.ref_allele == vcf_record.alt and marker.alt_allele == vcf_record.ref):
            return True
        
        # Handle multi-allelic sites (basic check)
        vcf_alts = vcf_record.alt.split(',')
        if marker.alt_allele in vcf_alts and marker.ref_allele == vcf_record.ref:
            return True
        
        return False


class MarkerAnalyzer:
    """Main analyzer for genetic markers and VCF data."""
    
    def __init__(self, vcf_processor: Optional[VCFProcessor] = None,
                 excel_processor: Optional[ExcelProcessor] = None):
        """
        Initialize marker analyzer.
        
        Args:
            vcf_processor: VCF processor instance
            excel_processor: Excel processor instance
        """
        self.config = get_config()
        self.logger = get_logger(__name__)
        
        self.vcf_processor = vcf_processor or VCFProcessor()
        self.excel_processor = excel_processor or ExcelProcessor()
        self.proximity_matcher = ProximityMatcher(self.config.proximity_bp)
        
        # Thread safety for concurrent processing
        self._lock = Lock()
        self._progress_callback = None
    
    def set_progress_callback(self, callback):
        """Set callback function for progress updates."""
        self._progress_callback = callback
    
    def _update_progress(self, current: int, total: int, message: str = ""):
        """Update progress if callback is set."""
        if self._progress_callback:
            self._progress_callback(current, total, message)
    
    def analyze_single_marker(self, marker: MarkerInfo, vcf_filepath: str) -> AnalysisResult:
        """
        Analyze a single marker against VCF data.
        
        Args:
            marker: MarkerInfo object to analyze
            vcf_filepath: Path to VCF file
            
        Returns:
            AnalysisResult object
        """
        try:
            self.logger.debug(f"Analyzing marker: {marker.name}")
            
            # Get VCF records for the marker's chromosome region
            start_pos = max(1, marker.position - self.config.proximity_bp)
            end_pos = marker.position + self.config.proximity_bp
            
            # Extract genotypes for the region
            region_genotypes = self.vcf_processor.extract_region_genotypes(
                vcf_filepath, marker.chromosome, start_pos, end_pos
            )
            
            if not region_genotypes:
                return AnalysisResult(
                    marker=marker,
                    matches=[],
                    best_match=None,
                    total_variants_nearby=0,
                    analysis_status='no_variants',
                    error_message="No variants found in proximity"
                )
            
            # Convert to VCF records for proximity matching
            vcf_records = []
            for record in self.vcf_processor.read_vcf_records(vcf_filepath):
                if (record.chrom == marker.chromosome and 
                    start_pos <= record.pos <= end_pos):
                    vcf_records.append(record)
            
            # Find nearby variants
            nearby_variants = self.proximity_matcher.find_nearby_variants(marker, vcf_records)
            
            if not nearby_variants:
                return AnalysisResult(
                    marker=marker,
                    matches=[],
                    best_match=None,
                    total_variants_nearby=0,
                    analysis_status='no_matches',
                    error_message="No variants found within proximity range"
                )
            
            # Create matches
            matches = []
            for vcf_record, distance in nearby_variants:
                # Calculate genotype statistics
                genotype_stats = self.vcf_processor._calculate_genotype_stats(
                    vcf_record.genotypes, vcf_record.ref, vcf_record.alt
                )
                
                # Determine match characteristics
                match_type = self.proximity_matcher.determine_match_type(marker, vcf_record, distance)
                allele_match = self.proximity_matcher.check_allele_match(marker, vcf_record)
                
                match = MarkerMatch(
                    marker=marker,
                    vcf_record=vcf_record,
                    distance=distance,
                    match_type=match_type,
                    allele_match=allele_match,
                    genotype_stats=genotype_stats
                )
                
                matches.append(match)
            
            # Determine best match (closest with highest call rate)
            best_match = self._select_best_match(matches)
            
            analysis_status = 'success'
            if best_match and best_match.distance == 0:
                analysis_status = 'exact_match'
            elif best_match and best_match.allele_match:
                analysis_status = 'allele_match'
            
            return AnalysisResult(
                marker=marker,
                matches=matches,
                best_match=best_match,
                total_variants_nearby=len(nearby_variants),
                analysis_status=analysis_status
            )
            
        except Exception as e:
            self.logger.error(f"Analysis failed for marker {marker.name}: {e}")
            return AnalysisResult(
                marker=marker,
                matches=[],
                best_match=None,
                total_variants_nearby=0,
                analysis_status='error',
                error_message=str(e)
            )
    
    def _select_best_match(self, matches: List[MarkerMatch]) -> Optional[MarkerMatch]:
        """
        Select the best match from a list of matches.
        
        Args:
            matches: List of MarkerMatch objects
            
        Returns:
            Best MarkerMatch or None
        """
        if not matches:
            return None
        
        # Priority scoring:
        # 1. Exact position match (distance = 0)
        # 2. Allele match
        # 3. Shortest distance
        # 4. Highest call rate
        # 5. Highest MAF (within reasonable range)
        
        def score_match(match: MarkerMatch) -> Tuple[int, int, int, float, float]:
            exact_pos = 1 if match.distance == 0 else 0
            allele_match = 1 if match.allele_match else 0
            distance_score = -match.distance  # Negative for sorting (closer is better)
            call_rate = match.genotype_stats.call_rate
            maf = match.genotype_stats.maf
            
            return (exact_pos, allele_match, distance_score, call_rate, maf)
        
        # Sort by scoring criteria (descending order)
        sorted_matches = sorted(matches, key=score_match, reverse=True)
        return sorted_matches[0]
    
    def analyze_markers_batch(self, markers: List[MarkerInfo], vcf_filepath: str,
                            use_parallel: bool = True) -> BatchAnalysisResults:
        """
        Analyze multiple markers in batch.
        
        Args:
            markers: List of MarkerInfo objects
            vcf_filepath: Path to VCF file
            use_parallel: Whether to use parallel processing
            
        Returns:
            BatchAnalysisResults object
        """
        import time
        start_time = time.time()
        
        self.logger.info(f"Starting batch analysis of {len(markers)} markers")
        
        results = []
        successful_matches = 0
        failed_analyses = 0
        
        if use_parallel and self.config.max_workers > 1:
            # Parallel processing
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
                # Submit all tasks
                future_to_marker = {
                    executor.submit(self.analyze_single_marker, marker, vcf_filepath): marker
                    for marker in markers
                }
                
                # Collect results
                for i, future in enumerate(concurrent.futures.as_completed(future_to_marker)):
                    try:
                        result = future.result()
                        results.append(result)
                        
                        if result.analysis_status in ['success', 'exact_match', 'allele_match']:
                            successful_matches += 1
                        else:
                            failed_analyses += 1
                        
                        # Update progress
                        self._update_progress(i + 1, len(markers), f"Analyzed {result.marker.name}")
                        
                    except Exception as e:
                        marker = future_to_marker[future]
                        self.logger.error(f"Parallel analysis failed for {marker.name}: {e}")
                        failed_analyses += 1
        else:
            # Sequential processing
            for i, marker in enumerate(markers):
                result = self.analyze_single_marker(marker, vcf_filepath)
                results.append(result)
                
                if result.analysis_status in ['success', 'exact_match', 'allele_match']:
                    successful_matches += 1
                else:
                    failed_analyses += 1
                
                # Update progress
                self._update_progress(i + 1, len(markers), f"Analyzed {marker.name}")
        
        processing_time = time.time() - start_time
        
        # Calculate summary statistics
        summary_stats = self._calculate_summary_stats(results)
        
        batch_results = BatchAnalysisResults(
            results=results,
            summary_stats=summary_stats,
            processing_time=processing_time,
            total_markers=len(markers),
            successful_matches=successful_matches,
            failed_analyses=failed_analyses
        )
        
        self.logger.info(f"Batch analysis complete: {successful_matches}/{len(markers)} successful matches "
                        f"in {processing_time:.2f} seconds")
        
        return batch_results
    
    def _calculate_summary_stats(self, results: List[AnalysisResult]) -> Dict[str, Any]:
        """
        Calculate summary statistics from analysis results.
        
        Args:
            results: List of AnalysisResult objects
            
        Returns:
            Dictionary of summary statistics
        """
        stats = {
            'total_markers': len(results),
            'exact_matches': 0,
            'allele_matches': 0,
            'proximity_matches': 0,
            'no_matches': 0,
            'errors': 0,
            'avg_distance': 0.0,
            'avg_call_rate': 0.0,
            'avg_maf': 0.0,
            'chromosomes': set(),
            'match_types': defaultdict(int),
            'distance_distribution': defaultdict(int)
        }
        
        successful_results = []
        total_distance = 0
        total_call_rate = 0
        total_maf = 0
        
        for result in results:
            stats['chromosomes'].add(result.marker.chromosome)
            
            if result.analysis_status == 'exact_match':
                stats['exact_matches'] += 1
            elif result.analysis_status == 'allele_match':
                stats['allele_matches'] += 1
            elif result.analysis_status == 'success':
                stats['proximity_matches'] += 1
            elif result.analysis_status == 'error':
                stats['errors'] += 1
            else:
                stats['no_matches'] += 1
            
            if result.best_match:
                successful_results.append(result)
                total_distance += result.best_match.distance
                total_call_rate += result.best_match.genotype_stats.call_rate
                total_maf += result.best_match.genotype_stats.maf
                
                stats['match_types'][result.best_match.match_type] += 1
                
                # Distance distribution (binned)
                distance_bin = (result.best_match.distance // 10) * 10
                stats['distance_distribution'][distance_bin] += 1
        
        # Calculate averages
        if successful_results:
            stats['avg_distance'] = total_distance / len(successful_results)
            stats['avg_call_rate'] = total_call_rate / len(successful_results)
            stats['avg_maf'] = total_maf / len(successful_results)
        
        # Convert sets to lists for JSON serialization
        stats['chromosomes'] = sorted(list(stats['chromosomes']))
        stats['match_types'] = dict(stats['match_types'])
        stats['distance_distribution'] = dict(stats['distance_distribution'])
        
        return stats
    
    def filter_results_by_quality(self, results: List[AnalysisResult],
                                min_call_rate: Optional[float] = None,
                                min_maf: Optional[float] = None,
                                max_distance: Optional[int] = None) -> List[AnalysisResult]:
        """
        Filter analysis results by quality criteria.
        
        Args:
            results: List of AnalysisResult objects
            min_call_rate: Minimum call rate threshold
            min_maf: Minimum MAF threshold
            max_distance: Maximum distance threshold
            
        Returns:
            Filtered list of results
        """
        if min_call_rate is None:
            min_call_rate = self.config.min_call_rate
        if min_maf is None:
            min_maf = self.config.min_maf
        if max_distance is None:
            max_distance = self.config.proximity_bp
        
        filtered_results = []
        
        for result in results:
            if result.best_match is None:
                continue
            
            # Check quality criteria
            meets_call_rate = result.best_match.genotype_stats.call_rate >= min_call_rate
            meets_maf = result.best_match.genotype_stats.maf >= min_maf
            meets_distance = result.best_match.distance <= max_distance
            
            if meets_call_rate and meets_maf and meets_distance:
                filtered_results.append(result)
        
        self.logger.info(f"Quality filter: {len(filtered_results)}/{len(results)} results passed")
        return filtered_results
    
    def get_genotype_matrix(self, results: List[AnalysisResult], 
                          samples: Optional[List[str]] = None) -> Tuple[np.ndarray, List[str], List[str]]:
        """
        Create genotype matrix from analysis results.
        
        Args:
            results: List of AnalysisResult objects
            samples: List of sample names to include (None for all)
            
        Returns:
            Tuple of (genotype_matrix, marker_names, sample_names)
        """
        # Collect all samples if not specified
        if samples is None:
            all_samples = set()
            for result in results:
                if result.best_match:
                    all_samples.update(result.best_match.vcf_record.genotypes.keys())
            samples = sorted(list(all_samples))
        
        # Filter results with successful matches
        valid_results = [r for r in results if r.best_match is not None]
        marker_names = [r.marker.name for r in valid_results]
        
        if not valid_results or not samples:
            return np.array([]), [], []
        
        # Create genotype matrix
        # Encoding: 0=homozygous ref, 1=heterozygous, 2=homozygous alt, -1=missing
        matrix = np.full((len(valid_results), len(samples)), -1, dtype=int)
        
        for i, result in enumerate(valid_results):
            genotypes = result.best_match.vcf_record.genotypes
            
            for j, sample in enumerate(samples):
                if sample in genotypes:
                    gt = genotypes[sample]
                    if gt == Constants.HOMOZYGOUS_REF:
                        matrix[i, j] = 0
                    elif gt == Constants.HETEROZYGOUS:
                        matrix[i, j] = 1
                    elif gt == Constants.HOMOZYGOUS_ALT:
                        matrix[i, j] = 2
                    # else: remains -1 (missing)
        
        self.logger.info(f"Created genotype matrix: {len(marker_names)} markers × {len(samples)} samples")
        return matrix, marker_names, samples
    
    def calculate_allele_frequencies(self, results: List[AnalysisResult]) -> Dict[str, Dict[str, float]]:
        """
        Calculate allele frequencies for analyzed markers.
        
        Args:
            results: List of AnalysisResult objects
            
        Returns:
            Dictionary mapping marker names to allele frequencies
        """
        allele_frequencies = {}
        
        for result in results:
            if not result.best_match:
                continue
            
            marker_name = result.marker.name
            stats = result.best_match.genotype_stats
            
            # Calculate frequencies
            total_alleles = (stats.total_samples - stats.missing_count) * 2
            
            if total_alleles > 0:
                ref_freq = (stats.homozygous_ref_count * 2 + stats.heterozygous_count) / total_alleles
                alt_freq = (stats.homozygous_alt_count * 2 + stats.heterozygous_count) / total_alleles
                
                allele_frequencies[marker_name] = {
                    'ref_frequency': ref_freq,
                    'alt_frequency': alt_freq,
                    'maf': min(ref_freq, alt_freq),
                    'call_rate': stats.call_rate,
                    'total_samples': stats.total_samples
                }
        
        self.logger.info(f"Calculated allele frequencies for {len(allele_frequencies)} markers")
        return allele_frequencies