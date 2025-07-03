"""
VCF file processing module with caching and compression support.

This module handles VCF file reading, parsing, and genotype data extraction
with optimizations for large files including pre-caching and memory management.
"""

import gzip
import os
import pickle
import hashlib
from typing import Dict, List, Tuple, Optional, Iterator, Set, Any
from dataclasses import dataclass
from pathlib import Path
import re
from collections import defaultdict

from .config import get_config, get_logger, Constants, validate_file_path, check_file_size


@dataclass
class VCFHeader:
    """VCF header information."""
    info_lines: List[str]
    format_lines: List[str]
    samples: List[str]
    column_headers: List[str]


@dataclass
class VCFRecord:
    """Single VCF record representation."""
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: str
    filter: str
    info: str
    format: str
    genotypes: Dict[str, str]


@dataclass
class GenotypeStats:
    """Genotype statistics for a variant."""
    total_samples: int
    missing_count: int
    homozygous_ref_count: int
    heterozygous_count: int
    homozygous_alt_count: int
    call_rate: float
    maf: float
    allele_counts: Dict[str, int]


class VCFCache:
    """Caching system for VCF data to improve performance."""
    
    def __init__(self, cache_dir: Optional[str] = None, max_size_mb: int = 100):
        """
        Initialize VCF cache.
        
        Args:
            cache_dir: Directory for cache files (default: temp directory)
            max_size_mb: Maximum cache size in MB
        """
        self.logger = get_logger(__name__)
        self.max_size_mb = max_size_mb
        
        if cache_dir is None:
            cache_dir = os.path.join(os.path.expanduser("~"), ".vcf_cache")
        
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        
        self._cache_index = self._load_cache_index()
    
    def _load_cache_index(self) -> Dict[str, Dict[str, Any]]:
        """Load cache index from disk."""
        index_file = self.cache_dir / "cache_index.pkl"
        if index_file.exists():
            try:
                with open(index_file, 'rb') as f:
                    return pickle.load(f)
            except Exception as e:
                self.logger.warning(f"Failed to load cache index: {e}")
        return {}
    
    def _save_cache_index(self) -> None:
        """Save cache index to disk."""
        index_file = self.cache_dir / "cache_index.pkl"
        try:
            with open(index_file, 'wb') as f:
                pickle.dump(self._cache_index, f)
        except Exception as e:
            self.logger.warning(f"Failed to save cache index: {e}")
    
    def _get_file_hash(self, filepath: str) -> str:
        """Generate hash for file identification."""
        file_path = Path(filepath)
        stat = file_path.stat()
        
        # Use file path, size, and modification time for hash
        content = f"{filepath}_{stat.st_size}_{stat.st_mtime}"
        return hashlib.md5(content.encode()).hexdigest()
    
    def _cleanup_cache(self) -> None:
        """Remove old cache files if cache size exceeds limit."""
        try:
            total_size = sum(
                (self.cache_dir / filename).stat().st_size 
                for filename in os.listdir(self.cache_dir)
                if filename.endswith('.pkl')
            ) / (1024 * 1024)  # Convert to MB
            
            if total_size > self.max_size_mb:
                # Sort by access time and remove oldest
                cache_files = [
                    (f, self._cache_index.get(f.stem, {}).get('last_access', 0))
                    for f in self.cache_dir.glob('*.pkl')
                    if f.name != 'cache_index.pkl'
                ]
                cache_files.sort(key=lambda x: x[1])
                
                # Remove oldest files until under limit
                for cache_file, _ in cache_files:
                    cache_file.unlink()
                    if cache_file.stem in self._cache_index:
                        del self._cache_index[cache_file.stem]
                    
                    # Recalculate size
                    remaining_size = sum(
                        f.stat().st_size for f in self.cache_dir.glob('*.pkl')
                        if f.name != 'cache_index.pkl'
                    ) / (1024 * 1024)
                    
                    if remaining_size <= self.max_size_mb * 0.8:  # 80% of limit
                        break
                
                self._save_cache_index()
                
        except Exception as e:
            self.logger.warning(f"Cache cleanup failed: {e}")
    
    def get_cached_data(self, filepath: str, data_type: str) -> Optional[Any]:
        """
        Retrieve cached data for a file.
        
        Args:
            filepath: Path to VCF file
            data_type: Type of cached data ('header', 'genotypes', 'stats')
            
        Returns:
            Cached data or None if not found
        """
        file_hash = self._get_file_hash(filepath)
        cache_key = f"{file_hash}_{data_type}"
        
        if cache_key not in self._cache_index:
            return None
        
        cache_file = self.cache_dir / f"{cache_key}.pkl"
        if not cache_file.exists():
            # Remove stale index entry
            del self._cache_index[cache_key]
            return None
        
        try:
            with open(cache_file, 'rb') as f:
                data = pickle.load(f)
            
            # Update access time
            import time
            self._cache_index[cache_key]['last_access'] = time.time()
            self._save_cache_index()
            
            self.logger.debug(f"Cache hit for {filepath} ({data_type})")
            return data
            
        except Exception as e:
            self.logger.warning(f"Failed to load cached data: {e}")
            # Remove corrupted cache file
            cache_file.unlink()
            if cache_key in self._cache_index:
                del self._cache_index[cache_key]
            return None
    
    def cache_data(self, filepath: str, data_type: str, data: Any) -> None:
        """
        Cache data for a file.
        
        Args:
            filepath: Path to VCF file
            data_type: Type of data to cache
            data: Data to cache
        """
        file_hash = self._get_file_hash(filepath)
        cache_key = f"{file_hash}_{data_type}"
        cache_file = self.cache_dir / f"{cache_key}.pkl"
        
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(data, f)
            
            # Update index
            import time
            self._cache_index[cache_key] = {
                'filepath': filepath,
                'data_type': data_type,
                'created': time.time(),
                'last_access': time.time(),
                'size': cache_file.stat().st_size
            }
            
            self._save_cache_index()
            self._cleanup_cache()
            
            self.logger.debug(f"Cached {data_type} data for {filepath}")
            
        except Exception as e:
            self.logger.warning(f"Failed to cache data: {e}")


class VCFProcessor:
    """Main VCF file processor with caching and optimization features."""
    
    def __init__(self, enable_caching: bool = True, cache_dir: Optional[str] = None):
        """
        Initialize VCF processor.
        
        Args:
            enable_caching: Whether to enable caching
            cache_dir: Directory for cache files
        """
        self.config = get_config()
        self.logger = get_logger(__name__)
        self.enable_caching = enable_caching
        
        if self.enable_caching:
            self.cache = VCFCache(cache_dir, self.config.cache_size_mb)
        else:
            self.cache = None
        
        self._header_cache = {}
        self._genotype_cache = {}
    
    def _open_vcf_file(self, filepath: str):
        """
        Open VCF file, handling both regular and gzipped files.
        
        Args:
            filepath: Path to VCF file
            
        Returns:
            File handle
        """
        path = validate_file_path(filepath)
        
        if str(path).endswith('.gz'):
            return gzip.open(path, 'rt', encoding='utf-8')
        else:
            return open(path, 'r', encoding='utf-8')
    
    def _parse_vcf_header(self, file_handle) -> VCFHeader:
        """
        Parse VCF header information.
        
        Args:
            file_handle: Open file handle
            
        Returns:
            VCFHeader object
        """
        info_lines = []
        format_lines = []
        column_headers = []
        samples = []
        
        for line in file_handle:
            line = line.strip()
            
            if line.startswith(Constants.VCF_HEADER_PREFIX):
                if line.startswith("##INFO="):
                    info_lines.append(line)
                elif line.startswith("##FORMAT="):
                    format_lines.append(line)
            elif line.startswith(Constants.VCF_COLUMN_HEADER):
                column_headers = line.split('\t')
                # Extract sample names (columns after FORMAT)
                if len(column_headers) > 9:
                    samples = column_headers[9:]
                break
        
        if not column_headers:
            raise ValueError("Invalid VCF format: missing column headers")
        
        # Validate required columns
        for required_col in Constants.VCF_REQUIRED_COLUMNS:
            if required_col not in column_headers:
                raise ValueError(f"Missing required VCF column: {required_col}")
        
        return VCFHeader(info_lines, format_lines, samples, column_headers)
    
    def get_vcf_header(self, filepath: str) -> VCFHeader:
        """
        Get VCF header information with caching.
        
        Args:
            filepath: Path to VCF file
            
        Returns:
            VCFHeader object
        """
        # Check cache first
        if self.cache and self.enable_caching:
            cached_header = self.cache.get_cached_data(filepath, 'header')
            if cached_header:
                return cached_header
        
        # Check memory cache
        if filepath in self._header_cache:
            return self._header_cache[filepath]
        
        check_file_size(filepath)
        
        try:
            with self._open_vcf_file(filepath) as f:
                header = self._parse_vcf_header(f)
            
            # Cache the result
            self._header_cache[filepath] = header
            if self.cache and self.enable_caching:
                self.cache.cache_data(filepath, 'header', header)
            
            self.logger.info(f"Parsed VCF header: {len(header.samples)} samples")
            return header
            
        except Exception as e:
            error_msg = Constants.ERROR_MESSAGES['processing_error'].format(
                filepath=filepath, error=str(e)
            )
            self.logger.error(error_msg)
            raise RuntimeError(error_msg) from e
    
    def _parse_genotype(self, gt_string: str) -> str:
        """
        Parse genotype string to standardized format.
        
        Args:
            gt_string: Raw genotype string from VCF
            
        Returns:
            Standardized genotype string
        """
        if not gt_string or gt_string.startswith('.'):
            return Constants.MISSING_GENOTYPE
        
        # Extract GT field (first field in FORMAT)
        gt_field = gt_string.split(':')[0]
        
        # Handle different separators
        if '/' in gt_field:
            alleles = gt_field.split('/')
        elif '|' in gt_field:
            alleles = gt_field.split('|')
        else:
            return Constants.MISSING_GENOTYPE
        
        try:
            # Convert to integers and sort for consistency
            allele_ints = sorted([int(a) for a in alleles if a != '.'])
            if len(allele_ints) != len(alleles):
                return Constants.MISSING_GENOTYPE
            
            return f"{allele_ints[0]}/{allele_ints[1]}" if len(allele_ints) == 2 else Constants.MISSING_GENOTYPE
            
        except (ValueError, IndexError):
            return Constants.MISSING_GENOTYPE
    
    def _calculate_genotype_stats(self, genotypes: Dict[str, str], ref: str, alt: str) -> GenotypeStats:
        """
        Calculate statistics for genotype data.
        
        Args:
            genotypes: Dictionary of sample -> genotype
            ref: Reference allele
            alt: Alternative allele
            
        Returns:
            GenotypeStats object
        """
        total_samples = len(genotypes)
        missing_count = 0
        homozygous_ref_count = 0
        heterozygous_count = 0
        homozygous_alt_count = 0
        allele_counts = defaultdict(int)
        
        for genotype in genotypes.values():
            if genotype == Constants.MISSING_GENOTYPE:
                missing_count += 1
            elif genotype == Constants.HOMOZYGOUS_REF:
                homozygous_ref_count += 1
                allele_counts[ref] += 2
            elif genotype == Constants.HETEROZYGOUS:
                heterozygous_count += 1
                allele_counts[ref] += 1
                allele_counts[alt] += 1
            elif genotype == Constants.HOMOZYGOUS_ALT:
                homozygous_alt_count += 1
                allele_counts[alt] += 2
        
        # Calculate call rate
        called_samples = total_samples - missing_count
        call_rate = called_samples / total_samples if total_samples > 0 else 0.0
        
        # Calculate minor allele frequency
        total_alleles = called_samples * 2
        maf = 0.0
        if total_alleles > 0:
            ref_freq = allele_counts[ref] / total_alleles
            alt_freq = allele_counts[alt] / total_alleles
            maf = min(ref_freq, alt_freq)
        
        return GenotypeStats(
            total_samples=total_samples,
            missing_count=missing_count,
            homozygous_ref_count=homozygous_ref_count,
            heterozygous_count=heterozygous_count,
            homozygous_alt_count=homozygous_alt_count,
            call_rate=call_rate,
            maf=maf,
            allele_counts=dict(allele_counts)
        )
    
    def read_vcf_records(self, filepath: str, chunk_size: Optional[int] = None) -> Iterator[VCFRecord]:
        """
        Read VCF records in chunks for memory efficiency.
        
        Args:
            filepath: Path to VCF file
            chunk_size: Number of records to read at once
            
        Yields:
            VCFRecord objects
        """
        if chunk_size is None:
            chunk_size = self.config.chunk_size
        
        header = self.get_vcf_header(filepath)
        
        try:
            with self._open_vcf_file(filepath) as f:
                # Skip header lines
                for line in f:
                    if line.startswith(Constants.VCF_COLUMN_HEADER):
                        break
                
                records_read = 0
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) < 9:
                        continue
                    
                    # Parse basic fields
                    chrom, pos, id_field, ref, alt, qual, filter_field, info, format_field = fields[:9]
                    
                    # Parse genotypes
                    genotypes = {}
                    if len(fields) > 9 and len(header.samples) > 0:
                        for i, sample in enumerate(header.samples):
                            if i + 9 < len(fields):
                                genotypes[sample] = self._parse_genotype(fields[i + 9])
                    
                    record = VCFRecord(
                        chrom=chrom,
                        pos=int(pos),
                        id=id_field,
                        ref=ref,
                        alt=alt,
                        qual=qual,
                        filter=filter_field,
                        info=info,
                        format=format_field,
                        genotypes=genotypes
                    )
                    
                    yield record
                    
                    records_read += 1
                    if records_read % chunk_size == 0:
                        self.logger.debug(f"Processed {records_read} VCF records")
                
                self.logger.info(f"Finished reading {records_read} VCF records from {filepath}")
                
        except Exception as e:
            error_msg = Constants.ERROR_MESSAGES['processing_error'].format(
                filepath=filepath, error=str(e)
            )
            self.logger.error(error_msg)
            raise RuntimeError(error_msg) from e
    
    def get_samples(self, filepath: str) -> List[str]:
        """
        Get list of sample names from VCF file.
        
        Args:
            filepath: Path to VCF file
            
        Returns:
            List of sample names
        """
        header = self.get_vcf_header(filepath)
        return header.samples
    
    def get_variant_count(self, filepath: str) -> int:
        """
        Get total number of variants in VCF file.
        
        Args:
            filepath: Path to VCF file
            
        Returns:
            Number of variants
        """
        count = 0
        try:
            with self._open_vcf_file(filepath) as f:
                # Skip header
                for line in f:
                    if line.startswith(Constants.VCF_COLUMN_HEADER):
                        break
                
                # Count data lines
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        count += 1
            
            self.logger.info(f"VCF file {filepath} contains {count} variants")
            return count
            
        except Exception as e:
            self.logger.error(f"Failed to count variants in {filepath}: {e}")
            return 0
    
    def extract_region_genotypes(self, filepath: str, chrom: str, start: int, end: int) -> Dict[str, Dict[str, str]]:
        """
        Extract genotypes for a specific genomic region.
        
        Args:
            filepath: Path to VCF file
            chrom: Chromosome name
            start: Start position
            end: End position
            
        Returns:
            Dictionary mapping variant IDs to sample genotypes
        """
        region_genotypes = {}
        
        for record in self.read_vcf_records(filepath):
            if (record.chrom == chrom and start <= record.pos <= end):
                variant_id = f"{record.chrom}:{record.pos}:{record.ref}:{record.alt}"
                region_genotypes[variant_id] = record.genotypes.copy()
        
        self.logger.info(f"Extracted genotypes for {len(region_genotypes)} variants in region {chrom}:{start}-{end}")
        return region_genotypes
    
    def get_variant_stats(self, filepath: str) -> Dict[str, GenotypeStats]:
        """
        Calculate statistics for all variants in VCF file.
        
        Args:
            filepath: Path to VCF file
            
        Returns:
            Dictionary mapping variant IDs to statistics
        """
        # Check cache first
        if self.cache and self.enable_caching:
            cached_stats = self.cache.get_cached_data(filepath, 'stats')
            if cached_stats:
                return cached_stats
        
        variant_stats = {}
        
        for record in self.read_vcf_records(filepath):
            variant_id = f"{record.chrom}:{record.pos}:{record.ref}:{record.alt}"
            stats = self._calculate_genotype_stats(record.genotypes, record.ref, record.alt)
            variant_stats[variant_id] = stats
        
        # Cache the results
        if self.cache and self.enable_caching:
            self.cache.cache_data(filepath, 'stats', variant_stats)
        
        self.logger.info(f"Calculated statistics for {len(variant_stats)} variants")
        return variant_stats
    
    def filter_variants(self, filepath: str, min_call_rate: Optional[float] = None, 
                       min_maf: Optional[float] = None) -> Set[str]:
        """
        Filter variants based on quality criteria.
        
        Args:
            filepath: Path to VCF file
            min_call_rate: Minimum call rate threshold
            min_maf: Minimum minor allele frequency threshold
            
        Returns:
            Set of variant IDs that pass filters
        """
        if min_call_rate is None:
            min_call_rate = self.config.min_call_rate
        if min_maf is None:
            min_maf = self.config.min_maf
        
        variant_stats = self.get_variant_stats(filepath)
        passed_variants = set()
        
        for variant_id, stats in variant_stats.items():
            if stats.call_rate >= min_call_rate and stats.maf >= min_maf:
                passed_variants.add(variant_id)
        
        self.logger.info(f"Filter results: {len(passed_variants)}/{len(variant_stats)} variants passed quality filters")
        return passed_variants