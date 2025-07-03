"""
Main orchestration script that coordinates all components of the genetic marker analysis system.

This module provides the primary interface for running complete analyses,
managing workflow, and coordinating between all system components.
"""

import os
import sys
import argparse
import time
from typing import Optional, List, Dict, Any, Callable
from pathlib import Path
import traceback

from .config import (
    get_config_manager, get_config, get_logger, 
    validate_file_path, check_file_size, Constants
)
from .vcf_processor import VCFProcessor
from .excel_processor import ExcelProcessor
from .marker_analyzer import MarkerAnalyzer
from .report_generator import ReportGenerator


class ProgressTracker:
    """Tracks and reports progress of analysis operations."""
    
    def __init__(self, total_steps: int = 100, show_progress: bool = True):
        """
        Initialize progress tracker.
        
        Args:
            total_steps: Total number of steps in the process
            show_progress: Whether to display progress messages
        """
        self.total_steps = total_steps
        self.current_step = 0
        self.show_progress = show_progress
        self.start_time = time.time()
        self.logger = get_logger(__name__)
    
    def update(self, step: int, message: str = ""):
        """
        Update progress.
        
        Args:
            step: Current step number
            message: Optional progress message
        """
        self.current_step = step
        
        if self.show_progress:
            progress_pct = (step / self.total_steps) * 100
            elapsed_time = time.time() - self.start_time
            
            if step > 0:
                estimated_total = elapsed_time * (self.total_steps / step)
                remaining_time = estimated_total - elapsed_time
                time_msg = f" (ETA: {remaining_time:.1f}s)"
            else:
                time_msg = ""
            
            progress_msg = f"Progress: {step}/{self.total_steps} ({progress_pct:.1f}%){time_msg}"
            if message:
                progress_msg += f" - {message}"
            
            print(f"\r{progress_msg}", end="", flush=True)
            
            if step == self.total_steps:
                print()  # New line when complete
    
    def finish(self, message: str = "Complete"):
        """Mark progress as finished."""
        self.update(self.total_steps, message)


class GeneticMarkerAnalysisController:
    """Main controller for the genetic marker analysis system."""
    
    def __init__(self, config_file: Optional[str] = None, 
                 enable_caching: bool = True, show_progress: bool = True):
        """
        Initialize the analysis controller.
        
        Args:
            config_file: Optional path to configuration file
            enable_caching: Whether to enable VCF caching
            show_progress: Whether to show progress messages
        """
        # Initialize configuration
        self.config_manager = get_config_manager(config_file)
        self.config = get_config()
        self.logger = get_logger(__name__)
        self.show_progress = show_progress
        
        # Initialize processors
        self.vcf_processor = VCFProcessor(enable_caching=enable_caching)
        self.excel_processor = ExcelProcessor()
        self.marker_analyzer = MarkerAnalyzer(self.vcf_processor, self.excel_processor)
        self.report_generator = ReportGenerator()
        
        # Progress tracking
        self.progress_tracker = None
        
        self.logger.info("Genetic Marker Analysis Controller initialized")
    
    def validate_inputs(self, vcf_filepath: str, excel_filepath: str,
                       output_dir: str) -> Dict[str, Any]:
        """
        Validate input files and output directory.
        
        Args:
            vcf_filepath: Path to VCF file
            excel_filepath: Path to Excel file
            output_dir: Output directory path
            
        Returns:
            Dictionary with validation results
        """
        validation_results = {
            'vcf_valid': False,
            'excel_valid': False,
            'output_dir_valid': False,
            'errors': [],
            'warnings': []
        }
        
        try:
            # Validate VCF file
            self.logger.info(f"Validating VCF file: {vcf_filepath}")
            vcf_path = validate_file_path(vcf_filepath)
            
            if not any(str(vcf_path).lower().endswith(ext) for ext in Constants.SUPPORTED_VCF_EXTENSIONS):
                validation_results['errors'].append(f"Unsupported VCF format: {vcf_filepath}")
            else:
                check_file_size(vcf_filepath)
                
                # Basic VCF header validation
                header = self.vcf_processor.get_vcf_header(vcf_filepath)
                if not header.samples:
                    validation_results['warnings'].append("No samples found in VCF file")
                
                validation_results['vcf_valid'] = True
                self.logger.info(f"VCF validation passed: {len(header.samples)} samples")
        
        except Exception as e:
            validation_results['errors'].append(f"VCF validation failed: {e}")
        
        try:
            # Validate Excel file
            self.logger.info(f"Validating Excel file: {excel_filepath}")
            excel_validation = self.excel_processor.validate_excel_file(excel_filepath)
            
            if not excel_validation.is_valid:
                validation_results['errors'].extend(excel_validation.errors)
            else:
                validation_results['excel_valid'] = True
                validation_results['warnings'].extend(excel_validation.warnings)
                self.logger.info(f"Excel validation passed: {excel_validation.marker_count} markers")
        
        except Exception as e:
            validation_results['errors'].append(f"Excel validation failed: {e}")
        
        try:
            # Validate output directory
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            
            # Test write permissions
            test_file = output_path / "test_write.tmp"
            with open(test_file, 'w') as f:
                f.write("test")
            test_file.unlink()
            
            validation_results['output_dir_valid'] = True
            self.logger.info(f"Output directory validated: {output_dir}")
        
        except Exception as e:
            validation_results['errors'].append(f"Output directory validation failed: {e}")
        
        return validation_results
    
    def run_analysis(self, vcf_filepath: str, excel_filepath: str, 
                    output_dir: str, output_prefix: str = "analysis",
                    validate_inputs: bool = True, 
                    create_detailed_report: bool = True,
                    create_summary_report: bool = True,
                    create_qc_report: bool = True,
                    export_json: bool = False,
                    use_parallel: bool = True) -> Dict[str, Any]:
        """
        Run complete genetic marker analysis.
        
        Args:
            vcf_filepath: Path to VCF file
            excel_filepath: Path to Excel file  
            output_dir: Output directory
            output_prefix: Prefix for output files
            validate_inputs: Whether to validate inputs first
            create_detailed_report: Whether to create detailed Excel report
            create_summary_report: Whether to create summary report
            create_qc_report: Whether to create QC report
            export_json: Whether to export results to JSON
            use_parallel: Whether to use parallel processing
            
        Returns:
            Dictionary with analysis results and file paths
        """
        start_time = time.time()
        self.logger.info("=== Starting Genetic Marker Analysis ===")
        
        analysis_results = {
            'success': False,
            'processing_time': 0.0,
            'output_files': [],
            'summary_stats': {},
            'errors': [],
            'warnings': []
        }
        
        try:
            # Step 1: Input validation
            if validate_inputs:
                self.logger.info("Step 1: Validating inputs...")
                if self.show_progress:
                    self.progress_tracker = ProgressTracker(6, True)
                    self.progress_tracker.update(1, "Validating inputs")
                
                validation = self.validate_inputs(vcf_filepath, excel_filepath, output_dir)
                
                if validation['errors']:
                    analysis_results['errors'] = validation['errors']
                    self.logger.error(f"Input validation failed: {validation['errors']}")
                    return analysis_results
                
                analysis_results['warnings'].extend(validation['warnings'])
            
            # Step 2: Load marker data
            self.logger.info("Step 2: Loading marker data...")
            if self.progress_tracker:
                self.progress_tracker.update(2, "Loading marker data")
            
            markers = self.excel_processor.read_markers(excel_filepath, validate=False)
            self.logger.info(f"Loaded {len(markers)} markers from Excel file")
            
            if not markers:
                analysis_results['errors'].append("No valid markers found in Excel file")
                return analysis_results
            
            # Step 3: Set up progress tracking for marker analysis
            if self.progress_tracker:
                self.progress_tracker.update(3, "Preparing analysis")
            
            def analysis_progress_callback(current: int, total: int, message: str = ""):
                if self.progress_tracker:
                    # Map analysis progress to overall progress (steps 3-4)
                    analysis_progress = 3 + (current / total)
                    self.progress_tracker.update(int(analysis_progress), f"Analyzing markers: {message}")
            
            self.marker_analyzer.set_progress_callback(analysis_progress_callback)
            
            # Step 4: Run marker analysis
            self.logger.info("Step 4: Running marker analysis...")
            batch_results = self.marker_analyzer.analyze_markers_batch(
                markers, vcf_filepath, use_parallel=use_parallel
            )
            
            self.logger.info(f"Analysis completed: {batch_results.successful_matches}/{batch_results.total_markers} successful matches")
            
            # Step 5: Generate reports
            self.logger.info("Step 5: Generating reports...")
            if self.progress_tracker:
                self.progress_tracker.update(5, "Generating reports")
            
            output_files = []
            
            # Detailed Excel report
            if create_detailed_report:
                detailed_path = Path(output_dir) / f"{output_prefix}_detailed_report.xlsx"
                self.report_generator.create_detailed_report(
                    batch_results, str(detailed_path), vcf_filepath, excel_filepath
                )
                output_files.append(str(detailed_path))
                self.logger.info(f"Created detailed report: {detailed_path}")
            
            # Summary report
            if create_summary_report:
                summary_path = Path(output_dir) / f"{output_prefix}_summary.xlsx"
                self.report_generator.create_summary_report(batch_results, str(summary_path))
                output_files.append(str(summary_path))
                self.logger.info(f"Created summary report: {summary_path}")
            
            # QC report
            if create_qc_report:
                qc_path = Path(output_dir) / f"{output_prefix}_qc_report.xlsx"
                self.report_generator.create_qc_report(batch_results, str(qc_path))
                output_files.append(str(qc_path))
                self.logger.info(f"Created QC report: {qc_path}")
            
            # JSON export
            if export_json:
                json_path = Path(output_dir) / f"{output_prefix}_results.json"
                self.report_generator.export_to_json(batch_results, str(json_path))
                output_files.append(str(json_path))
                self.logger.info(f"Exported JSON: {json_path}")
            
            # Step 6: Finalize
            if self.progress_tracker:
                self.progress_tracker.update(6, "Analysis complete")
            
            processing_time = time.time() - start_time
            
            analysis_results.update({
                'success': True,
                'processing_time': processing_time,
                'output_files': output_files,
                'summary_stats': batch_results.summary_stats,
                'total_markers': batch_results.total_markers,
                'successful_matches': batch_results.successful_matches,
                'failed_analyses': batch_results.failed_analyses
            })
            
            self.logger.info(f"=== Analysis Complete ===")
            self.logger.info(f"Total processing time: {processing_time:.2f} seconds")
            self.logger.info(f"Successful matches: {batch_results.successful_matches}/{batch_results.total_markers}")
            self.logger.info(f"Output files: {len(output_files)}")
            
            return analysis_results
            
        except Exception as e:
            error_msg = f"Analysis failed: {e}"
            self.logger.error(error_msg)
            self.logger.error(traceback.format_exc())
            analysis_results['errors'].append(error_msg)
            return analysis_results
    
    def run_quick_analysis(self, vcf_filepath: str, excel_filepath: str,
                          output_filepath: str) -> Dict[str, Any]:
        """
        Run a quick analysis with minimal output.
        
        Args:
            vcf_filepath: Path to VCF file
            excel_filepath: Path to Excel file
            output_filepath: Output file path
            
        Returns:
            Dictionary with basic results
        """
        self.logger.info("Running quick analysis...")
        
        try:
            # Load markers
            markers = self.excel_processor.read_markers(excel_filepath)
            
            # Run analysis
            batch_results = self.marker_analyzer.analyze_markers_batch(
                markers, vcf_filepath, use_parallel=True
            )
            
            # Create summary report
            self.report_generator.create_summary_report(batch_results, output_filepath)
            
            return {
                'success': True,
                'total_markers': batch_results.total_markers,
                'successful_matches': batch_results.successful_matches,
                'output_file': output_filepath,
                'processing_time': batch_results.processing_time
            }
            
        except Exception as e:
            self.logger.error(f"Quick analysis failed: {e}")
            return {
                'success': False,
                'error': str(e)
            }
    
    def get_file_info(self, vcf_filepath: str, excel_filepath: str) -> Dict[str, Any]:
        """
        Get information about input files without running analysis.
        
        Args:
            vcf_filepath: Path to VCF file
            excel_filepath: Path to Excel file
            
        Returns:
            Dictionary with file information
        """
        info = {
            'vcf_info': {},
            'excel_info': {},
            'compatibility_check': {}
        }
        
        try:
            # VCF file info
            header = self.vcf_processor.get_vcf_header(vcf_filepath)
            variant_count = self.vcf_processor.get_variant_count(vcf_filepath)
            
            info['vcf_info'] = {
                'samples': header.samples,
                'sample_count': len(header.samples),
                'variant_count': variant_count,
                'file_size_mb': os.path.getsize(vcf_filepath) / (1024 * 1024)
            }
            
            # Excel file info
            excel_validation = self.excel_processor.validate_excel_file(excel_filepath)
            markers = self.excel_processor.read_markers(excel_filepath, validate=False)
            chromosomes = self.excel_processor.get_marker_chromosomes(markers)
            
            info['excel_info'] = {
                'total_rows': excel_validation.row_count,
                'valid_markers': excel_validation.marker_count,
                'chromosomes': sorted(list(chromosomes)),
                'chromosome_count': len(chromosomes),
                'detected_columns': excel_validation.detected_columns
            }
            
            # Compatibility check
            vcf_samples = set(header.samples)
            if len(vcf_samples) > 0:
                info['compatibility_check'] = {
                    'vcf_has_samples': True,
                    'markers_have_positions': excel_validation.marker_count > 0,
                    'estimated_analysis_time': self._estimate_analysis_time(
                        excel_validation.marker_count, len(header.samples)
                    )
                }
            
        except Exception as e:
            self.logger.error(f"Failed to get file info: {e}")
            info['error'] = str(e)
        
        return info
    
    def _estimate_analysis_time(self, marker_count: int, sample_count: int) -> str:
        """
        Estimate analysis time based on data size.
        
        Args:
            marker_count: Number of markers
            sample_count: Number of samples
            
        Returns:
            Estimated time string
        """
        # Rough estimation based on typical performance
        base_time = marker_count * 0.1  # 0.1 seconds per marker base
        sample_factor = max(1, sample_count / 100)  # Factor for large sample sets
        
        estimated_seconds = base_time * sample_factor
        
        if estimated_seconds < 60:
            return f"{estimated_seconds:.0f} seconds"
        elif estimated_seconds < 3600:
            return f"{estimated_seconds/60:.1f} minutes"
        else:
            return f"{estimated_seconds/3600:.1f} hours"


def create_argument_parser() -> argparse.ArgumentParser:
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description="Genetic Marker Analysis System",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  python main_controller.py data.vcf markers.xlsx output/
  
  # Quick analysis with summary only
  python main_controller.py data.vcf markers.xlsx output/ --quick
  
  # Custom configuration
  python main_controller.py data.vcf markers.xlsx output/ --config config.json
  
  # Disable parallel processing
  python main_controller.py data.vcf markers.xlsx output/ --no-parallel
        """
    )
    
    parser.add_argument('vcf_file', help='Path to VCF file (.vcf or .vcf.gz)')
    parser.add_argument('excel_file', help='Path to Excel file with marker data')
    parser.add_argument('output_dir', help='Output directory for results')
    
    parser.add_argument('--config', '-c', help='Configuration file path')
    parser.add_argument('--prefix', '-p', default='analysis', 
                       help='Prefix for output files (default: analysis)')
    
    parser.add_argument('--quick', '-q', action='store_true',
                       help='Run quick analysis (summary only)')
    parser.add_argument('--no-detailed', action='store_true',
                       help='Skip detailed Excel report')
    parser.add_argument('--no-summary', action='store_true',
                       help='Skip summary report')
    parser.add_argument('--no-qc', action='store_true',
                       help='Skip QC report')
    parser.add_argument('--json', '-j', action='store_true',
                       help='Export results to JSON')
    
    parser.add_argument('--no-parallel', action='store_true',
                       help='Disable parallel processing')
    parser.add_argument('--no-cache', action='store_true',
                       help='Disable VCF caching')
    parser.add_argument('--no-progress', action='store_true',
                       help='Hide progress messages')
    
    parser.add_argument('--info', '-i', action='store_true',
                       help='Show file information only (no analysis)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    
    return parser


def main():
    """Main entry point for command line interface."""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Initialize controller
    try:
        controller = GeneticMarkerAnalysisController(
            config_file=args.config,
            enable_caching=not args.no_cache,
            show_progress=not args.no_progress
        )
        
        # Set logging level if verbose
        if args.verbose:
            controller.config_manager.update_config(log_level='DEBUG')
        
        # File info mode
        if args.info:
            print("Getting file information...")
            info = controller.get_file_info(args.vcf_file, args.excel_file)
            
            print("\n=== File Information ===")
            print(f"VCF Samples: {info['vcf_info'].get('sample_count', 'unknown')}")
            print(f"VCF Variants: {info['vcf_info'].get('variant_count', 'unknown')}")
            print(f"VCF Size: {info['vcf_info'].get('file_size_mb', 0):.1f} MB")
            print(f"Excel Markers: {info['excel_info'].get('valid_markers', 'unknown')}")
            print(f"Excel Chromosomes: {info['excel_info'].get('chromosome_count', 'unknown')}")
            
            if 'estimated_analysis_time' in info.get('compatibility_check', {}):
                print(f"Estimated Analysis Time: {info['compatibility_check']['estimated_analysis_time']}")
            
            return
        
        # Quick analysis mode
        if args.quick:
            output_file = Path(args.output_dir) / f"{args.prefix}_summary.xlsx"
            result = controller.run_quick_analysis(
                args.vcf_file, args.excel_file, str(output_file)
            )
            
            if result['success']:
                print(f"\nQuick analysis complete!")
                print(f"Matches: {result['successful_matches']}/{result['total_markers']}")
                print(f"Output: {result['output_file']}")
                print(f"Time: {result['processing_time']:.2f} seconds")
            else:
                print(f"Quick analysis failed: {result['error']}")
                sys.exit(1)
            
            return
        
        # Full analysis
        result = controller.run_analysis(
            vcf_filepath=args.vcf_file,
            excel_filepath=args.excel_file,
            output_dir=args.output_dir,
            output_prefix=args.prefix,
            create_detailed_report=not args.no_detailed,
            create_summary_report=not args.no_summary,
            create_qc_report=not args.no_qc,
            export_json=args.json,
            use_parallel=not args.no_parallel
        )
        
        if result['success']:
            print(f"\n=== Analysis Complete ===")
            print(f"Successful matches: {result['successful_matches']}/{result['total_markers']}")
            print(f"Processing time: {result['processing_time']:.2f} seconds")
            print(f"Output files ({len(result['output_files'])}):")
            for file_path in result['output_files']:
                print(f"  - {file_path}")
            
            if result['warnings']:
                print(f"\nWarnings ({len(result['warnings'])}):")
                for warning in result['warnings'][:5]:  # Show first 5
                    print(f"  - {warning}")
                if len(result['warnings']) > 5:
                    print(f"  ... and {len(result['warnings']) - 5} more")
        else:
            print(f"\n=== Analysis Failed ===")
            print("Errors:")
            for error in result['errors']:
                print(f"  - {error}")
            sys.exit(1)
    
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"Fatal error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()