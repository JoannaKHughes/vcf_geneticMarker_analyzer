import argparse
import glob
import builtins
import pandas as pd
import argparse
import sys
import os
from pathlib import Path
import logging
import traceback
from datetime import datetime
import warnings
import time

# Command line example:
# python "C:\Users\joanna\Documents\Python_VSCode\Pangenome\src\100PioneerHorses_ISAG.py" --vcf_folder "C:\Users\joanna\Documents\Python_VSCode\Pangenome\vcf" --excel_file "C:\Users\joanna\Documents\Python_VSCode\Pangenome\ISAG_Marker_List.xlsx" --output_file "C:\Users\joanna\Documents\Python_VSCode\Pangenome\output\marker_comparison_multVCF.xlsx" --mapping_file "C:\Users\joanna\Documents\Python_VSCode\Pangenome\ISAG_2025_CT_Samples_from_Pioneer100.xlsx"
version = "3"

# Utility: auto-detect or prompt for column name
def resolve_column(df, logical_name, alternatives):
    """
    Try alternatives for a logical column name, prompt user if not found.
    Caches result in-memory for session.
    """
    if not hasattr(resolve_column, "_cache"):
        resolve_column._cache = {}
    cache = resolve_column._cache
    key = (id(df), logical_name)
    # Use cache if available
    if key in cache:
        return cache[key]
    # Try alternatives
    for alt in alternatives:
        if alt in df.columns:
            cache[key] = alt
            return alt
    # Prompt user
    print(f"\nAvailable columns: {list(df.columns)}")
    user_col = builtins.input(f"Column for '{logical_name}' not found. Please enter the correct column name: ")
    if user_col in df.columns:
        cache[key] = user_col
        return user_col
    raise KeyError(f"Column for '{logical_name}' not found and user input did not match.")

def main():
    ## Command line example:

    parser = argparse.ArgumentParser(description='Compare multiple VCF files in a folder to Excel marker file')
    parser.add_argument('--vcf_folder', required=True, help='Path to folder containing VCF files')
    parser.add_argument('--excel_file', required=True, help='Path to Excel marker file')
    parser.add_argument('--output_file', required=True, help='Path to output Excel file')
    parser.add_argument('--mapping_file', required=False, help='Path to sample mapping Excel file')
    args = parser.parse_args()

    vcf_folder = args.vcf_folder
    excel_file = args.excel_file
    output_file = args.output_file
    mapping_file = args.mapping_file

    # Discover all VCF files in the folder
    vcf_paths = sorted(glob.glob(os.path.join(vcf_folder, '*.vcf')))
    if not vcf_paths:
        logger.error(f"No VCF files found in folder: {vcf_folder}")
        sys.exit(1)

    # Performance timing
    script_start_time = time.time()
    
    # Pre-load all VCF genotype data for performance optimization
    vcf_genotype_cache = preload_vcf_genotype_data(vcf_paths)
    
    # Read all VCFs (for basic variant data only, genotypes are pre-cached)
    vcf_list = []
    vcf_sample_names_map = {}
    for vcf_path in vcf_paths:
        vcf_df, sample_names = read_vcf_file(vcf_path)
        vcf_list.append((os.path.basename(vcf_path), vcf_df))
        vcf_sample_names_map[os.path.basename(vcf_path)] = sample_names

    # Read Excel marker file
    excel_data = read_excel_file(excel_file)

    # Read sample mapping file if provided
    sample_mapping_df = None
    vcf_sample_mapping_df = None
    if mapping_file:
        sample_mapping_df = read_sample_mapping(mapping_file)
        # Use the first VCF's sample names for mapping
        first_vcf_name = os.path.basename(vcf_paths[0])
        sample_names = vcf_sample_names_map[first_vcf_name]
        sampleMap_output = build_vcf_sample_mapping(sample_mapping_df, sample_names)
        vcf_sample_mapping_df = pd.DataFrame(sampleMap_output)

    # Compare markers across all VCFs (with optimized genotype cache)
    comparison_start_time = time.time()
    results, detailed_stats, missing_allele_markers = compare_markers_multi(vcf_list, excel_data, vcf_genotype_cache)
    comparison_time = time.time() - comparison_start_time
    
    logger.info(f"Marker comparison completed in {comparison_time:.2f} seconds")

    # Collect marker positions for subset VCF (from first VCF)
    marker_positions = set()
    for sheet_name, df in results.items():
        for idx, row in df.iterrows():
            # Use the first VCF's match columns for marker extraction
            vcf1_match_col = [col for col in df.columns if col.endswith('_Match')][0]
            if row.get(vcf1_match_col, '') == "Match":
                chrom = row['Excel_CHROM']
                pos = row['Excel_POS']
                if chrom and not pd.isna(pos):
                    marker_positions.add((chrom, pos))

    # Save results
    save_results(results, detailed_stats, vcf_list[0][1], excel_data, output_file, marker_positions, vcf_sample_mapping_df, missing_allele_markers)

    # Final performance summary
    total_time = time.time() - script_start_time
    logger.info(f"Total script execution time: {total_time:.2f} seconds")
    logger.info("Comparison completed successfully with performance optimizations!")
 

def save_results(results, detailed_stats, vcf_df, excel_data, output_path, marker_positions, vcf_sample_mapping_df=None, missing_allele_markers=None):
    """
    Save comprehensive comparison results to Excel file with detailed statistics.
    Removes all columns ending with '_AlleleMatch' from results.
    Adds a sheet listing all markers missing REF or ALT in any VCF.
    """
    logger.info(f"Saving results to: {output_path}")

    # Defensive check for vcf_df type
    if not isinstance(vcf_df, pd.DataFrame):
        raise TypeError(f"vcf_df must be a pandas DataFrame, got {type(vcf_df)}. Check upstream code for errors in VCF loading.")

    try:
        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            # Save individual sheet results with all details
            for sheet_name, df in results.items():
                # Remove all columns ending with '_AlleleMatch'
                df = df.loc[:, ~df.columns.str.endswith('_AlleleMatch')]
                df.to_excel(writer, sheet_name=f"{sheet_name}_results", index=False)

            chrom_col = resolve_column(vcf_df, 'chromosome', ['CHROM_STD', 'CHROM', '#CHROM'])

            # Create new summary sheet with required columns
            # Gather all VCF names from the first results DataFrame
            vcf_match_cols = [col for col in results[list(results.keys())[0]].columns if col.endswith('_Match') and col != 'Match All']
            vcf_names = [col[:-6] for col in vcf_match_cols]  # Remove '_Match'

            summary_data = []
            for sheet_name, df in results.items():
                total = len(df)
                is_y = df['Excel_CHROM'] == 'Y'
                non_y = ~is_y
                non_y_total = non_y.sum()
                # For each VCF, count non-Y matches
                non_y_matches_per_vcf = []
                non_y_success_per_vcf = []
                for vcf_name in vcf_names:
                    matches = ((df[f"{vcf_name}_Match"] == "Match") & non_y).sum()
                    non_y_matches_per_vcf.append(matches)
                    rate = matches / non_y_total * 100 if non_y_total > 0 else 0
                    non_y_success_per_vcf.append(rate)
                # Overlapping: non-Y markers matched in any VCF
                overlap = ((df[[f"{vcf_name}_Match" for vcf_name in vcf_names]].eq("Match").any(axis=1)) & non_y).sum()
                overlap_rate = overlap / non_y_total * 100 if non_y_total > 0 else 0

                # Build row
                row = {
                    'ISAG_Sheet_Name': sheet_name,
                    'Total_Markers': total,
                    'Non_Y_Markers': non_y_total,
                }
                # Add per-VCF columns (assume 2 VCFs for this summary)
                for i, vcf_name in enumerate(vcf_names):
                    row[f'Non_Y_Matches to {vcf_name}'] = non_y_matches_per_vcf[i]
                for i, vcf_name in enumerate(vcf_names):
                    row[f'Non_Y_Success_Rate of {vcf_name}'] = f"{non_y_success_per_vcf[i]:.1f}%"
                row['Overlapping success rate'] = f"{overlap_rate:.1f}%"
                summary_data.append(row)

            summary_df = pd.DataFrame(summary_data)
            summary_df.to_excel(writer, sheet_name='Summary', index=False)

            # Create error report sheet
            all_errors = []
            for sheet_name, stats in detailed_stats.items():
                for error in stats.get('error_details', []):
                    all_errors.append({
                        'Sheet': sheet_name,
                        'Error': error,
                        'Timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                    })

            if all_errors:
                error_df = pd.DataFrame(all_errors)
                error_df.to_excel(writer, sheet_name='Error_Report', index=False)

            # Create unmatched markers sheet (use "Match All" column)
            unmatched_data = []
            for sheet_name, df in results.items():
                if "Match All" in df.columns:
                    unmatched_in_sheet = df[df['Match All'] != "Match"]
                else:
                    # fallback for legacy
                    unmatched_in_sheet = df[df['Match'] != "Match"]
                for idx, row in unmatched_in_sheet.iterrows():
                    chrom = row['Excel_CHROM']
                    pos = row['Excel_POS']

                    # Check if there's a VCF variant within 100bp (using first VCF)
                    nearby_variant = False
                    if chrom and not pd.isna(pos):
                        # Look for variants on the same chromosome within 100bp
                        chrom_col = resolve_column(vcf_df, 'chromosome', ['CHROM_STD', 'CHROM', '#CHROM'])
                        pos_col = resolve_column(vcf_df, 'position', ['POS', 'Position'])
                        nearby = vcf_df[
                            (vcf_df[chrom_col] == chrom) &
                            (vcf_df[pos_col] >= pos - 100) &
                            (vcf_df[pos_col] <= pos + 100)
                        ]
                        nearby_variant = not nearby.empty

                    # Find which VCFs this marker is matched in
                    matched_vcfs = []
                    for col in row.index:
                        if col.endswith('_Match') and col != 'Match All':
                            if row[col] == "Match":
                                matched_vcfs.append(col[:-6])
                    matched_vcfs_str = ", ".join(matched_vcfs) if matched_vcfs else "None"

                    unmatched_data.append({
                        'Sheet': sheet_name,
                        'Marker_ID': row['Marker_ID'],
                        'CHROM': chrom,
                        'POS': pos,
                        'Nearby_Variant': 'Yes' if nearby_variant else 'No',
                        'Matched_in_VCFs': matched_vcfs_str
                    })

            unmatched_df = pd.DataFrame(unmatched_data)
            if not unmatched_df.empty:
                unmatched_df.to_excel(writer, sheet_name='Unmatched_Markers', index=False)
            else:
                # Write empty sheet if no unmatched markers
                pd.DataFrame(columns=['Sheet','Marker_ID','CHROM','POS','Nearby_Variant','Matched_in_other_VCF?']).to_excel(
                    writer, sheet_name='Unmatched_Markers', index=False)

            # Create analysis metadata sheet with number of VCF files and overlapping success rate
            total_non_y_all_sheets = sum(stats['non_y_total'] for stats in detailed_stats.values())

            # Calculate overall overlapping success rate (non-Y markers matched in any VCF across all sheets)
            overlap_any_all_sheets = 0
            for sheet_name, df in results.items():
                non_y = df['Excel_CHROM'] != 'Y'
                vcf_match_cols = [col for col in df.columns if col.endswith('_Match') and col != 'Match All']
                if non_y.sum() > 0 and vcf_match_cols:
                    overlap_any_all_sheets += ((df[vcf_match_cols].eq("Match").any(axis=1)) & non_y).sum()
            overall_overlap_success_rate = overlap_any_all_sheets / total_non_y_all_sheets * 100 if total_non_y_all_sheets > 0 else 0

            metadata = {
                'Analysis_Date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'Script_Version': version,
                'Number_of_VCF_Files': len(results[list(results.keys())[0]].filter(like='_Match').columns) if results else 0,
                'VCF_File_Variants': len(vcf_df),
                'Excel_Sheets_Processed': len(excel_data),
                'Total_Markers_Analyzed': sum(len(df) for df in results.values()),
                'Non_Y_Markers_Analyzed': total_non_y_all_sheets,
                'Overall_Non_Y_Overlapping_Success_Rate': f"{overall_overlap_success_rate:.1f}%",
                'Author': 'Joanna Hughes'
            }

            metadata_data = [{'Parameter': k, 'Value': v} for k, v in metadata.items()]
            metadata_df = pd.DataFrame(metadata_data)
            metadata_df.to_excel(writer, sheet_name='Analysis_Metadata', index=False)


            # Add missing allele markers sheet
            if missing_allele_markers:
                missing_df = pd.DataFrame(missing_allele_markers)
                missing_df.to_excel(writer, sheet_name='Missing_Allele_Markers', index=False)
            else:
                pd.DataFrame(columns=['Sheet', 'Marker_ID', 'CHROM', 'POS', 'Missing_in_VCFs']).to_excel(
                    writer, sheet_name='Missing_Allele_Markers', index=False
                )

            # Add comprehensive README as a sheet
            readme_content = [
                ["README - 100PioneerHorses_ISAG.py Output Documentation"],
                [""],
                ["Overview:"],
                ["This Excel file is generated by the 100PioneerHorses_ISAG.py program, which compares genetic marker data from an Excel file (multiple sheets) with variant data from one or more VCF files for 100 pioneer horses."],
                [""],
                ["Workflow Summary:"],
                ["1. Reads all VCF files in a folder, each containing variant data for multiple horse samples."],
                ["2. Reads an Excel file with ISAG marker lists (multiple sheets)."],
                ["3. Reads a sample mapping file to relate VCF sample IDs to 'Finno Number' and ISAG CT 2025 IDs."],
                ["4. Compares markers from the Excel sheets to all VCF data, checking for exact matches by chromosome and position."],
                ["5. Outputs a comprehensive Excel file with detailed results, statistics, and metadata."],
                [""],
                ["Input Files:"],
                ["- VCF files: All .vcf files in the specified folder."],
                ["- Excel marker file: Contains marker information in multiple sheets (e.g., ISAG Core, ISAG Back up, X_Y)."],
                ["- Sample mapping file: Contains mapping between VCF sample IDs, Finno Number, and ISAG CT 2025 ID."],
                [""],
                ["Chromosome Standardization:"],
                ["- Chromosome names are standardized for comparison:"],
                ["  - 'eMSYv3' is mapped to 'Y'."],
                ["  - 'chr'/'Chr'/'CHR' prefixes are removed (e.g., 'chr1' becomes '1')."],
                ["  - Markers on the Y chromosome may not match due to the reference being female."],
                [""],
                ["Output Sheets:"],
                ["1. [SheetName]_results: Detailed comparison results for each marker in each input Excel sheet."],
                ["   - Columns include marker ID, chromosome, position, alleles, and for each VCF: <VCFNAME>_CHROM_exact, <VCFNAME>_POS_exact, <VCFNAME>_Match, plus 'Match All', plus allele frequency columns."],
                ["   - Columns 'freq ref of vcfX', 'freq alt of vcfX', 'freq ref of all vcf', 'freq alt of all vcf' are included for each marker."],
                ["2. Summary: For each input sheet, includes:"],
                ["   - ISAG_Sheet_Name, Total_Markers, Non_Y_Markers, Non_Y_Matches to VCF1, Non_Y_Matches to VCF2, Non_Y_Success_Rate of VCF1, Non_Y_Success_Rate of VCF2, Overlapping success rate (percent of non-Y markers matched in all VCFs)."],
                ["3. Unmatched_Markers: Markers from the Excel file that did not have an exact match in all VCFs."],
                ["   - For each unmatched marker, shows sheet name, marker ID, chromosome, position, whether a nearby variant (within 100bp) exists in the first VCF, and a column 'Matched_in_VCFs' listing the VCF(s) where the marker is matched (comma-separated, or 'None' if not matched in any)."],
                ["4. Error_Report: Any errors encountered during processing, such as invalid positions or missing data."],
                ["5. Analysis_Metadata: Metadata about the analysis, including date, script version, number of variants, and overall statistics."],
                ["6. Sample_Mapping: Lists all VCF sample names and, if available, the corresponding 'Finno Number' and 'ISAG CT 2025 ID' from the mapping file."],
                ["7. Missing_Allele_Markers: Markers missing REF or ALT allele in any VCF."],
                ["8. README: This documentation sheet."],
                [""],
                ["Details on Key Sheets:"],
                ["- [SheetName]_results: For each marker, shows whether an exact match was found in each VCF at the same chromosome and position. Includes reference and variant alleles, strand, and other marker details, plus allele frequencies."],
                ["- Summary: For each input sheet, shows total markers, non-Y markers, non-Y matches and success rates for each VCF, and the percent of non-Y markers matched in all VCFs."],
                ["- Unmatched_Markers: The 'Matched_in_VCFs' column lists which VCF(s) (by file name) the marker is matched in, or 'None' if not matched in any."],
                ["- Missing_Allele_Markers: Lists markers missing REF or ALT allele in any VCF."],
                ["- Sample_Mapping: Uses 'Finno Number' for matching VCF sample IDs. If a VCF sample matches a Finno Number in the mapping file, the corresponding ISAG CT 2025 ID is shown."],
                [""],
                ["Notes:"],
                ["- Unmatched markers are checked for nearby variants (within 100bp) in the first VCF to help identify potential close matches."],
                ["- All processing steps are logged, and errors are reported in the Error_Report sheet."],
                ["- The script is designed for flexibility and can be adapted for other marker lists or VCF files with similar structure."],
                [""],
                ["Contact:"],
                ["For questions or issues regarding this analysis or script, contact the program author or the VGL bioinformatics team."]
            ]
            readme_df = pd.DataFrame(readme_content)
            readme_df.to_excel(writer, sheet_name='README', index=False, header=False)

        logger.info("Results saved successfully with comprehensive statistics")

        # all VCF data are now included in the main Excel output as sheets.

        # Print final summary
        print("\n" + "="*60)
        print("ANALYSIS COMPLETE - FINAL SUMMARY")
        print("="*60)

        # Print new summary based on non-Y success rate
        for sheet_name, stats in detailed_stats.items():
            total = stats['total_markers']
            non_y_total = stats['non_y_total']
            non_y_matches = stats['non_y_matches']
            non_y_success_rate = non_y_matches / non_y_total * 100 if non_y_total > 0 else 0
            print(f"{sheet_name}: Non-Y markers: {non_y_matches}/{non_y_total} matched ({non_y_success_rate:.1f}%) [Total markers: {total}]")

        # Print overall non-Y success rate
        total_non_y_all_sheets = sum(stats['non_y_total'] for stats in detailed_stats.values())
        total_non_y_matches_all_sheets = sum(stats['non_y_matches'] for stats in detailed_stats.values())
        overall_non_y_success_rate = total_non_y_matches_all_sheets / total_non_y_all_sheets * 100 if total_non_y_all_sheets > 0 else 0
        print(f"Overall Non-Y Success Rate: {overall_non_y_success_rate:.1f}%")
        print(f"Results saved to: {output_path}")
        print("="*60)

    except Exception as e:
        logger.error(f"Error saving results: {e}")
        logger.error(traceback.format_exc())
        sys.exit(1)#!/usr/bin/env python3
"""
VCF and Excel Marker Comparison Script
Compares parentage markers from Excel file with VCF data from 100 pioneer horses
"""

# Suppress pandas warnings for cleaner output
warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def validate_file_exists(file_path, file_type):
    """
    Validate that file exists and is readable
    """
    if not Path(file_path).exists():
        logger.error(f"{file_type} file not found: {file_path}")
        return False
    
    if not os.access(file_path, os.R_OK):
        logger.error(f"{file_type} file is not readable: {file_path}")
        return False
    
    file_size = Path(file_path).stat().st_size
    if file_size == 0:
        logger.error(f"{file_type} file is empty: {file_path}")
        return False
    
    logger.info(f"{file_type} file validated: {file_path} (Size: {file_size} bytes)")
    return True

def preload_vcf_genotype_data(vcf_paths):
    """
    Pre-load all VCF genotype data and create position-based index for fast lookups.
    This eliminates the need to reopen VCF files for each marker during allele frequency calculation.
    
    Returns:
        dict: {vcf_filename: {(chrom, pos): genotype_data}}
    """
    logger.info("Pre-loading VCF genotype data for performance optimization...")
    start_time = time.time()
    
    vcf_genotype_cache = {}
    
    for vcf_path in vcf_paths:
        vcf_filename = os.path.basename(vcf_path)
        logger.info(f"Pre-loading genotype data from: {vcf_filename}")
        
        vcf_genotype_cache[vcf_filename] = {}
        
        try:
            with open(vcf_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 10:  # Need at least FORMAT and one sample
                        continue
                    
                    try:
                        chrom = fields[0]
                        pos = int(fields[1])
                        ref = fields[3]
                        alt = fields[4]
                        
                        # Standardize chromosome for indexing
                        std_chrom = standardize_chromosome(chrom)
                        
                        # Store genotype data for this position
                        genotypes = fields[9:]  # Sample genotype columns
                        
                        vcf_genotype_cache[vcf_filename][(std_chrom, pos)] = {
                            'ref': ref,
                            'alt': alt,
                            'genotypes': genotypes,
                            'line_num': line_num
                        }
                        
                    except (ValueError, IndexError) as e:
                        # Skip malformed lines
                        continue
                        
        except Exception as e:
            logger.error(f"Error pre-loading genotype data from {vcf_filename}: {e}")
            continue
    
    load_time = time.time() - start_time
    total_positions = sum(len(cache) for cache in vcf_genotype_cache.values())
    logger.info(f"Pre-loaded genotype data for {total_positions} positions from {len(vcf_paths)} VCF files in {load_time:.2f} seconds")
    
    return vcf_genotype_cache

def read_vcf_file(vcf_path):
    """
    Read VCF file and extract relevant information with validation
    Returns a DataFrame with CHROM, POS, REF, ALT columns and list of sample names
    """
    logger.info(f"Reading VCF file: {vcf_path}")
    
    if not validate_file_exists(vcf_path, "VCF"):
        sys.exit(1)
    
    vcf_data = []
    line_count = 0
    header_lines = 0
    error_lines = []
    sample_names = []
    
    try:
        # Log file path for debugging
        logger.info(f"Reading VCF file: {vcf_path} (Size: {os.path.getsize(vcf_path)} bytes)")
        
        # DIAGNOSTIC: Test file encoding and readability
        logger.info("DIAGNOSTIC: Testing file encoding and initial content...")
        
        # Test 1: Check if file can be opened with different encodings
        test_encodings = ['utf-8', 'utf-8-sig', 'latin-1', 'cp1252']
        working_encoding = None
        
        for enc in test_encodings:
            try:
                with open(vcf_path, 'r', encoding=enc) as test_f:
                    # Try to read first 1000 characters
                    test_content = test_f.read(1000)
                    logger.info(f"DIAGNOSTIC: Encoding '{enc}' - SUCCESS. First 100 chars: {repr(test_content[:100])}")
                    working_encoding = enc
                    break
            except UnicodeDecodeError as e:
                logger.info(f"DIAGNOSTIC: Encoding '{enc}' - FAILED: {e}")
            except Exception as e:
                logger.error(f"DIAGNOSTIC: Encoding '{enc}' - ERROR: {e}")
        
        if not working_encoding:
            logger.error("DIAGNOSTIC: No compatible encoding found!")
            # Try binary mode to check for binary content
            with open(vcf_path, 'rb') as bin_f:
                first_bytes = bin_f.read(100)
                logger.error(f"DIAGNOSTIC: First 100 bytes (binary): {first_bytes}")
            raise ValueError("Unable to determine file encoding")
        
        # Test 2: Check line endings and file structure
        logger.info(f"DIAGNOSTIC: Using encoding '{working_encoding}' for detailed analysis...")
        with open(vcf_path, 'rb') as bin_f:
            sample_bytes = bin_f.read(1000)
            crlf_count = sample_bytes.count(b'\r\n')
            lf_count = sample_bytes.count(b'\n') - crlf_count
            cr_count = sample_bytes.count(b'\r') - crlf_count
            logger.info(f"DIAGNOSTIC: Line endings - CRLF: {crlf_count}, LF: {lf_count}, CR: {cr_count}")
            
            # Check for null bytes or other binary markers
            null_bytes = sample_bytes.count(b'\x00')
            if null_bytes > 0:
                logger.warning(f"DIAGNOSTIC: Found {null_bytes} null bytes - possible binary content!")
        
        # Use the working encoding for actual file processing
        logger.info(f"DIAGNOSTIC: Proceeding with encoding '{working_encoding}'")
        with open(vcf_path, 'r', encoding=working_encoding, errors='replace') as f:
            for line_num, line in enumerate(f, 1):
                # DIAGNOSTIC: Log problematic lines
                if line_num <= 5:
                    logger.info(f"DIAGNOSTIC: Line {line_num} length={len(line)}, ends_with={repr(line[-10:]) if len(line) >= 10 else repr(line)}")
                
                # Check for problematic characters every 1000 lines
                if line_num % 1000 == 0:
                    try:
                        line.encode('utf-8')
                    except UnicodeEncodeError as e:
                        logger.warning(f"DIAGNOSTIC: Line {line_num} has encoding issues: {e}")
                line_count += 1
                
                # Debug: log first 5 lines
                if line_num <= 5:
                    logger.debug(f"Line {line_num}: {line.strip()}")
                
                if line.startswith('#'): 
                    # Capture sample names from header
                    if line.startswith('#CHROM'):
                        fields = line.strip().split('\t')
                        if len(fields) >= 10:
                            sample_names = fields[9:]
                            logger.info(f"Found {len(sample_names)} samples in VCF header")
                    header_lines += 1
                    continue
                
                # Parse VCF line
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    error_lines.append(f"Line {line_num}: Insufficient fields ({len(fields)} < 5)")
                    continue
                
                try:
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    
                    # Basic validation
                    if not chrom:
                        error_lines.append(f"Line {line_num}: Empty chromosome")
                        continue
                    
                    if pos <= 0:
                        error_lines.append(f"Line {line_num}: Invalid position ({pos})")
                        continue
                    
                    if not ref or not alt:
                        error_lines.append(f"Line {line_num}: Empty alleles")
                        continue
                    
                    vcf_data.append({
                        'CHROM': chrom,
                        'POS': pos,
                        'REF': ref,
                        'ALT': alt,
                        'LINE_NUM': line_num
                    })
                    
                except ValueError as e:
                    error_lines.append(f"Line {line_num}: {e}")
                    continue
    
    except FileNotFoundError:
        logger.error(f"VCF file not found: {vcf_path}")
        sys.exit(1)
    except OSError as e:
        logger.error(f"OS error reading VCF file: {vcf_path}")
        logger.error(f"Error details: {e.errno} - {e.strerror}")
        logger.error(f"File info: Size={os.path.getsize(vcf_path)} bytes, Exists={os.path.exists(vcf_path)}")
        
        # DIAGNOSTIC: Enhanced handling for Errno 22 (Invalid argument)
        if e.errno == 22:  # ERROR_INVALID_PARAMETER on Windows
            logger.error("DIAGNOSTIC: [Errno 22] Invalid argument detected!")
            logger.error("DIAGNOSTIC: This typically indicates:")
            logger.error("DIAGNOSTIC: 1. File encoding issues (non-UTF-8 characters)")
            logger.error("DIAGNOSTIC: 2. File corruption or binary data in text file")
            logger.error("DIAGNOSTIC: 3. Invalid file path characters")
            logger.error("DIAGNOSTIC: 4. File being accessed by another process")
            
            # Additional diagnostics for Errno 22
            try:
                import stat
                file_stat = os.stat(vcf_path)
                logger.error(f"DIAGNOSTIC: File mode: {stat.filemode(file_stat.st_mode)}")
                logger.error(f"DIAGNOSTIC: File size: {file_stat.st_size} bytes")
                logger.error(f"DIAGNOSTIC: Last modified: {file_stat.st_mtime}")
                
                # Check if file is actually readable
                with open(vcf_path, 'rb') as test_file:
                    first_kb = test_file.read(1024)
                    logger.error(f"DIAGNOSTIC: First 1KB binary read successful")
                    logger.error(f"DIAGNOSTIC: Contains null bytes: {b'\\x00' in first_kb}")
                    logger.error(f"DIAGNOSTIC: First 50 bytes: {first_kb[:50]}")
                    
            except Exception as diag_e:
                logger.error(f"DIAGNOSTIC: Additional diagnosis failed: {diag_e}")
        
        sys.exit(1)
    except UnicodeDecodeError as e:
        logger.error(f"DIAGNOSTIC: Unicode decode error in VCF file: {vcf_path}")
        logger.error(f"DIAGNOSTIC: Error details: {e}")
        logger.error(f"DIAGNOSTIC: Error position: {e.start}-{e.end}")
        logger.error(f"DIAGNOSTIC: Problematic bytes: {e.object[e.start:e.end] if hasattr(e, 'object') else 'N/A'}")
        
        # Try to provide file encoding suggestions
        logger.error("DIAGNOSTIC: Suggested solutions:")
        logger.error("DIAGNOSTIC: 1. Check file encoding - may not be UTF-8")
        logger.error("DIAGNOSTIC: 2. File may contain binary data")
        logger.error("DIAGNOSTIC: 3. Try opening with different encoding (latin-1, cp1252)")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error reading VCF file: {vcf_path}")
        logger.error(f"Error: {type(e).__name__} - {str(e)}")
        logger.error(traceback.format_exc())
        sys.exit(1)
    
    vcf_df = pd.DataFrame(vcf_data)
    
    # Report parsing statistics
    logger.info(f"VCF parsing statistics:")
    logger.info(f"  Total lines: {line_count}")
    logger.info(f"  Header lines: {header_lines}")
    logger.info(f"  Data lines parsed: {len(vcf_df)}")
    logger.info(f"  Error lines: {len(error_lines)}")
    
    if error_lines:
        logger.warning(f"First 5 parsing errors: {error_lines[:5]}")
    
    if len(vcf_df) == 0:
        logger.error("No valid variants found in VCF file")
        sys.exit(1)
    
    # Validate data integrity
    logger.info("Validating VCF data integrity...")
    unique_positions = vcf_df.groupby(['CHROM', 'POS']).size()
    duplicates = unique_positions[unique_positions > 1]
    if len(duplicates) > 0:
        logger.warning(f"Found {len(duplicates)} duplicate positions")
    
    logger.info(f"Loaded {len(vcf_df)} variants from VCF file")
    
    return vcf_df, sample_names

def read_excel_file(excel_path):
    """
    Read Excel file with multiple sheets (core, backup, x_y) with validation
    Returns a dictionary of DataFrames
    """
    logger.info(f"Reading Excel file: {excel_path}")
    
    if not validate_file_exists(excel_path, "Excel"):
        sys.exit(1)
    
    try:
        # First, check what sheets are available
        xl_file = pd.ExcelFile(excel_path)
        available_sheets = xl_file.sheet_names
        logger.info(f"Available sheets: {available_sheets}")
        
        expected_sheets = ['ISAG Core', 'ISAG Back up', 'X_Y']
        missing_sheets = [sheet for sheet in expected_sheets if sheet not in available_sheets]
        
        if missing_sheets:
            logger.warning(f"Expected sheets not found: {missing_sheets}")
            # Try to find sheets with similar names
            for missing in missing_sheets:
                similar = [s for s in available_sheets if missing.lower() in s.lower() or s.lower() in missing.lower()]
                if similar:
                    logger.info(f"Similar sheets found for '{missing}': {similar}")
        
        # Read available expected sheets
        sheets_to_read = [sheet for sheet in expected_sheets if sheet in available_sheets]
        if not sheets_to_read:
            logger.error("No expected sheets found in Excel file")
            sys.exit(1)
        
        excel_data = pd.read_excel(excel_path, sheet_name=sheets_to_read)
        
        # If only one sheet was read, convert to dict format
        if not isinstance(excel_data, dict):
            excel_data = {sheets_to_read[0]: excel_data}
        
        # Validate each sheet
        required_columns = ['CHROM', 'POS', 'IDs', 'List', 'Reference Allele', 'Variant Allele', 'Strand', 'Sequence']
        
        for sheet_name, df in excel_data.items():
            logger.info(f"Validating sheet '{sheet_name}': {len(df)} rows, {len(df.columns)} columns")
            
            # Check for empty dataframe
            if df.empty:
                logger.warning(f"Sheet '{sheet_name}' is empty")
                continue
            
            # Check columns
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                logger.warning(f"Sheet '{sheet_name}' missing columns: {missing_columns}")
                logger.info(f"Available columns: {list(df.columns)}")
            
            # Check for required data
            essential_columns = ['CHROM', 'POS']
            for col in essential_columns:
                if col in df.columns:
                    null_count = df[col].isnull().sum()
                    if null_count > 0:
                        logger.warning(f"Sheet '{sheet_name}' has {null_count} null values in column '{col}'")
            
            # Data type validation
            if 'POS' in df.columns:
                try:
                    df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
                    invalid_positions = df['POS'].isnull().sum()
                    if invalid_positions > 0:
                        logger.warning(f"Sheet '{sheet_name}' has {invalid_positions} invalid positions")
                except Exception as e:
                    logger.warning(f"Could not convert POS to numeric in sheet '{sheet_name}': {e}")
        
        return excel_data
    
    except FileNotFoundError:
        logger.error(f"Excel file not found: {excel_path}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error reading Excel file: {e}")
        logger.error(traceback.format_exc())
        sys.exit(1)

def calculate_allele_frequencies_from_cache(vcf_genotype_cache, vcf_filename, chrom, pos):
    """
    Calculate allele frequencies from pre-cached genotype data.
    This replaces the file reopening bottleneck with fast dictionary lookup.
    
    Returns:
        tuple: (freq_ref, freq_alt, total_alleles, missing_allele)
    """
    genotype_data = vcf_genotype_cache.get(vcf_filename, {}).get((chrom, pos))
    
    if not genotype_data:
        return 0.0, 0.0, 0, True
    
    vcf_ref = genotype_data['ref']
    vcf_alt = genotype_data['alt']
    genotypes = genotype_data['genotypes']
    
    # Check if both REF and ALT are present
    if not vcf_ref or not vcf_alt:
        return 0.0, 0.0, 0, True
    
    # Calculate allele frequencies from genotypes
    ref_count = 0
    alt_count = 0
    non_missing = 0
    
    for gt in genotypes:
        gt_field = gt.split(':')[0].replace('|', '/')
        if gt_field in ['0/0']:
            ref_count += 2
            non_missing += 1
        elif gt_field in ['1/1']:
            alt_count += 2
            non_missing += 1
        elif gt_field in ['0/1', '1/0']:
            ref_count += 1
            alt_count += 1
            non_missing += 1
        # skip missing (e.g., './.')
    
    total_alleles = 2 * non_missing
    freq_ref = ref_count / total_alleles if total_alleles > 0 else 0.0
    freq_alt = alt_count / total_alleles if total_alleles > 0 else 0.0
    
    return freq_ref, freq_alt, total_alleles, False

def standardize_chromosome(chrom):
    """
    Standardize chromosome naming:
    - Map 'eMSYv3' to 'Y'
    - Remove 'chr' prefix if present for other chromosomes
    """
    if isinstance(chrom, str):
        if chrom == 'eMSYv3':
            return 'Y'
        return chrom.replace('chr', '').replace('Chr', '').replace('CHR', '')
    return str(chrom)

def compare_markers_multi(vcf_list, excel_data, vcf_genotype_cache):
    """
    Compare markers between multiple VCFs and Excel data, including allele matching and allele frequency calculation.
    Returns comparison results for each sheet and a list of markers missing REF or ALT in any VCF.
    vcf_list: list of (vcf_filename, vcf_df)
    excel_data: dict of DataFrames (sheet_name -> df)
    vcf_genotype_cache: dict with pre-loaded genotype data {vcf_filename: {(chrom, pos): genotype_data}}
    """
    logger.info("Starting multi-VCF marker comparison...")

    # Standardize chromosome names in all VCFs
    vcf_std_list = []
    for vcf_filename, vcf_df in vcf_list:
        vcf_df = vcf_df.copy()
        vcf_df['CHROM_STD'] = vcf_df['CHROM'].apply(standardize_chromosome)
        vcf_std_list.append((vcf_filename, vcf_df))

    results = {}
    detailed_stats = {}

    missing_allele_markers = []  # Track markers missing REF or ALT in any VCF

    for sheet_name, excel_df in excel_data.items():
        logger.info(f"Comparing sheet: {sheet_name}")

        if excel_df.empty:
            logger.warning(f"Sheet '{sheet_name}' is empty, skipping")
            continue

        # Standardize chromosome names in Excel
        excel_df = excel_df.copy()
        excel_df['CHROM_STD'] = excel_df['CHROM'].apply(standardize_chromosome)

        # Initialize detailed statistics
        stats = {
            'total_markers': len(excel_df),
            'valid_markers': 0,
            'invalid_markers': 0,
            'position_matches': 0,
            'complete_matches': 0,
            'missing_markers': 0,
            'allele_conflicts': 0,
            'processing_errors': 0,
            'chromosome_distribution': {},
            'position_range': {'min': None, 'max': None},
            'error_details': [],
            'y_markers': 0,
            'non_y_matches': 0,
            'non_y_total': 0
        }

        # Create comparison results
        comparison_results = []

        for idx, row in excel_df.iterrows():
            try:
                marker_id = row.get('IDs', f"Unknown_{idx}")
                chrom = row.get('CHROM_STD', '')
                pos = row.get('POS', None)
                ref_allele = str(row.get('Reference Allele', '')).strip()
                var_allele = str(row.get('Variant Allele', '')).strip()
                strand = row.get('Strand', '')
                sequence = row.get('Sequence', '')
                list_info = row.get('List', '')

                # For frequency calculation
                per_vcf_freq = []
                missing_in_vcf = []
                total_ref_all_vcf = 0
                total_alt_all_vcf = 0
                total_alleles_all_vcf = 0

                # Validate essential fields
                if pd.isna(pos) or pos == '' or pos <= 0:
                    stats['invalid_markers'] += 1
                    stats['error_details'].append(f"Marker {marker_id}: Invalid position")
                    result_row = {
                        'Marker_ID': marker_id,
                        'Excel_CHROM': chrom,
                        'Excel_POS': pos,
                        'Excel_Ref': ref_allele,
                        'Excel_Var': var_allele,
                        'Strand': strand,
                        'List': list_info,
                    }
                    for vcf_filename, _ in vcf_std_list:
                        result_row[f"{vcf_filename}_CHROM_exact"] = ''
                        result_row[f"{vcf_filename}_POS_exact"] = ''
                        result_row[f"{vcf_filename}_Match"] = 'Error'
                        result_row[f"{vcf_filename}_AlleleMatch"] = ''
                    result_row['Match All'] = 'Error'
                    comparison_results.append(result_row)
                    continue

                if not chrom:
                    stats['invalid_markers'] += 1
                    stats['error_details'].append(f"Marker {marker_id}: Missing chromosome")
                    result_row = {
                        'Marker_ID': marker_id,
                        'Excel_CHROM': chrom,
                        'Excel_POS': pos,
                        'Excel_Ref': ref_allele,
                        'Excel_Var': var_allele,
                        'Strand': strand,
                        'List': list_info,
                    }
                    for vcf_filename, _ in vcf_std_list:
                        result_row[f"{vcf_filename}_CHROM_exact"] = ''
                        result_row[f"{vcf_filename}_POS_exact"] = ''
                        result_row[f"{vcf_filename}_Match"] = 'Error'
                        result_row[f"{vcf_filename}_AlleleMatch"] = ''
                    result_row['Match All'] = 'Error'
                    comparison_results.append(result_row)
                    continue

                stats['valid_markers'] += 1

                # Track Y markers separately
                is_y_marker = (chrom == 'Y')
                if is_y_marker:
                    stats['y_markers'] += 1
                else:
                    stats['non_y_total'] += 1

                # Update chromosome distribution
                if chrom not in stats['chromosome_distribution']:
                    stats['chromosome_distribution'][chrom] = 0
                stats['chromosome_distribution'][chrom] += 1

                # Update position range
                if stats['position_range']['min'] is None or pos < stats['position_range']['min']:
                    stats['position_range']['min'] = pos
                if stats['position_range']['max'] is None or pos > stats['position_range']['max']:
                    stats['position_range']['max'] = pos

                # For each VCF, check for exact match and allele match
                vcf_match_cols = []
                vcf_allele_cols = []
                all_matched = True
                all_allele_matched = True
                for vcf_filename, vcf_df in vcf_std_list:
                    vcf_exact_match = vcf_df[
                        (vcf_df['CHROM_STD'] == chrom) &
                        (vcf_df['POS'] == pos)
                    ]
                    # Frequency calculation for this VCF
                    freq_ref = 0.0
                    freq_alt = 0.0
                    total_alleles = 0
                    missing_allele = False

                    if len(vcf_exact_match) > 0:
                        exact_match = vcf_exact_match.iloc[0]
                        vcf_chr = exact_match['CHROM']
                        vcf_pos = exact_match['POS']
                        vcf_ref = str(exact_match['REF']).strip()
                        vcf_alt = str(exact_match['ALT']).strip()
                        # Compare alleles
                        if (vcf_ref == ref_allele) and (vcf_alt == var_allele):
                            allele_status = "AlleleMatch"
                        else:
                            allele_status = "AlleleMismatch"
                            all_allele_matched = False
                        match_status = "Match" if allele_status == "AlleleMatch" else "No Match"

                        # Count REF/ALT alleles in all samples for this marker in this VCF
                        # VCF sample columns start at index 9
                        # Only count if both REF and ALT are present
                        # OPTIMIZED: Use pre-cached genotype data instead of reopening VCF files
                        freq_ref, freq_alt, total_alleles, missing_allele = calculate_allele_frequencies_from_cache(
                            vcf_genotype_cache, vcf_filename, chrom, pos
                        )
                        
                        # Accumulate totals for all VCFs frequency calculation
                        if not missing_allele and total_alleles > 0:
                            ref_count = int(freq_ref * total_alleles)
                            alt_count = int(freq_alt * total_alleles)
                            total_ref_all_vcf += ref_count
                            total_alt_all_vcf += alt_count
                            total_alleles_all_vcf += total_alleles
                    else:
                        vcf_chr = ''
                        vcf_pos = ''
                        vcf_ref = ''
                        vcf_alt = ''
                        allele_status = ''
                        match_status = "No Match"
                        all_matched = False
                        all_allele_matched = False
                        missing_allele = True
                        freq_ref = 0.0
                        freq_alt = 0.0
                        total_alleles = 0

                    # Track missing alleles for this marker/VCF
                    if missing_allele:
                        missing_in_vcf.append(vcf_filename)

                    #shorten vcf filename to 15 chars
                    short_col_name = vcf_filename[:15]
                    vcf_match_cols.append((short_col_name, vcf_chr, vcf_pos, match_status))
                    vcf_allele_cols.append((short_col_name, allele_status))
                    per_vcf_freq.append((short_col_name, freq_ref, freq_alt))

                # For non-Y markers, track matches (all VCFs must match)
                if not is_y_marker and all([m[3] == "Match" for m in vcf_match_cols]):
                    stats['non_y_matches'] += 1

                # Build result row
                result_row = {
                    'Marker_ID': marker_id,
                    'Excel_CHROM': chrom,
                    'Excel_POS': pos,
                    'Excel_Ref': ref_allele,
                    'Excel_Var': var_allele,
                    'Strand': strand,
                    'List': list_info,
                }
                for (vcf_filename, vcf_chr, vcf_pos, match_status), (_, allele_status) in zip(vcf_match_cols, vcf_allele_cols):
                    result_row[f"{vcf_filename}_CHROM_exact"] = vcf_chr
                    result_row[f"{vcf_filename}_POS_exact"] = vcf_pos
                    result_row[f"{vcf_filename}_Match"] = match_status
                    result_row[f"{vcf_filename}_AlleleMatch"] = allele_status
                # "Match All" only if all VCFs match by position and alleles
                result_row['Match All'] = "Match" if all([m[3] == "Match" for m in vcf_match_cols]) else "No Match"

                # Add frequency columns for each VCF
                for (vcf_short, freq_ref, freq_alt) in per_vcf_freq:
                    result_row[f'freq ref of {vcf_short}'] = freq_ref
                    result_row[f'freq alt of {vcf_short}'] = freq_alt

                # Add overall frequency columns
                freq_ref_all = total_ref_all_vcf / total_alleles_all_vcf if total_alleles_all_vcf > 0 else 0.0
                freq_alt_all = total_alt_all_vcf / total_alleles_all_vcf if total_alleles_all_vcf > 0 else 0.0
                result_row['freq ref of all vcf'] = freq_ref_all
                result_row['freq alt of all vcf'] = freq_alt_all

                # Track markers missing REF or ALT in any VCF
                if missing_in_vcf:
                    missing_allele_markers.append({
                        'Sheet': sheet_name,
                        'Marker_ID': marker_id,
                        'CHROM': chrom,
                        'POS': pos,
                        'Missing_in_VCFs': ', '.join(missing_in_vcf)
                    })

                comparison_results.append(result_row)

            except Exception as e:
                stats['processing_errors'] += 1
                stats['error_details'].append(f"Marker {idx}: Processing error - {str(e)}")
                logger.error(f"Error processing marker {idx} in sheet {sheet_name}: {e}")

                result_row = {
                    'Marker_ID': f"Error_{idx}",
                    'Excel_CHROM': '',
                    'Excel_POS': '',
                    'Excel_Ref': '',
                    'Excel_Var': '',
                    'Strand': '',
                    'List': '',
                }
                for vcf_filename, _ in vcf_std_list:
                    result_row[f"{vcf_filename}_CHROM_exact"] = ''
                    result_row[f"{vcf_filename}_POS_exact"] = ''
                    result_row[f"{vcf_filename}_Match"] = 'Error'
                    result_row[f"{vcf_filename}_AlleleMatch"] = ''
                result_row['Match All'] = 'Error'
                comparison_results.append(result_row)

        results[sheet_name] = pd.DataFrame(comparison_results)
        detailed_stats[sheet_name] = stats

        # Log detailed statistics
        logger.info(f"Sheet '{sheet_name}' detailed statistics:")
        logger.info(f"  Total markers: {stats['total_markers']}")
        logger.info(f"  Valid markers: {stats['valid_markers']}")
        logger.info(f"  Invalid markers: {stats['invalid_markers']}")
        logger.info(f"  Y chromosome markers: {stats['y_markers']} (excluded from success rate)")
        if stats['non_y_total'] > 0:
            non_y_success_rate = stats['non_y_matches'] / stats['non_y_total'] * 100
            logger.info(f"  Non-Y markers: {stats['non_y_total']}, Matches: {stats['non_y_matches']} ({non_y_success_rate:.1f}%)")
        else:
            logger.info("  No non-Y markers to calculate success rate")
        logger.info(f"  Processing errors: {stats['processing_errors']}")

        if stats['chromosome_distribution']:
            logger.info(f"  Chromosome distribution: {dict(sorted(stats['chromosome_distribution'].items()))}")

    return results, detailed_stats, missing_allele_markers

# Removed duplicate/old save_results function to avoid confusion.

   

def create_isag_marker_subset_vcf(excel_data, vcf_path, output_path):
    """
    Create VCF subset for markers found in ISAG Excel sheets (core, backup, X_Y)
    excluding Y chromosome markers
    """
    logger.info("Extracting ISAG marker positions (excluding Y chromosome)")
    marker_positions = set()
    
    for sheet_name, df in excel_data.items():
        # Standardize chromosome names and filter out Y
        df['CHROM_STD'] = df['CHROM'].apply(standardize_chromosome)
        non_y_df = df[df['CHROM_STD'] != 'Y']
        
        for idx, row in non_y_df.iterrows():
            chrom = row['CHROM_STD']
            pos = row['POS']
            if not pd.isna(pos) and pos != '' and pos > 0:
                marker_positions.add((chrom, pos))
    
    logger.info(f"Found {len(marker_positions)} non-Y marker positions across all sheets")
    
    # Create the subset VCF
    logger.info(f"Creating subset VCF at: {output_path}")
    try:
        with open(vcf_path, 'r') as f_in, open(output_path, 'w') as f_out:
            # Write header lines
            for line in f_in:
                if line.startswith('#'):
                    f_out.write(line)
                else:
                    break
            
            # Write matching variants
            f_in.seek(0)
            for line in f_in:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 2:
                    continue
                    
                chrom = fields[0]
                pos = int(fields[1])
                std_chrom = standardize_chromosome(chrom)
                
                if (std_chrom, pos) in marker_positions:
                    f_out.write(line)
                    
        logger.info(f"Subset VCF created successfully with {len(marker_positions)} markers")
        
    except Exception as e:
        logger.error(f"Error creating subset VCF: {e}")
        logger.error(traceback.format_exc())
        sys.exit(1)

def create_marker_subset_vcf(vcf_path, marker_positions, output_path):
    """
    Create subset VCF containing only specified marker positions
    """
    logger.info(f"Creating subset VCF for {len(marker_positions)} markers")
    
    try:
        with open(vcf_path, 'r') as f_in, open(output_path, 'w') as f_out:
            # Write header lines
            for line in f_in:
                if line.startswith('#'):
                    f_out.write(line)
                else:
                    break
            
            # Write matching variants
            f_in.seek(0)
            for line in f_in:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 2:
                    continue
                    
                chrom = fields[0]
                pos = int(fields[1])
                std_chrom = standardize_chromosome(chrom)
                
                if (std_chrom, pos) in marker_positions:
                    f_out.write(line)
                    
    except Exception as e:
        logger.error(f"Error creating subset VCF: {e}")
        logger.error(traceback.format_exc())
        sys.exit(1)
        
    # Build VCF sample mapping DataFrame for output sheet using a hashmap for fast lookup
def build_vcf_sample_mapping(sample_mapping_df, sample_names):
    # Step 1: Build the hashmap from the sample mapping file
    mapping_dict = {}
    if not sample_mapping_df.empty:
        for _, row in sample_mapping_df.iterrows():
            key = str(row['VCF_Sample_ID']).strip()
            mapping_dict[key] = {
                'Finno Number': row['Finno Number'],
                'ISAG_CT_2025_ID': row['ISAG_CT_2025_ID']
            }
    # Step 2: Get VCF sample names
    vcf_sample_names = sample_names
    # Step 3: Build output data
    output_rows = []
    for sample in vcf_sample_names:
        lookup_key = str(sample).strip()
        info = mapping_dict.get(lookup_key, {'Finno Number': None, 'ISAG_CT_2025_ID': None})
        output_rows.append({
            'VCF_Sample_ID': sample,
            'Finno Number': info['Finno Number'],
            'ISAG_CT_2025_ID': info['ISAG_CT_2025_ID']
        })
    return output_rows
            


if __name__ == "__main__":
    main()
