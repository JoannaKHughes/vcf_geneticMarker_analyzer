# 100PioneerHorses_ISAG.py
# Joanna Hughes June 26 2025
# This script compares genetic marker data from Excel sheets with variant data from multiple VCF files.
# It outputs a comprehensive Excel file with detailed results, statistics, and metadata for 100 pioneer horses.

import argparse   # For parsing command-line arguments
import glob       # For file pattern matching (finding VCF files)
import builtins   # For using built-in input() in custom prompts
import pandas as pd  # For data manipulation and Excel I/O
import argparse
import sys
import os
from pathlib import Path  # For robust file path handling
import logging   # For logging progress and errors
import traceback # For detailed error reporting
from datetime import datetime  # For timestamps in logs and output
import warnings   # For suppressing pandas warnings

# Set up logging at the very top before any logger usage
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Command line example:
# python "C:\Users\joanna\Documents\Python_VSCode\Pangenome\src\100PioneerHorses_ISAG.py" --vcf_folder "C:\Users\joanna\Documents\Python_VSCode\Pangenome\vcf" --excel_file "C:\Users\joanna\Documents\Python_VSCode\Pangenome\ISAG_Marker_List.xlsx" --output_file "C:\Users\joanna\Documents\Python_VSCode\Pangenome\output\marker_comparison_multVCF.xlsx" --mapping_file "C:\Users\joanna\Documents\Python_VSCode\Pangenome\ISAG_2025_CT_Samples_from_Pioneer100.xlsx"

version = "3"  # Script version for metadata

# Utility: auto-detect or prompt for column name
def resolve_column(df, logical_name, alternatives):
    """
    Utility function to resolve a logical column name in a DataFrame.
    - Tries a list of alternative names.
    - If not found, prompts the user for the correct column name.
    - Caches the result for the session to avoid repeated prompts.
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
    # Prompt user if no match found
    print(f"\nAvailable columns: {list(df.columns)}")
    user_col = builtins.input(f"Column for '{logical_name}' not found. Please enter the correct column name: ")
    if user_col in df.columns:
        cache[key] = user_col
        return user_col
    raise KeyError(f"Column for '{logical_name}' not found and user input did not match.")

def main():
    """
    Main entry point for the script.
    - Parses command-line arguments.
    - Loads VCF and Excel files.
    - Compares markers across all VCFs.
    - Saves results to output Excel file.
    """
    # Parse command-line arguments for input/output files and folders
    parser = argparse.ArgumentParser(description='Compare multiple VCF files in a folder to Excel marker file')
    parser.add_argument('--vcf_folder', required=True, help='Path to folder containing VCF files')
    parser.add_argument('--excel_file', required=True, help='Path to Excel marker file')
    parser.add_argument('--output_file', required=True, help='Path to output Excel file')
    args = parser.parse_args()

    vcf_folder = args.vcf_folder
    excel_file = args.excel_file
    output_file = args.output_file

    # Collect all .vcf files in the specified folder as full paths
    vcf_files = [os.path.join(vcf_folder, f) for f in os.listdir(vcf_folder) if f.lower().endswith('.vcf') and os.path.isfile(os.path.join(vcf_folder, f))]

    try:
        logger.info("Starting marker comparison across VCF files...")
        combined_core_df, combined_backup_df, combined_xy_df, matched_markers = compare_markers_multi(vcf_files, excel_file)
        logger.info("Marker comparison completed successfully.")

        logger.info("Saving results to output file...")
        save_results(
            combined_core_df,
            combined_backup_df,
            combined_xy_df,
            matched_markers,
            output_file
        )
        logger.info(f"Results saved successfully to {output_file}")

    except Exception as e:
        logger.error(f"An error occurred during processing: {e}")
        logger.error(traceback.format_exc())
        sys.exit(1)

    # Suppress pandas warnings for cleaner output
    warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)


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
        required_columns = ['CHROM', 'POS', 'IDs',  'Reference Allele', 'Variant Allele']
        
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
        
        #convert excel_data into a map of {pos,[chrom, ids, ref_allele, var_allele]}
        coreMap = {} #first sheet
        backupMap = {} #second sheet
        xyMap = {} #third sheet
        allIds = set() # to track all IDs across sheets
        for sheet_name, df in excel_data.items(): # creates map {pos: [chrom, ids, ref_allele, var_allele]}}
            logger.info(f"Processing sheet '{sheet_name}' with {len(df)} rows")
            allIds.update(df['IDs'].dropna().unique())  # Collect all IDs from this sheet
            if sheet_name == 'ISAG Core':
                coreMap = df.set_index('POS').T.to_dict('list')
            elif sheet_name == 'ISAG Back up':
                backupMap = df.set_index('POS').T.to_dict('list')
            elif sheet_name == 'X_Y':
                xyMap = df.set_index('POS').T.to_dict('list')
        logger.info(f"Excel data loaded successfully with {len(coreMap)} core markers, {len(backupMap)} backup markers, and {len(xyMap)} X_Y markers")
        return coreMap, backupMap, xyMap, allIds
    
    except FileNotFoundError:
        logger.error(f"Excel file not found: {excel_path}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error reading Excel file: {e}")
        logger.error(traceback.format_exc())
        sys.exit(1)

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

def compare_markers_multi(vcf_folderPath, excel_filePath):
    """
    Marker comparison: for each marker, query each line in each VCF for any variants in the excel.
    """
    logger.info("Starting multi-VCF marker comparison (per-marker VCF query)...")
    
    #get excel data
    coreMap, backupMap, xyMap, allVariants = read_excel_file(excel_filePath)
    
    # Store the marker results per VCF
    vcf_coreMap = {}  # Map to store VCF core markers
    vcf_backupMap = {}  # Map to store VCF backup markers
    vcf_xyMap = {}  # Map to store VCF X_Y markers
    
    # Keep track of which excel markers got matched
    matched_markers = {}  # marker name, vcf matched in
    
    # for each vcf, read vcf line by line
    for vcf_path in vcf_folderPath:
        vcf_filename = os.path.basename(vcf_path)
        logger.info(f"Processing VCF file: {vcf_filename}")
        
        # Read VCF header to get sample names
        sample_names = []
        
        # Prepare DataFrames for results of single vcf
        core_df = pd.DataFrame(columns=['ID', 'Excel Chrom', 'Excel Pos', 'Excel Alt', 'Excel Ref',
                                        '{vcf_filename} Chrom', '{vcf_filename} Pos', '{vcf_filename} Ref', '{vcf_filename} Alt',
                                        '{vcf_filename} homo ref', '{vcf_filename} homo alt', '{vcf_filename} hets',
                                        '{vcf_filename} Ref Allele Freq', '{vcf_filename} Alt Allele Freq'])  # DataFrame to store core results
        backup_df = pd.DataFrame(columns=['ID', 'Excel Chrom', 'Excel Pos', 'Excel Alt', 'Excel Ref',
                                        '{vcf_filename} Chrom', '{vcf_filename} Pos', '{vcf_filename} Ref', '{vcf_filename} Alt',
                                        '{vcf_filename} homo ref', '{vcf_filename} homo alt', '{vcf_filename} hets',
                                        '{vcf_filename} Ref Allele Freq', '{vcf_filename} Alt Allele Freq'])
        xy_df = pd.DataFrame(columns=['ID', 'Excel Chrom', 'Excel Pos', 'Excel Alt', 'Excel Ref',
                                        '{vcf_filename} Chrom', '{vcf_filename} Pos', '{vcf_filename} Ref', '{vcf_filename} Alt',
                                        '{vcf_filename} homo ref', '{vcf_filename} homo alt', '{vcf_filename} hets',
                                        '{vcf_filename} Ref Allele Freq', '{vcf_filename} Alt Allele Freq'])
        
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    fields = line.strip().split('\t')
                    if len(fields) >= 10:
                        sample_names = fields[9:]
                    break
                if line.startswith('#'):
                    continue # Skip header lines
                else: #not a header, compare markers
                    fields = line.strip().split('\t')
                    if len(fields) >= 10:
                        markerInExcel = ''
                        pos = int(fields[1])
                        if pos in coreMap:
                            compareList = calculateComparison(fields, coreMap)
                            if not compareList:
                                core_df.loc[len(core_df)] = compareList
                                markerInExcel = coreMap.get(pos)[1]  # Get marker ID from coreMap
                        elif pos in backupMap:
                            compareList = calculateComparison(fields, backupMap)
                            if not compareList:
                                backup_df.loc[len(backup_df)] = compareList
                                markerInExcel = backupMap.get(pos)[1]  # Get marker ID from backupMap
                        elif pos in xyMap:
                            compareList = calculateComparison(fields, xyMap)
                            if not compareList:
                                xy_df.loc[len(xy_df)] = compareList
                                markerInExcel = xyMap.get(pos)[1]  # Get marker ID from xyMap
                        if markerInExcel != '':
                            if markerInExcel not in matched_markers:
                                matched_markers[markerInExcel] = [] # Initialize list for this marker
                            matched_markers[markerInExcel].append(vcf_filename)  # Append filename for each match
            
        vcf_coreMap[vcf_filename] = core_df
        vcf_backupMap[vcf_filename] = backup_df
        vcf_xyMap[vcf_filename] = xy_df
        
    # check if any variants were missed
    for marker in allVariants:
        if marker not in matched_markers:
            matched_markers[marker] = []
        
    # Combine all VCF DataFrames into a single DataFrame for each category
    combined_core_df = combine_map_dfs(vcf_coreMap)
    combined_backup_df = combine_map_dfs(vcf_backupMap)
    combined_xy_df = combine_map_dfs(vcf_xyMap)
    
    return combined_core_df, combined_backup_df, combined_xy_df, matched_markers
                  
def combine_map_dfs(vcf_map):
    """
    Combine all DataFrames in a VCF map into a single DataFrame using outer join.
    After combining, remove columns matching '{vcf_filename} homo ref', '{vcf_filename} homo alt', and '{vcf_filename} hets'
    for every VCF in the map, keeping marker columns and all other VCF-specific columns.
    """
    dfs = list(vcf_map.values())
    if not dfs:
        return pd.DataFrame()
    # Start with the first DataFrame
    combined_df = dfs[0]
    join_cols = ['ID', 'Excel Chrom', 'Excel Pos', 'Excel Alt', 'Excel Ref']
    for df in dfs[1:]:
        combined_df = pd.merge(
            combined_df,
            df,
            on=join_cols,
            how='outer',
            suffixes=('', '_dup')
        )
        # Remove any duplicate columns from merge
        dup_cols = [col for col in combined_df.columns if col.endswith('_dup')]
        combined_df = combined_df.drop(columns=dup_cols)

    return combined_df

def calculateComparison(fields, sheetMap):
    pos = int(fields[1])
    chrom = standardize_chromosome(fields[0])
    if chrom != sheetMap.get(pos): #check chrom matches
        return 
    
    excel_list = sheetMap.get(pos)
    ref = fields[3]
    alt = fields[4]
    
    homo_ref = 0
    homo_alt = 0
    het = 0
    total_allele = 0
    missing_genotypes = 0
    for i in range(9,len(fields)):
        total_allele +=1
        g = fields[i].split(':')[0]
        if g == '0/0':
            homo_ref +=1
            total_allele += 2
        elif g == '0/1' or fields[i] == '1/0':
            het += 1
            total_allele += 2
        elif g == '1/1':
            homo_alt += 1
            total_allele += 2
        else:
            missing_genotypes += 1
            
    #check alleles match??
    
    # calculate freq
    afreq = (homo_alt * 2 + het)/total_allele
    rfreq = (homo_ref * 2 + het)/total_allele
            
    return [
        excel_list[1], pos, excel_list[0], excel_list[2], excel_list[3],
        chrom, pos, ref, alt,
        homo_ref, homo_alt, het, rfreq, afreq
        ]

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

# Function to save results to Excel with error handling and logging
def save_results(combined_core_df, combined_backup_df, combined_xy_df, matched_markers, output_file):
    """
    Save results to an Excel file with multiple sheets using pandas ExcelWriter.

    Parameters:
        combined_core_df (pd.DataFrame): DataFrame for core results.
        combined_backup_df (pd.DataFrame): DataFrame for backup results.
        combined_xy_df (pd.DataFrame): DataFrame for XY results.
        matched_markers (dict): Mapping of marker IDs to lists of matched VCF filenames.
        output_file (str): Path to the output Excel file.

    The function writes:
        - combined_core_df to "Core Results"
        - combined_backup_df to "Backup Results"
        - combined_xy_df to "XY Results"
        - matched_markers (as DataFrame) to "Matched Markers"
    Includes error handling and logs success or failure.
    """
    try:
        # Convert matched_markers dict to DataFrame
        matched_df = pd.DataFrame([
            {"Marker": marker, "Matched VCFs": ", ".join(sorted(set(vcf_list)))}
            for marker, vcf_list in matched_markers.items()
        ])

        # Write all DataFrames to separate sheets in the Excel file
        with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
            combined_core_df.to_excel(writer, sheet_name="Core Results", index=False)
            combined_backup_df.to_excel(writer, sheet_name="Backup Results", index=False)
            combined_xy_df.to_excel(writer, sheet_name="XY Results", index=False)
            matched_df.to_excel(writer, sheet_name="Matched Markers", index=False)
        logger.info(f"Results successfully saved to {output_file}")
    except Exception as e:
        logger.error(f"Failed to save results to {output_file}: {e}")
        logger.error(traceback.format_exc())
        raise

if __name__ == "__main__":
    main()
