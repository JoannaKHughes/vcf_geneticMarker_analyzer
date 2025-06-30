# Genetic Marker Analyzer

A Python script for analyzing genetic markers from Excel workbooks against large collections of VCF files. This tool efficiently processes 200+ VCF files (totaling 200+ million KB) to identify exact matches, proximity matches, and calculate allele frequencies.

## Features

- **Multi-sheet Excel processing**: Handles "ISAG Core", "ISAG Backup", and "X_Y" sheets
- **Large-scale VCF analysis**: Efficiently processes hundreds of VCF files
- **Exact and proximity matching**: Finds exact matches and variants within ±100 base pairs
- **Comprehensive frequency calculations**: Uses precise allele frequency formulas
- **Detailed output reports**: Generates multiple Excel sheets with analysis results
- **Progress tracking**: Real-time progress bars for long-running operations
- **Robust error handling**: Handles malformed files and missing data gracefully
- **Memory-efficient processing**: Single-pass VCF reading with smart caching

## Requirements

- Python 3.7+
- Required packages (install via `pip install -r requirements.txt`):
  - pandas>=1.5.0
  - numpy>=1.21.0
  - openpyxl>=3.0.9
  - tqdm>=4.64.0

## Installation

1. Clone or download the script files
2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

```bash
python genetic_marker_analyzer.py --excel_file <excel_file> --vcf_folder <vcf_folder> --output_file <output_file>
```

### Examples

```bash
# Basic analysis
python genetic_marker_analyzer.py --excel_file markers.xlsx --vcf_folder ./vcf_data --output_file results.xlsx

# With verbose logging
python genetic_marker_analyzer.py --excel_file markers.xlsx --vcf_folder ./vcf_data --output_file results.xlsx --verbose

# Custom log file
python genetic_marker_analyzer.py --excel_file markers.xlsx --vcf_folder ./vcf_data --output_file results.xlsx --log_file analysis.log

# Real-world example with Windows paths
python genetic_marker_analyzer.py --vcf_folder "C:\Users\joanna\Documents\Python_VSCode\Pangenome\tests\mini_vcf" --excel_file "C:\Users\joanna\Documents\Python_VSCode\Pangenome\ISAG_Marker_List.xlsx" --output_file "C:\Users\joanna\Documents\Python_VSCode\Pangenome\output\TEST_marker_comparison_miniVCF.xlsx"
```

### Command Line Arguments

**Required Arguments:**
- `--excel_file`: Path to Excel workbook with genetic marker data
- `--vcf_folder`: Path to folder containing VCF files
- `--output_file`: Path for output Excel file

**Optional Arguments:**
- `--verbose, -v`: Enable verbose logging
- `--log_file`: Specify custom log file path (default: genetic_marker_analysis.log)

## Input File Formats

### Excel Workbook Structure

The Excel file must contain three sheets with the following columns:

**Required Sheets:**
- "ISAG Core"
- "ISAG Backup" 
- "X_Y"

**Required Columns:**
- `CHROM`: Chromosome identifier
- `POS`: Position (numeric)
- `IDs`: Marker identifier
- `List`: List/category information
- `Reference Allele`: Reference allele
- `Variant Allele`: Variant allele
- `Strand`: Strand information
- `Sequence`: Sequence data

### VCF Files

Standard VCF format files (.vcf or .vcf.gz) with columns:
- `#CHROM`: Chromosome
- `POS`: Position
- `ID`: Variant ID
- `REF`: Reference allele
- `ALT`: Alternative allele
- `QUAL`: Quality score
- `FILTER`: Filter status
- `INFO`: Additional information
- `FORMAT`: Format specification
- Sample columns (starting from column 10)

## Output Structure

The output Excel file contains multiple sheets:

### Comparison Sheets
- `ISAG_Core_Results`
- `ISAG_Backup_Results`
- `X_Y_Results`

**Columns include:**
- Marker information (ID, CHROM, POS, REF, VAR, Strand, List)
- Per-VCF match results (`[identifier]_CHROM_exact`, `[identifier]_POS_exact`, `[identifier]_Match`)
- `Match_All`: Boolean indicating if marker found in ALL VCF files
- Frequency data (`freq_ref_of_[identifier]`, `freq_alt_of_[identifier]`)
- Combined frequencies (`freq_ref_of_all_vcf`, `freq_alt_of_all_vcf`)

### Summary Sheet
- Total markers analyzed
- Total VCF files processed
- Match statistics
- VCF file information (filename, identifier, file size)

### Missing_markers Sheet
- Markers missing from any VCF file
- Lists which VCFs contain/don't contain each marker
- Proximity match information

### Allele_representation Sheet
- Allele frequency analysis
- Unusual frequency pattern flagging
- REF/ALT allele representation across all VCFs

## Matching Logic

### Exact Matching
- CHROM and POS must match exactly between Excel and VCF data
- REF and ALT alleles are also verified

### Proximity Matching
- For unmatched markers, searches for variants within ±100 base pairs
- Reports closest matches with distance information

### Match_All Logic
- Marker must be found in ALL VCF files to be considered "Match_All"
- Missing from even one VCF results in Match_All = False

## Frequency Calculation

Uses the exact formula specified in requirements:

```python
for i in range(9, len(fields)):  # Sample data starts at column 9
    total_allele += 1
    g = fields[i].split(':')[0]
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

afreq = (homo_alt * 2 + het) / total_allele  # Alt allele frequency
rfreq = (homo_ref * 2 + het) / total_allele  # Ref allele frequency
```

## VCF File Identification

- Uses first 10 characters of filename as identifier
- Example: "minihorse1" for file "minihorse101s.vcf"
- Each VCF file tracked separately by its 10-character identifier

## Performance Optimizations

- **Single-pass processing**: Each VCF file read only once
- **Efficient data structures**: Uses dictionaries and sets for fast lookups
- **Memory management**: Smart caching with sequential processing
- **Progress tracking**: Visual progress bars for user feedback
- **Error resilience**: Continues processing despite individual file errors

## Logging

The script provides comprehensive logging:
- Analysis progress and statistics
- File processing status
- Error handling and warnings
- Performance metrics
- Final summary statistics

Log files include timestamps and severity levels for debugging large dataset issues.

## Error Handling

- **Input validation**: Checks for file existence and format
- **Malformed data**: Handles invalid VCF lines and Excel data
- **Missing files**: Graceful handling of inaccessible files
- **Memory constraints**: Efficient processing for large datasets
- **Data validation**: Validates required columns and data types

## Troubleshooting

### Common Issues

1. **Missing Excel sheets**: Ensure workbook contains required sheets with exact names
2. **VCF format errors**: Check VCF files are properly formatted
3. **Memory issues**: For very large datasets, consider processing in smaller batches
4. **File permissions**: Ensure read access to input files and write access to output location

### Performance Tips

- Use SSD storage for faster file I/O
- Ensure sufficient RAM for large VCF collections
- Close other applications to free system resources
- Consider processing subsets for initial validation

## Technical Details

### Architecture
- Object-oriented design with main `GeneticMarkerAnalyzer` class
- Modular functions for each processing step
- Efficient data structures for fast marker lookups
- Comprehensive error handling and logging

### Dependencies
- **pandas**: Excel I/O and data manipulation
- **numpy**: Numerical operations
- **openpyxl**: Excel file handling
- **tqdm**: Progress bar display
- **Standard library**: argparse, logging, pathlib, typing, collections, time, gzip

## License

This script is provided as-is for genetic marker analysis purposes. Please ensure appropriate citations and permissions when using with research data.

## Support

For technical issues or questions:
1. Check the log file for detailed error information
2. Verify input file formats match requirements
3. Ensure all dependencies are properly installed
4. Review the troubleshooting section above

## Version History

- **1.0**: Initial release with full functionality
  - Multi-sheet Excel processing
  - Large-scale VCF analysis
  - Comprehensive output reporting
  - Performance optimizations for large datasets