# 100PioneerHorses_ISAG.py: Marker Comparison Pipeline

## Overview

This program compares genetic marker data from an Excel file (with multiple sheets) against variant data from one or more VCF files for 100 pioneer horses. It produces a comprehensive Excel report with detailed results, per-sheet statistics, unmatched marker analysis, and metadata for downstream analysis.

---

## Workflow

1. Reads all VCF files in a specified folder (each may contain multiple horse samples).
2. Reads an Excel marker file with multiple sheets (e.g., ISAG Core, ISAG Back up, X_Y).
3. Optionally reads a sample mapping file to relate VCF sample IDs to 'Finno Number' and ISAG CT 2025 IDs.
4. Standardizes chromosome names for comparison.
5. Compares markers from each Excel sheet to all VCFs, checking for exact matches by chromosome and position.
6. Generates a multi-sheet Excel output with results, statistics, unmatched marker details, and metadata.

---

## Input Files

- **VCF files:** All `.vcf` files in the specified folder. Each file contains variant data for one or more horse samples.
- **Excel marker file:** Contains marker information in multiple sheets (e.g., `ISAG_Marker_List.xlsx` with sheets like "ISAG Core", "ISAG Back up", "X_Y").
- **Sample mapping file (optional):** Maps VCF sample IDs to 'Finno Number' and 'ISAG CT 2025 ID' (e.g., `ISAG_2025_CT_Samples_from_Pioneer100.xlsx`).

---

## Chromosome Standardization

- Chromosome names are standardized for matching:
  - `'eMSYv3'` is mapped to `'Y'`
  - `'chr'`, `'Chr'`, `'CHR'` prefixes are removed (e.g., `'chr1'` → `'1'`)
  - Markers on the Y chromosome may not match due to the reference being female

---

## Output Excel File

The main output is an Excel file with the following sheets:

### 1. `[SheetName]_results`
- Detailed comparison results for each marker in each input Excel sheet.
- Columns include: marker ID, chromosome, position, alleles, and for each VCF: `<VCFNAME>_CHROM_exact`, `<VCFNAME>_POS_exact`, `<VCFNAME>_Match`, plus a "Match All" column.

### 2. `Summary`
- For each input sheet, includes:
  - `ISAG_Sheet_Name`
  - `Total_Markers`
  - `Non_Y_Markers`
  - `Non_Y_Matches to VCF1`
  - `Non_Y_Matches to VCF2`
  - `Non_Y_Success_Rate of VCF1`
  - `Non_Y_Success_Rate of VCF2`
  - `Overlapping success rate` (percent of non-Y markers matched in all VCFs)
- The number of VCF columns will match the number of VCFs compared.

### 3. `Unmatched_Markers`
- Markers from the Excel file that did not have an exact match in all VCFs.
- For each unmatched marker: sheet name, marker ID, chromosome, position, whether a nearby variant (within 100bp) exists in the first VCF, and a column `Matched_in_VCFs` listing the VCF(s) where the marker is matched (comma-separated, or "None" if not matched in any).

### 4. `Error_Report`
- Any errors encountered during processing (e.g., invalid positions, missing data).

### 5. `Analysis_Metadata`
- Metadata about the analysis: date, script version, number of variants, number of sheets, marker counts, and overall statistics.

### 6. `Sample_Mapping`
- Lists all VCF sample names and, if available, the corresponding 'Finno Number' and 'ISAG CT 2025 ID' from the mapping file.

### 7. `README`
- This documentation sheet, included for reference in the Excel output.

---

## Example Usage

```bash
python 100PioneerHorses_ISAG.py --vcf_folder /path/to/vcf_folder --excel_file /path/to/marker_file.xlsx --output_file /path/to/output.xlsx --mapping_file /path/to/mapping_file.xlsx
```

---

## Notes and Processing Details

- Unmatched markers are checked for nearby variants (within 100bp) in the first VCF to help identify potential close matches.
- The `Matched_in_VCFs` column in `Unmatched_Markers` lists which VCF(s) (by file name) the marker is matched in, or "None" if not matched in any.
- All processing steps are logged, and errors are reported in the `Error_Report` sheet.
- The script is designed for flexibility and can be adapted for other marker lists or VCF files with similar structure.
- Uses 'Finno Number' for matching VCF sample IDs in the mapping file.

---

## Contact

For questions or issues regarding this analysis or script, contact the program author or the VGL bioinformatics team.