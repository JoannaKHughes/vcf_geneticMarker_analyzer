#!/usr/bin/env python3
"""
Example usage script for Genetic Marker Analyzer

This script demonstrates how to use the GeneticMarkerAnalyzer class programmatically.
"""

from genetic_marker_analyzer import GeneticMarkerAnalyzer
import sys
import os

def run_example_analysis():
    """
    Example function showing how to use the analyzer programmatically.
    """
    
    # Example file paths (modify these for your actual data)
    excel_file = "example_markers.xlsx"  # Your Excel file with marker data
    vcf_folder = "./vcf_files"          # Folder containing VCF files
    output_file = "analysis_results.xlsx"  # Output Excel file
    
    print("Genetic Marker Analyzer - Example Usage")
    print("=" * 50)
    
    # Check if example files exist
    if not os.path.exists(excel_file):
        print(f"⚠️  Example Excel file '{excel_file}' not found.")
        print("   Please provide your own Excel file with the required format.")
        return False
    
    if not os.path.exists(vcf_folder):
        print(f"⚠️  VCF folder '{vcf_folder}' not found.")
        print("   Please provide a folder containing VCF files.")
        return False
    
    try:
        # Create analyzer instance
        print(f"📊 Initializing analyzer...")
        analyzer = GeneticMarkerAnalyzer(
            excel_path=excel_file,
            vcf_folder=vcf_folder,
            output_path=output_file
        )
        
        # Run the complete analysis
        print(f"🚀 Starting analysis...")
        analyzer.run_analysis()
        
        print(f"✅ Analysis completed successfully!")
        print(f"📄 Results saved to: {output_file}")
        print(f"📋 Check log file: genetic_marker_analysis.log")
        
        return True
        
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        return False

def create_sample_excel():
    """
    Create a sample Excel file with the required format for testing.
    """
    import pandas as pd
    
    # Sample data for each required sheet
    isag_core_data = {
        'CHROM': ['1', '2', '3', '1', '2'],
        'POS': [1000, 2000, 3000, 1500, 2500],
        'IDs': ['marker1', 'marker2', 'marker3', 'marker4', 'marker5'],
        'List': ['Core', 'Core', 'Core', 'Core', 'Core'],
        'Reference Allele': ['A', 'G', 'C', 'T', 'A'],
        'Variant Allele': ['T', 'C', 'G', 'A', 'G'],
        'Strand': ['+', '+', '-', '+', '-'],
        'Sequence': ['ATCG', 'GCTA', 'CGAT', 'TACG', 'AGCT']
    }
    
    isag_backup_data = {
        'CHROM': ['3', '4', '5'],
        'POS': [3500, 4000, 5000],
        'IDs': ['backup1', 'backup2', 'backup3'],
        'List': ['Backup', 'Backup', 'Backup'],
        'Reference Allele': ['G', 'A', 'C'],
        'Variant Allele': ['A', 'T', 'T'],
        'Strand': ['+', '-', '+'],
        'Sequence': ['GATC', 'ATCG', 'CTAG']
    }
    
    x_y_data = {
        'CHROM': ['X', 'Y', 'X'],
        'POS': [1000, 2000, 3000],
        'IDs': ['X_marker1', 'Y_marker1', 'X_marker2'],
        'List': ['X_Y', 'X_Y', 'X_Y'],
        'Reference Allele': ['A', 'G', 'C'],
        'Variant Allele': ['G', 'A', 'T'],
        'Strand': ['+', '+', '-'],
        'Sequence': ['AGCT', 'GCAT', 'CTGA']
    }
    
    # Create DataFrames
    isag_core_df = pd.DataFrame(isag_core_data)
    isag_backup_df = pd.DataFrame(isag_backup_data)
    x_y_df = pd.DataFrame(x_y_data)
    
    # Write to Excel file
    sample_file = "sample_markers.xlsx"
    with pd.ExcelWriter(sample_file, engine='openpyxl') as writer:
        isag_core_df.to_excel(writer, sheet_name='ISAG Core', index=False)
        isag_backup_df.to_excel(writer, sheet_name='ISAG Backup', index=False)
        x_y_df.to_excel(writer, sheet_name='X_Y', index=False)
    
    print(f"📄 Sample Excel file created: {sample_file}")
    return sample_file

def create_sample_vcf():
    """
    Create a sample VCF file for testing.
    """
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
1	1000	.	A	T	60	PASS	AF=0.5	GT	0/0	0/1	1/1
2	2000	.	G	C	60	PASS	AF=0.3	GT	0/0	0/0	0/1
1	1500	.	T	A	60	PASS	AF=0.7	GT	1/1	0/1	1/1
3	3000	.	C	G	60	PASS	AF=0.4	GT	0/1	0/0	0/1
"""
    
    # Create VCF folder if it doesn't exist
    os.makedirs("sample_vcf_files", exist_ok=True)
    
    # Create sample VCF files
    vcf_files = ["sample001.vcf", "sample002.vcf"]
    
    for vcf_file in vcf_files:
        vcf_path = os.path.join("sample_vcf_files", vcf_file)
        with open(vcf_path, 'w') as f:
            f.write(vcf_content)
        print(f"📄 Sample VCF file created: {vcf_path}")
    
    return "sample_vcf_files"

def main():
    """
    Main function to demonstrate different usage scenarios.
    """
    
    print("Genetic Marker Analyzer - Example Usage")
    print("=" * 50)
    print()
    
    # Check command line arguments
    if len(sys.argv) > 1:
        if sys.argv[1] == "--create-samples":
            print("Creating sample data files...")
            excel_file = create_sample_excel()
            vcf_folder = create_sample_vcf()
            
            print(f"\n🎯 Sample files created!")
            print(f"   Excel file: {excel_file}")
            print(f"   VCF folder: {vcf_folder}")
            print(f"\nTo run analysis with sample data:")
            print(f"   python genetic_marker_analyzer.py {excel_file} {vcf_folder} sample_results.xlsx")
            return
    
    # Show usage examples
    print("Usage Examples:")
    print("-" * 20)
    print()
    
    print("1. Basic Analysis:")
    print("   python genetic_marker_analyzer.py markers.xlsx ./vcf_data results.xlsx")
    print()
    
    print("2. With Verbose Logging:")
    print("   python genetic_marker_analyzer.py markers.xlsx ./vcf_data results.xlsx --verbose")
    print()
    
    print("3. Create Sample Data:")
    print("   python example_usage.py --create-samples")
    print()
    
    print("4. Programmatic Usage:")
    print("   See the run_example_analysis() function in this file")
    print()
    
    # Try to run example analysis if files exist
    print("Checking for example data files...")
    
    # Look for any Excel files and VCF folders
    excel_files = [f for f in os.listdir('.') if f.endswith('.xlsx')]
    vcf_folders = [d for d in os.listdir('.') if os.path.isdir(d) and 'vcf' in d.lower()]
    
    if excel_files and vcf_folders:
        print(f"Found Excel file: {excel_files[0]}")
        print(f"Found VCF folder: {vcf_folders[0]}")
        
        response = input("\nWould you like to run analysis with these files? (y/N): ")
        if response.lower() == 'y':
            # Update the analyzer with found files
            analyzer = GeneticMarkerAnalyzer(
                excel_path=excel_files[0],
                vcf_folder=vcf_folders[0],
                output_path="example_results.xlsx"
            )
            
            try:
                analyzer.run_analysis()
                print("✅ Analysis completed! Check example_results.xlsx")
            except Exception as e:
                print(f"❌ Analysis failed: {e}")
    else:
        print("No suitable data files found.")
        print("Use --create-samples to generate sample data files.")

if __name__ == "__main__":
    main()