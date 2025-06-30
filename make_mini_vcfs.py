# make_mini_vcfs.py
#
# Usage:
#   Run this script to generate "mini" VCF files for testing.
#   For each .vcf file in C:/Users/joanna/Documents/Python_VSCode/Pangenome/vcf,
#   a new file named mini_{original_filename} will be created in
#   C:/Users/joanna/Documents/Python_VSCode/Pangenome/tests containing:
#     - All header lines (lines starting with '#')
#     - The first 200 variant lines (lines not starting with '#')
#
# Requirements:
#   - Python 3.x
#   - No external dependencies

import os

# Input and output directories
VCF_DIR = r"C:/Users/joanna/Documents/Python_VSCode/Pangenome/vcf"
OUT_DIR = r"C:/Users/joanna/Documents/Python_VSCode/Pangenome/tests"

def make_mini_vcf(input_path, output_path, max_variants=200):
    """
    Create a mini VCF file with all header lines and the first `max_variants` variant lines.
    """
    header_lines = []
    variant_lines = []
    with open(input_path, "r", encoding="utf-8") as infile:
        for line in infile:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                if len(variant_lines) < max_variants:
                    variant_lines.append(line)
                # Stop reading variant lines after reaching the limit
                if len(variant_lines) >= max_variants:
                    break
        # If file is not sorted (headers not all at top), continue reading for more headers
        for line in infile:
            if line.startswith("#"):
                header_lines.append(line)
    # Write output file
    with open(output_path, "w", encoding="utf-8") as outfile:
        outfile.writelines(header_lines)
        outfile.writelines(variant_lines)

def main():
    """
    Main function to process all .vcf files in the input directory.
    """
    # Ensure output directory exists
    os.makedirs(OUT_DIR, exist_ok=True)

    # List all .vcf files in the input directory
    for filename in os.listdir(VCF_DIR):
        if filename.lower().endswith(".vcf"):
            input_path = os.path.join(VCF_DIR, filename)
            output_filename = f"mini_{filename}"
            output_path = os.path.join(OUT_DIR, output_filename)
            print(f"Processing {filename} -> {output_filename}")
            make_mini_vcf(input_path, output_path)

if __name__ == "__main__":
    main()