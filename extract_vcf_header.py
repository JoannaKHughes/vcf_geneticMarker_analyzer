# extract_vcf_header.py
# Tiny script to find the header line in a VCF file and write it to a text file

def extract_vcf_header(vcf_path, output_path):
    """
    Finds the header line (starting with #CHROM) in a VCF file and writes it to output_path.
    """
    with open(vcf_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            if line.startswith('#CHROM'):
                outfile.write(line)
                print(f"Header written to {output_path}")
                return
    print("No #CHROM header found in the VCF file.")

if __name__ == "__main__":
    vcf = "C:/Users/joanna/Documents/Python_VSCode/Pangenome/vcf/Pioneer100CEH_03212022.vcf"
    out = "C:/Users/joanna/Documents/Python_VSCode/Pangenome/output/vcf_header.txt"
    extract_vcf_header(vcf, out)