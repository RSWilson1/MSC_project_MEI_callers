"""
For MELT 3 VCFs are created theses need to be merged into one VCF using
bcftools concat and then bgzipped and tabix indexed.
"""
import subprocess

def main():
    """
    Main function to merge MELT VCFs.

    This function iterates over a list of IDs and performs the following steps for each ID:
    1. Changes the working directory to the MELT output directory for the ID.
    2. Concatenates the ALU.final.vcf, SVA.final.vcf,
        and LINE1.final.vcf files into a single VCF file.
    3. Compresses the merged VCF file using bgzip.
    4. Indexes the compressed VCF file using tabix.
    5. Prints a message indicating that the process for the ID is finished.
    """
    # id_list = ['HG00096', 'HG00268', 'HG00419', 'HG00759', 'HG01051', 'HG01112',
    #            'HG01500', 'HG01565', 'HG01583', 'HG01595', 'HG01879', 'HG02568',
    #            'HG02922', 'HG03006', 'HG03052', 'HG03642', 'HG03742', 'NA18525',
    #            'NA18939', 'NA19017', 'NA19625', 'NA19648', 'NA20502', 'NA20845',
    #            'NA12878', 'NA19238', 'NA19239', 'NA19240']

    id_list = ["HG00096","HG00268","HG00419","HG00759",
              "HG01051","HG01112","HG01500","HG01565"]

    for hg_id in id_list:
        # ALU.final.vcf, SVA.final.vcf, and LINE1.final.vcf are the inputs.
        file_path = f'/project/003_230901_MSc_MEI_detection/benchmarking_output/{hg_id}/MELT/'
        subprocess.run(['bcftools', 'concat', '-a',
                        f'{file_path}ALU.final_comp.vcf.gz',
                        f'{file_path}SVA.final_comp.vcf.gz',
                        f'{file_path}LINE1.final_comp.vcf.gz',
                        '-o', f'{file_path}{hg_id}_MELT_concat.vcf'])
        subprocess.run(['bgzip', f'{file_path}{hg_id}_MELT_concat.vcf'])
        print(f'{hg_id} Finished')

if __name__ == "__main__":
    main()
