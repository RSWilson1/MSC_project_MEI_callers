"""
This script compares two VCF files and reports the number of variants shared between them.
It bases this on different criteria:
 - Matches based on start position of variants.
 - +/- a specified overlap range (default is 50bp).
 - Matches based on the presence of ALU, LINE1 or SVA in the ALT field of the variants.

It also reports the percentage of variants shared between the two VCFs.
Returns a summary of the number of variants shared between the two VCFs.
CSV??
"""
#import vcf
from cyvcf2 import VCF
import argparse
import pandas as pd

LST_OF_CHRS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                        'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                        'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                        'chr20', 'chr21', 'chr22', 'chrX', 'chrY'
                        '1', '2', '3', '4', '5', '6', '7',
                        '8', '9', '10', '11', '12', '13',
                        '14', '15', '16', '17', '18', '19',
                        '20', '21', '22', 'X', 'Y']

def parse_args():
    """
    Parses command line arguments.
    Example usage:
    python Match_SV_v2.py -vcf_baseline /path/to/baseline.vcf -test_vcf /path/to/test.vcf -range_limit 50
    """
    parser = argparse.ArgumentParser(description="Compare variants in two VCF files.")
    parser.add_argument("-vcf_baseline", required=True, help="Path to the baseline VCF file.")
    parser.add_argument("-test_vcf", required=True, help="Path to the test VCF file.")
    parser.add_argument("-range_limit", type=int, default=50, help="Position range within which variants are considered similar.")

    arguments = parser.parse_args()

    return arguments


def get_MEI_caller(file_path):
    """
    Get the MEI caller from the file path.
    Get MEI caller, i.e. MELT, Scramble or Mobster.
    Or truth set. i.e. skip.

    Parameters
    ----------
    file_path : str
        The file path containing the MEI caller information.
    Returns
    -------
    str or None
        The MEI caller extracted from the file path, or None if not found.
    """
    MEI_callers = ['melt', 'scramble', 'mobster']
    file_path_lower = file_path.lower()
    print(file_path_lower)
    for caller in MEI_callers:
        if caller in file_path_lower:
            return caller
    return None


def get_alu_line_sva(alt):
    keywords = ['ALU', 'LINE1', 'SVA', 'L1']
    for keyword in keywords:
        if keyword in alt:
            if keyword == 'L1':
                keyword = 'LINE1'
            return keyword
    return None


def within_range(pos1, pos2, range_limit):
    return abs(pos1 - pos2) <= range_limit


def search_vcfs(vcf_baseline, test_vcf, svlen, range_limit):
    """
    Search for matching variants in the two VCF files.

    Parameters
    ----------
    vcf_baseline : vcf.Reader
        The baseline VCF file.
    test_vcf : vcf.Reader
        The test VCF file.
    svlen : boolean
        Whether to use svlen or not.
    range_limit : int
        The range limit within which variants are considered similar.
    """
    total_variants = 0
    shared_variants = 0

    chrom_variants1 = {}  # Dict to store variants for each chromosome in vcf_baseline
    shared_variants_vcf = []  # List to store shared variants

    # Group variants by chromosome in vcf_baseline
    for record1 in vcf_baseline:
        #total_variants += 1
        chrom = record1.CHROM
        if chrom not in LST_OF_CHRS:
            #print(f"rejected - {chrom}")
            continue
        #remove scramble dels
        id1 = record1.ID
        if id1 == "DEL":
            continue
        total_variants += 1
        if chrom.startswith("chr"):
            #print(chrom)
            #chrom = chrom[3:]  # Remove "chr" prefix # why not replace?
            chrom = chrom.replace("chr", "")
            record1.CHROM = chrom
        if chrom not in chrom_variants1:
            chrom_variants1[f"{chrom}"] = []
        chrom_variants1[f"{chrom}"].append(record1)

    print("Finished reading baseline VCF")

    for i, record2 in enumerate(test_vcf):
        print(f"Processing record {i}")
        chrom2 = record2.CHROM
        pos2 = record2.POS
        id2 = record2.ID
        if svlen:
            svlen2 = record2.INFO.get('SVLEN') #record2.INFO['SVLEN']
        elif not svlen:
            svlen2 = len(record2.ALT[0])
        else:
            svlen2 = None
            # Truth VCF contains no len attribute
            # print("only truth set")
        if id2 == "DEL":
            continue

        if chrom2 not in LST_OF_CHRS:
            continue
        if chrom2.startswith("chr"):
            #chrom2 = chrom2[3:]  # Remove "chr" prefix
            chrom2 = chrom2.replace("chr", "")
        if chrom2 in chrom_variants1:
            variants_in_chrom1 = chrom_variants1[chrom2]
            for record1 in variants_in_chrom1:
                pos1 = record1.POS
                #end1 = record1.INFO.get('END')
                if within_range(pos1, pos2, range_limit):
                    alt1 = record1.ALT[0]
                    alt2 = record2.ALT[0]

                    keyword1 = get_alu_line_sva(str(alt1))
                    keyword2 = get_alu_line_sva(str(alt2))

                    if keyword1 is not None and keyword2 is not None and keyword1 == keyword2:
                        shared_variants += 1
                        entry = (
                            f"{record2.CHROM} {record2.POS} {record2.ID}",
                            f"{record2.REF} {record2.ALT} {record2.QUAL}",
                            f"{record2.FILTER} {record2.INFO} {record2.FORMAT}"
                        )
                        shared_variants_vcf.append(entry)
                        break  # Found a matching position, no need to continue
                    else:
                        print("No keyword match found") #  {record2}
                else:
                    print("No pos match found") # {record2}
    return total_variants, shared_variants, shared_variants_vcf


def compare_vcfs(vcf_baseline, test_vcf, range_limit):
    #Check caller or truth
    caller_base = get_MEI_caller(vcf_baseline)
    caller_test = get_MEI_caller(test_vcf)
    if caller_base is None:
        caller_base = "truth"
    else:
        caller_base = caller_base
    #print(caller_base)
    if caller_test is None:
        caller_test = "truth"
    else:
        caller_test = caller_test
    #print(caller_test)
    if ("melt" in caller_base or caller_test) or ("mobster" in caller_base or caller_test):
        print("Checking with SVLEN")
        svlen = True
    elif "scramble" in caller_base or caller_test:
        svlen = False
        #extract svlen from len(ALT) - len(REF)
    else:
        print("only truth set")

    vcf_reader_base = VCF(vcf_baseline) #vcf.Reader(open(vcf_baseline, 'r'))
    vcf_reader_test = VCF(test_vcf) #vcf.Reader(open(test_vcf, 'r'))

    total_variants, shared_variants, shared_variants_vcf = search_vcfs(
        vcf_reader_base, vcf_reader_test, svlen, range_limit
        )

    shared_percentage = (shared_variants / total_variants) * 100 if total_variants > 0 else 0

    print("Summary Statistics:")
    print(f"Total variants in baseline VCF: {total_variants}")
    print(f"Variants shared between the two VCFs: {shared_variants}")
    print(f"Percentage of shared variants: {shared_percentage:.2f}%")
    #print(shared_variants_vcf)
    return shared_variants_vcf, shared_percentage, shared_variants, total_variants


def run_for_multiple_samples():
    """
    Run the comparison for multiple samples.
    """
    # Assuming you have a list of sample IDs
    sample_ids = ["HG00096", "HG00097", "HG00098"]

    # Assuming you have a list of tools
    tools = ["MELT", "scramble", "mobster"]

    # Run the comparison for each sample and tool
    results = []
    for sample in sample_ids:
        for tool in tools:
            # Get the path to the baseline VCF
            # Example names:
            #HG00096_MELT_concat.vcf.gz
            #HG00096_scramble.vcf
            #HG00096_mobster_predictions.vcf
            base_path = "/project/003_230901_MSc_MEI_detection/benchmarking_output"
            if tool == "MELT":
                test_vcf = f"{base_path}/{sample}/{tool}/{sample}_{tool}_concat.vcf.gz"
            elif tool == "scramble":
                test_vcf = f"{base_path}/{sample}/{tool}/{sample}_{tool}.vcf"
            elif tool == "mobster":
                test_vcf = f"{base_path}/{sample}/{tool}/{sample}_{tool}_predictions.vcf"
            else:
                print("No tool found")
                exit(1)

            vcf_baseline = f"{sample}_{tool}_baseline.vcf"
            # Run the comparison
            shared_variants_vcf, shared_percentage, shared_variants, total_variants = \
                compare_vcfs(vcf_baseline, test_vcf, 50)
            # Create a dictionary with the results
            result_dict = {
                "Sample_ID": sample,
                "Tool": tool,
                "Total_Variants": total_variants,
                "Shared_Variants": shared_variants,
                "Shared_Percentage": shared_percentage,
                "Shared_Variants_VCF": shared_variants_vcf
            }

            # Append the dictionary to the results list
            results.append(result_dict)
    # Create a DataFrame from the results list
    df = pd.DataFrame(results)

    # Write the DataFrame to a CSV file
    csv_filename = "results.csv"
    df.to_csv(csv_filename, index=False)

if __name__ == "__main__":
    args = parse_args()
    compare_vcfs(args.vcf_baseline, args.test_vcf, args.range_limit)
