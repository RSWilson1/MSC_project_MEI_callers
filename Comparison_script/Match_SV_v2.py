import vcf
import argparse

def parse_args():
    """
    Parses command line arguments.
    Example usage:
    python Match_SV_v2.py -vcf_baseline /path/to/baseline.vcf -testvcf /path/to/test.vcf -range_limit 5
    """
    parser = argparse.ArgumentParser(description="Compare variants in two VCF files.")
    parser.add_argument("-vcf_baseline", required=True, help="Path to the baseline VCF file.")
    parser.add_argument("-testvcf", required=True, help="Path to the test VCF file.")
    parser.add_argument("-range_limit", type=int, default=5, help="Position range within which variants are considered similar.")

    arguments = parser.parse_args()

    return arguments


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


def compare_vcfs(vcf_baseline, testvcf, range_limit):
    vcf_reader_base = vcf.Reader(open(vcf_baseline, 'r'))
    vcf_reader_test = vcf.Reader(open(testvcf, 'r'))

    total_variants = 0
    shared_variants = 0

    chrom_variants1 = {}  # Dict to store variants for each chromosome in vcf_baseline

    # Group variants by chromosome in vcf_baseline
    for record1 in vcf_reader_base:
        total_variants += 1
        chrom = record1.CHROM
        if chrom.startswith("chr"):
            chrom = chrom[3:]  # Remove "chr" prefix
        if chrom not in chrom_variants1:
            chrom_variants1[chrom] = []
        chrom_variants1[chrom].append(record1)

    for record2 in vcf_reader_test:
        chrom2 = record2.CHROM
        pos2 = record2.POS

        if chrom2.startswith("chr"):
            chrom2 = chrom2[3:]  # Remove "chr" prefix

        if chrom2 in chrom_variants1:
            variants_in_chrom1 = chrom_variants1[chrom2]

            for record1 in variants_in_chrom1:
                pos1 = record1.POS

                if within_range(pos1, pos2, range_limit):
                    alt1 = record1.ALT[0]
                    alt2 = record2.ALT[0]

                    keyword1 = get_alu_line_sva(str(alt1))
                    keyword2 = get_alu_line_sva(str(alt2))

                    if keyword1 is not None and keyword2 is not None and keyword1 == keyword2:
                        print(record2)
                        shared_variants += 1
                        break  # Found a matching position, no need to continue


    shared_percentage = (shared_variants / total_variants) * 100 if total_variants > 0 else 0

    print(f"Summary Statistics:")
    print(f"Total variants in baseline VCF: {total_variants}")
    print(f"Variants shared between the two VCFs: {shared_variants}")
    print(f"Percentage of shared variants: {shared_percentage:.2f}%")

if __name__ == "__main__":
    args = parse_args()

    compare_vcfs(args.vcf_baseline, args.testvcf, args.range_limit)