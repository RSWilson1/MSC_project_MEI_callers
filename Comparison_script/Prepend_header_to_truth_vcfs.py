import os

header = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20240116
##reference=hs37d5
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=PCR_confirmed,Number=1,Type=Integer,Description="Confirmed MEI by PCR">
##INFO=<ID=FOUND,Number=1,Type=Integer,Description="FOUND">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele Depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth Used For Variant Calling">
##FORMAT=<ID=VF,Number=.,Type=Float,Description="Variant Frequency">
##FORMAT=<ID=NL,Number=1,Type=Integer,Description="Applied BaseCall Noise Level">
##FORMAT=<ID=SB,Number=1,Type=Float,Description="StrandBias Score">
##FORMAT=<ID=NC,Number=1,Type=Float,Description="Fraction of bases which were uncalled or with basecall quality below the minimum threshold">
##FORMAT=<ID=US,Number=.,Type=Integer,Description="Supporting read type counts">
##FORMAT=<ID=AQ,Number=1,Type=Float,Description="Variant artifact adjusted quality score">
##FORMAT=<ID=LQ,Number=1,Type=Float,Description="Likelihood ratio quality score">
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##contig=<ID=MT,length=16569>
##contig=<ID=GL000207.1,length=4262>
##contig=<ID=GL000226.1,length=15008>
##contig=<ID=GL000229.1,length=19913>
##contig=<ID=GL000231.1,length=27386>
##contig=<ID=GL000210.1,length=27682>
##contig=<ID=GL000239.1,length=33824>
##contig=<ID=GL000235.1,length=34474>
##contig=<ID=GL000201.1,length=36148>
##contig=<ID=GL000247.1,length=36422>
##contig=<ID=GL000245.1,length=36651>
##contig=<ID=GL000197.1,length=37175>
##contig=<ID=GL000203.1,length=37498>
##contig=<ID=GL000246.1,length=38154>
##contig=<ID=GL000249.1,length=38502>
##contig=<ID=GL000196.1,length=38914>
##contig=<ID=GL000248.1,length=39786>
##contig=<ID=GL000244.1,length=39929>
##contig=<ID=GL000238.1,length=39939>
##contig=<ID=GL000202.1,length=40103>
##contig=<ID=GL000234.1,length=40531>
##contig=<ID=GL000232.1,length=40652>
##contig=<ID=GL000206.1,length=41001>
##contig=<ID=GL000240.1,length=41933>
##contig=<ID=GL000236.1,length=41934>
##contig=<ID=GL000241.1,length=42152>
##contig=<ID=GL000243.1,length=43341>
##contig=<ID=GL000242.1,length=43523>
##contig=<ID=GL000230.1,length=43691>
##contig=<ID=GL000237.1,length=45867>
##contig=<ID=GL000233.1,length=45941>
##contig=<ID=GL000204.1,length=81310>
##contig=<ID=GL000198.1,length=90085>
##contig=<ID=GL000208.1,length=92689>
##contig=<ID=GL000191.1,length=106433>
##contig=<ID=GL000227.1,length=128374>
##contig=<ID=GL000228.1,length=129120>
##contig=<ID=GL000214.1,length=137718>
##contig=<ID=GL000221.1,length=155397>
##contig=<ID=GL000209.1,length=159169>
##contig=<ID=GL000218.1,length=161147>
##contig=<ID=GL000220.1,length=161802>
##contig=<ID=GL000213.1,length=164239>
##contig=<ID=GL000211.1,length=166566>
##contig=<ID=GL000199.1,length=169874>
##contig=<ID=GL000217.1,length=172149>
##contig=<ID=GL000216.1,length=172294>
##contig=<ID=GL000215.1,length=172545>
##contig=<ID=GL000205.1,length=174588>
##contig=<ID=GL000219.1,length=179198>
##contig=<ID=GL000224.1,length=179693>
##contig=<ID=GL000223.1,length=180455>
##contig=<ID=GL000195.1,length=182896>
##contig=<ID=GL000212.1,length=186858>
##contig=<ID=GL000222.1,length=186861>
##contig=<ID=GL000200.1,length=187035>
##contig=<ID=GL000193.1,length=189789>
##contig=<ID=GL000194.1,length=191469>
##contig=<ID=GL000225.1,length=211173>
##contig=<ID=GL000192.1,length=547496>
##contig=<ID=NC_007605,length=171823>
##contig=<ID=hs37d5,length=35477943>
"""

def add_header_to_files(directory_path, prefix="Truth"):
    for filename in os.listdir(directory_path):
        if filename.startswith(prefix) and filename.endswith(".vcf"):
            input_file_path = os.path.join(directory_path, filename)
            output_file_path = os.path.join(directory_path, f"valid_{filename}")

            with open(input_file_path, "r") as input_file, open(output_file_path, "w") as output_file:
                output_file.write(header)
                output_file.write(input_file.read())

if __name__ == "__main__":
    directory_path = "truth_vcfs/"
    add_header_to_files(directory_path)
