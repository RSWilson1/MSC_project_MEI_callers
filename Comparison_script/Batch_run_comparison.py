import subprocess

# Define the dictionary containing IDs and countries

links_dict = {
    "HG00096": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam",
    "HG00268": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00268/alignment/HG00268.mapped.ILLUMINA.bwa.FIN.low_coverage.20130415.bam",
    "HG00419": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00419/alignment/HG00419.mapped.ILLUMINA.bwa.CHS.low_coverage.20130415.bam",
    "HG00759": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00759/alignment/HG00759.mapped.ILLUMINA.bwa.CDX.low_coverage.20130415.bam",
    "HG01051": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01051/alignment/HG01051.mapped.ILLUMINA.bwa.PUR.low_coverage.20120522.bam",
    "HG01112": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01112/alignment/HG01112.mapped.ILLUMINA.bwa.CLM.low_coverage.20120522.bam",
    "HG01500": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01500/alignment/HG01500.mapped.ILLUMINA.bwa.IBS.low_coverage.20120522.bam",
    "HG01565": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01565/alignment/HG01565.mapped.ILLUMINA.bwa.PEL.low_coverage.20120522.bam",
    "HG01583": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01583/alignment/HG01583.mapped.ILLUMINA.bwa.PJL.low_coverage.20130415.bam",
    "HG01595": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01595/alignment/HG01595.mapped.ILLUMINA.bwa.KHV.low_coverage.20120522.bam",
    "HG01879": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01879/alignment/HG01879.mapped.ILLUMINA.bwa.ACB.low_coverage.20120522.bam",
    "HG02568": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02568/alignment/HG02568.mapped.ILLUMINA.bwa.GWD.low_coverage.20130415.bam",
    "HG02922": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02922/alignment/HG02922.mapped.ILLUMINA.bwa.ESN.low_coverage.20130415.bam",
    "HG03006": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG03006/alignment/HG03006.mapped.ILLUMINA.bwa.BEB.low_coverage.20130415.bam",
    "HG03052": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG03052/alignment/HG03052.mapped.ILLUMINA.bwa.MSL.low_coverage.20121211.bam",
    "HG03642": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG03642/alignment/HG03642.mapped.ILLUMINA.bwa.STU.low_coverage.20130415.bam",
    "HG03742": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG03742/alignment/HG03742.mapped.ILLUMINA.bwa.ITU.low_coverage.20121211.bam",
    "NA18525": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18525/alignment/NA18525.mapped.ILLUMINA.bwa.CHB.low_coverage.20130415.bam",
    "NA18939": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18939/alignment/NA18939.mapped.ILLUMINA.bwa.JPT.low_coverage.20130415.bam",
    "NA19017": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19017/alignment/NA19017.mapped.ILLUMINA.bwa.LWK.low_coverage.20121211.bam",
    "NA19625": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19625/alignment/NA19625.mapped.ILLUMINA.bwa.ASW.low_coverage.20120522.bam",
    "NA19648": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19648/alignment/NA19648.mapped.ILLUMINA.bwa.MXL.low_coverage.20120522.bam",
    "NA20502": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA20502/alignment/NA20502.mapped.ILLUMINA.bwa.TSI.low_coverage.20120522.bam",
    "NA20845": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA20845/alignment/NA20845.mapped.ILLUMINA.bwa.GIH.low_coverage.20120522.bam",
    "NA12878": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam",
    "NA19238": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19238/alignment/NA19238.mapped.ILLUMINA.bwa.YRI.low_coverage.20130415.bam",
    "NA19239": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19239/alignment/NA19239.mapped.ILLUMINA.bwa.YRI.low_coverage.20130415.bam",
    "NA19240": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19240/alignment/NA19240.mapped.ILLUMINA.bwa.YRI.low_coverage.20130415.bam",
}

# Loop over the dictionary
for ID, url in links_dict.items():
    log_filename = f"{ID}_comparsion_log.txt"  # Define the log file name based on ID
    split_url = url.split("/")  # Split the url by '/'
    filename = split_url[-1]  # Get the last element of the split url
    print(f"Processing {filename}")
    with open(log_filename, "w") as log_file:
        # Execute the bash script and redirect stdout and stderr to the log file
        subprocess.run(
            ["bash", "Match_SV_v3", str(ID), str(url), str(filename)],
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )
