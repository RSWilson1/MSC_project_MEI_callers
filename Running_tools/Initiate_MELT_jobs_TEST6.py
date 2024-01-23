import subprocess

# Define the dictionary containing IDs and countries
links_dict = {
    "HG00268": "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00268/alignment/HG00268.mapped.ILLUMINA.bwa.FIN.low_coverage.20130415.bam"
}

# Loop over the dictionary
for ID, url in links_dict.items():
    log_filename = f'{ID}_log.txt'  # Define the log file name based on ID
    split_url = url.split('/')  # Split the url by '/'
    filename = split_url[-1]  # Get the last element of the split url
    print(f'Processing {filename}')
    with open(log_filename, 'w') as log_file:
        # Execute the bash script and redirect stdout and stderr to the log file
        subprocess.run(['bash', 'MELT_script_v0.1.20_second_run.sh', str(ID), str(url), str(filename)], stdout=log_file, stderr=subprocess.STDOUT)
