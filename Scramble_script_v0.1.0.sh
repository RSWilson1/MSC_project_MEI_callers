#!/bin/bash

# Set up for CWIC
#dx run --instance-type mem2_ssd1_v2_x8 --cost-limit 4 app-cwic --ssh -icredentials=project-GPgpg404jB0F71qxYqV6B980:file-GX9GKf84jB01gK7zPKQQ7kk4 -iimage="rswilson1/melt:latest" -iproject_mount_options="-readWrite" -y

set -exo pipefail

# Change to the root directory
cd ..
ls
# Unzip reference and create index
gunzip -c /project/003_230901_MSc_MEI_detection/references/hg19/hs37d5.fa.gz > scratch/hs37d5.fa
samtools faidx scratch/hs37d5.fa -o scratch/hs37d5.fa.fai

# Download files using aria2
apt install time
/usr/bin/time --verbose apt-get install -y aria2 # fix this by adding this to the docker.
cd scratch/
/usr/bin/time --verbose aria2c -x 16 -s 16 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/"$1"/alignment/"$1".mapped.ILLUMINA.bwa.ACB.low_coverage.20120522.bam
/usr/bin/time --verbose aria2c -x 16 -s 16 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/"$1"/alignment/"$1".mapped.ILLUMINA.bwa.ACB.low_coverage.20120522.bam.bai


# Run MELT
cd ..

/usr/bin/time --verbose java -jar ./opt/MELTv2.2.2/MELT.jar Single -a -b hs37d5/NC_007605 \
  -h scratch/hs37d5.fa \
  -bamfile scratch/"$1".mapped.ILLUMINA.bwa.ACB.low_coverage.20120522.bam \
  -n ./opt/MELTv2.2.2/add_bed_files/1KGP_Hg19/hg19.genes.bed \
  -t project/003_MSc_MEI_detection/test_MELT_docker/mei_list_complete_dockerenv.txt \
  -w scratch/

# Move output files to the output directory
mkdir -p project/003_MSc_MEI_detection/output/"$1"
mv scratch/*.vcf project/003_MSc_MEI_detection/output/"$1"/

# save the output files

dx-save-project

# Cleanup downloaded files
rm scratch/"$1".mapped.ILLUMINA.bwa.ACB.low_coverage.20120522.bam
rm scratch/"$1".mapped.ILLUMINA.bwa.ACB.low_coverage.20120522.bam.bai

# Then run: bash run_melt.sh HG01879
