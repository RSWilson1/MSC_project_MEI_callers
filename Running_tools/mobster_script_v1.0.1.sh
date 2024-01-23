#!/bin/bash

# Set up for CWIC
#dx run --instance-type mem2_ssd1_v2_x8 --cost-limit 4 app-cwic --ssh -icredentials=project-GPgpg404jB0F71qxYqV6B980:file-GX9GKf84jB01gK7zPKQQ7kk4 -iimage="rswilson1/scramble:latest" -iproject_mount_options="-readWrite" -y

# Run scramble script as follows: bash scramble.sh HG01879 ACB HG01879.mapped.ILLUMINA.bwa.ACB.low_coverage.20120522.bam
set -exo pipefail

# Change to the root directory
pwd

# Download files using aria2
cd /

# To be run prior to setting off.
#gunzip -c project/003_230901_MSc_MEI_detection/references/hg19/hs37d5.fa.gz > scratch/hs37d5.fa
#samtools faidx scratch/hs37d5.fa -o scratch/hs37d5.fa.fai

# Download files using aria2
cd scratch/
# /usr/bin/time --verbose aria2c -x 16 -s 16 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/"$1"/alignment/"$1".mapped.ILLUMINA.bwa."$2".low_coverage.20120522.bam
# /usr/bin/time --verbose aria2c -x 16 -s 16 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/"$1"/alignment/"$1".mapped.ILLUMINA.bwa."$2".low_coverage.20120522.bam.bai
ID=$1
URL=$2
filename=$3

aria2c -x 16 -s 16 "${URL}" # -o "${ID}.bam"
aria2c -x 16 -s 16 "${URL}.bai" # -o "${ID}.bam.bai"

# Run Mobster
cd /
mkdir "scratch/${ID}/"
cd /home/cwic/mobster/target/
pwd
ls
java -Xmx8G -jar MobileInsertions-0.2.4.1.jar -properties Mobster.properties -in /scratch/"${filename}" -out "/scratch/${ID}/${ID}_mobster"

cd /home/cwic/mobster/resources/MobsterVCF/
java -Xmx8G -jar MobsterVCF-0.0.1-SNAPSHOT.jar -file "/scratch/${ID}/${ID}_mobster_predictions.txt" -out "/scratch/${ID}/${ID}_mobster_predictions.vcf"

#make directory for output files
mkdir -p /project/003_230901_MSc_MEI_detection/benchmarking_output/"${ID}"/mobster

# Move output files to the output directory
mv /scratch/"${ID}"/*.vcf /project/003_230901_MSc_MEI_detection/benchmarking_output/"${ID}"/mobster/
mv /scratch/"${ID}"/*.bam /project/003_230901_MSc_MEI_detection/benchmarking_output/"${ID}"/mobster/
mv /scratch/"${ID}"/*.bai /project/003_230901_MSc_MEI_detection/benchmarking_output/"${ID}"/mobster/
mv /scratch/"${ID}"/*stat* /project/003_230901_MSc_MEI_detection/benchmarking_output/"${ID}"/mobster/
mv /scratch/"${ID}"/*.pdf /project/003_230901_MSc_MEI_detection/benchmarking_output/"${ID}"/mobster/
mv /scratch/"${ID}"/*.dat /project/003_230901_MSc_MEI_detection/benchmarking_output/"${ID}"/mobster/
mv /scratch/"${ID}"/*.fq /project/003_230901_MSc_MEI_detection/benchmarking_output/"${ID}"/mobster/

rm /scratch/*.bam
rm -r /scratch/"${ID}"/

# save the output files

dx-save-project
