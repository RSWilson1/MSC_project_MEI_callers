#!/bin/bash

# Set up for CWIC
#dx run --instance-type mem2_ssd1_v2_x8 --cost-limit 4 app-cwic --ssh -icredentials=project-GPgpg404jB0F71qxYqV6B980:file-GX9GKf84jB01gK7zPKQQ7kk4 -iimage="rswilson1/scramble:latest" -iproject_mount_options="-readWrite" -y
#gunzip -c project/003_MSc_MEI_detection/references/hg19/hs37d5.fa.gz > scratch/hs37d5.fa
#samtools faidx scratch/hs37d5.fa -o scratch/hs37d5.fa.fai
#install aria2 using apt-get (make sure to apt update first)

# Run scramble script as follows: bash scramble.sh HG01879 ACB HG01879.mapped.ILLUMINA.bwa.ACB.low_coverage.20120522.bam
set -exo pipefail

# Change to the root directory
pwd

# Download files using aria2
cd /

# Download files using aria2
cd scratch/

ID=$1
URL=$2
filename=$3

aria2c -x 16 -s 16 "${URL}" # -o "${ID}.bam"
aria2c -x 16 -s 16 "${URL}.bai" # -o "${ID}.bam.bai"

# Run SCRAMBLE
cd /
mkdir -p "scratch/${ID}/"

cluster_identifier \
    "scratch/${filename}" > "scratch/${ID}/${ID}.clusters.txt"

Rscript --vanilla /app/cluster_analysis/bin/SCRAMble.R \
    --out-name "/scratch/${ID}/${ID}_scramble" \
    --cluster-file "/scratch/${ID}/${ID}.clusters.txt" \
    --install-dir /app/cluster_analysis/bin \
    --mei-refs /app/cluster_analysis/resources/MEI_consensus_seqs.fa \
    --ref /scratch/hs37d5.fa \
    --eval-dels \
    --eval-meis

cd /

#make directory for output files
mkdir -p /project/003_230901_MSc_MEI_detection/benchmarking_output/"${ID}"/scramble/

# Move output files to the output directory
mv /scratch/"${ID}"/*.vcf project/003_230901_MSc_MEI_detection/benchmarking_output/"${ID}"/scramble/
mv /scratch/"${ID}"/*.txt project/003_230901_MSc_MEI_detection/benchmarking_output/"${ID}"/scramble/

# save the output files

dx-save-project
