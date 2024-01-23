#!/bin/bash

# Set up for CWIC
#dx run --instance-type mem2_ssd1_v2_x8 --cost-limit 4 app-cwic --ssh -icredentials=project-GPgpg404jB0F71qxYqV6B980:file-GX9GKf84jB01gK7zPKQQ7kk4 -iimage="rswilson1/scramble:latest" -iproject_mount_options="-readWrite" -y

# Run scramble script as follows: bash scramble.sh HG01879 ACB HG01879.mapped.ILLUMINA.bwa.ACB.low_coverage.20120522.bam
set -exo pipefail

# Change to the root directory
pwd

# Download files using aria2
cd /

#gunzip -c project/003_MSc_MEI_detection/references/hg19/hs37d5.fa.gz > scratch/hs37d5.fa
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

# Run SCRAMBLE
cd /
mkdir "scratch/${ID}/"

cluster_identifier \
    "scratch/${ID}/${filename}" > "scratch/$ID/$ID.clusters.txt"

Rscript --vanilla /app/cluster_analysis/bin/SCRAMble.R \
    --out-name "/scratch/$ID/" \
    --cluster-file "/scratch/$ID/$ID.clusters.txt" \
    --install-dir /app/cluster_analysis/bin \
    --mei-refs /app/cluster_analysis/resources/MEI_consensus_seqs.fa \
    --ref scratch/hs37d5.fa \
    --eval-dels \
    --eval-meis

#make directory for output files
mkdir project/003_MSc_MEI_detection/benchmarking_output/"$1"

# Move output files to the output directory
mv scratch/*.vcf project/003_MSc_MEI_detection/benchmarking_output/"$1"/
mv scratch/*.list project/003_MSc_MEI_detection/benchmarking_output/"$1"/
mv scratch/*.tsv project/003_MSc_MEI_detection/benchmarking_output/"$1"/
mv scratch/*final* project/003_MSc_MEI_detection/benchmarking_output/"$1"/


# save the output files

dx-save-project
