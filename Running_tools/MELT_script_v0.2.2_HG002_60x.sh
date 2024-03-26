#!/bin/bash

# Set up for CWIC
#dx run --instance-type mem2_ssd1_v2_x8 --cost-limit 4 app-cwic --ssh -icredentials=project-GPgpg404jB0F71qxYqV6B980:file-GX9GKf84jB01gK7zPKQQ7kk4 -iimage="rswilson1/melt:latest" -iproject_mount_options="-readWrite" -y

# Run melt script as follows: bash run_melt.sh HG01879 ACB
set -exo pipefail

# Change to the root directory
pwd

# Download files using aria2
cd /

# Download files using aria2
cd scratch/

# Run MELT
cd /

mkdir -p /scratch/"$1"/

/usr/bin/time --verbose java -jar ./opt/MELTv2.2.2/MELT.jar Single -a -b hs37d5/NC_007605 \
  -h scratch/hs37d5.fa \
  -bamfile project/003_230901_MSc_MEI_detection/HG002/hs37d5_aligned_240320/HG002.hs37d5.60x.1.bam \
  -n ./opt/MELTv2.2.2/add_bed_files/1KGP_Hg19/hg19.genes.bed \
  -t project/003_230901_MSc_MEI_detection/test_MELT_docker/mei_list_complete_dockerenv.txt \
  -w scratch/"$1"/


#make directory for output files
mkdir -p /project/003_230901_MSc_MEI_detection/benchmarking_output/"$1"/MELT/

# Move output files to the output directory
mv scratch/"$1"/*.vcf project/003_230901_MSc_MEI_detection/benchmarking_output/"$1"/MELT/
mv scratch/"$1"/*.list project/003_230901_MSc_MEI_detection/benchmarking_output/"$1"/MELT/
mv scratch/"$1"/*.tsv project/003_230901_MSc_MEI_detection/benchmarking_output/"$1"/MELT/
mv scratch/"$1"/*final* project/003_230901_MSc_MEI_detection/benchmarking_output/"$1"/MELT/

# save the output files

dx-save-project
