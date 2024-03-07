#!/bin/bash

# Specify the TSV file containing file IDs, names, and folders
TSV_FILE="melt_file_ids.tsv"

# Specify the project ID
PROJECT_ID="project-GQ8p00041kgj8VXvGVQ64V21"

# Loop over the TSV file
while IFS=$'\t' read -r file_id name folder; do
    # Remove the ".vcf.gz" extension from the file name
    name_no_extension="${name%.vcf.gz}"

    # Construct the output file names without the extension
    test_vcf_filtered_ASSESS_1="${name_no_extension}_filtered_ASSESS_eqgt_1.vcf"
    test_vcf_filtered_ASSESS_2="${name_no_extension}_filtered_ASSESS_eqgt_2.vcf"
    test_vcf_filtered_ASSESS_3="${name_no_extension}_filtered_ASSESS_eqgt_3.vcf"
    test_vcf_filtered_ASSESS_4="${name_no_extension}_filtered_ASSESS_eqgt_4.vcf"
    test_vcf_filtered_ASSESS_5="${name_no_extension}_filtered_ASSESS_eqgt_5.vcf"

    # Construct the BCFtools command with the appropriate file ID
    bcftools_cmd_ASSESS_1="tabix -p vcf ${name}; bcftools view -i 'ASSESS>=1' ${name} -o ${test_vcf_filtered_ASSESS_1}; bgzip ${test_vcf_filtered_ASSESS_1}"
    bcftools_cmd_ASSESS_2="tabix -p vcf ${name}; bcftools view -i 'ASSESS>=2' ${name} -o ${test_vcf_filtered_ASSESS_2}; bgzip ${test_vcf_filtered_ASSESS_2}"
    bcftools_cmd_ASSESS_3="tabix -p vcf ${name}; bcftools view -i 'ASSESS>=3' ${name} -o ${test_vcf_filtered_ASSESS_3}; bgzip ${test_vcf_filtered_ASSESS_3}"
    bcftools_cmd_ASSESS_4="tabix -p vcf ${name}; bcftools view -i 'ASSESS>=4' ${name} -o ${test_vcf_filtered_ASSESS_4}; bgzip ${test_vcf_filtered_ASSESS_4}"
    bcftools_cmd_ASSESS_5="tabix -p vcf ${name}; bcftools view -i 'ASSESS>=5' ${name} -o ${test_vcf_filtered_ASSESS_5}; bgzip ${test_vcf_filtered_ASSESS_5}"
    #  && FILTER=\"PASS\"

    # Run the swiss_army_knife app with the constructed command
    dx run app-swiss-army-knife \
    -i in="${PROJECT_ID}":"${file_id}" \
    -i cmd="$bcftools_cmd_ASSESS_1" \
    --destination="${folder}" -y

    sleep 1

    dx run app-swiss-army-knife \
    -i in="${PROJECT_ID}":"${file_id}" \
    -i cmd="$bcftools_cmd_ASSESS_2" \
    --destination="${folder}" -y

    sleep 1

    dx run app-swiss-army-knife \
    -i in="${PROJECT_ID}":"${file_id}" \
    -i cmd="$bcftools_cmd_ASSESS_3" \
    --destination="${folder}" -y

    sleep 1

    dx run app-swiss-army-knife \
    -i in="${PROJECT_ID}":"${file_id}" \
    -i cmd="$bcftools_cmd_ASSESS_4" \
    --destination="${folder}" -y

    sleep 1

    dx run app-swiss-army-knife \
    -i in="${PROJECT_ID}":"${file_id}" \
    -i cmd="$bcftools_cmd_ASSESS_5" \
    --destination="${folder}" -y

    sleep 5

    # Print a message indicating the completion of processing for the current file
    echo "Processed file ${name}"
done < "$TSV_FILE"
