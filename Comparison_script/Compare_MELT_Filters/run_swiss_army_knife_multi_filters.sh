#!/bin/bash

# Specify the TSV file containing file IDs, names, and folders
TSV_FILE="melt_retry_HG01112.tsv" #"melt_file_ids.tsv"

# Specify the project ID
PROJECT_ID="project-GQ8p00041kgj8VXvGVQ64V21"

# Loop over the TSV file
while IFS=$'\t' read -r file_id name folder; do
    # Remove the ".vcf.gz" extension from the file name
    name_no_extension="${name%.vcf.gz}"

    # Construct the output file names without the extension
    test_vcf_filtered_comp="${name_no_extension}_filtered_comprehensive.vcf"
    test_vcf_filtered_ASSESS_ONLY="${name_no_extension}_filtered_ASSESS_only.vcf"
    test_vcf_filtered_PASS_ONLY="${name_no_extension}_filtered_PASS_only.vcf"
    test_vcf_filtered_STRICT="${name_no_extension}_filtered_strict.vcf"

    # Construct the BCFtools command with the appropriate file ID
    bctools_cmd_comp="tabix -p vcf ${name}; bcftools view -i 'FILTER==\"PASS\" & ASSESS>3' ${name} -o ${test_vcf_filtered_comp}; bgzip ${test_vcf_filtered_comp}"
    bctools_cmd_ASSESS_only="tabix -p vcf ${name}; bcftools view -i 'ASSESS==5' ${name} -o ${test_vcf_filtered_ASSESS_ONLY}; bgzip ${test_vcf_filtered_ASSESS_ONLY}"
    bctools_cmd_PASS_only="tabix -p vcf ${name}; bcftools view -i 'FILTER==\"PASS\"' ${name} -o ${test_vcf_filtered_PASS_ONLY}; bgzip ${test_vcf_filtered_PASS_ONLY}"
    bctools_cmd_strict="tabix -p vcf ${name}; bcftools view -i 'FILTER==\"PASS\" & ASSESS>3 & LP>=1 & RP>=1' ${name} -o ${test_vcf_filtered_STRICT}; bgzip ${test_vcf_filtered_STRICT}"

    # Run the swiss_army_knife app with the constructed command
    dx run app-swiss-army-knife \
    -i in="${PROJECT_ID}":"${file_id}" \
    -i cmd="$bctools_cmd_comp" \
    --destination="${folder}" -y

    dx run app-swiss-army-knife \
    -i in="${PROJECT_ID}":"${file_id}" \
    -i cmd="$bctools_cmd_ASSESS_only" \
    --destination="${folder}" -y

    dx run app-swiss-army-knife \
    -i in="${PROJECT_ID}":"${file_id}" \
    -i cmd="$bctools_cmd_PASS_only" \
    --destination="${folder}" -y

    dx run app-swiss-army-knife \
    -i in="${PROJECT_ID}":"${file_id}" \
    -i cmd="$bctools_cmd_strict" \
    --destination="${folder}" -y

    # Print a message indicating the completion of processing for the current file
    echo "Processed file ${name}"
done < "$TSV_FILE"
