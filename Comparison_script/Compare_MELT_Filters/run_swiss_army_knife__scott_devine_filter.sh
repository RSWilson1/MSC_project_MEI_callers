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
    test_vcf_filtered_SD="${name_no_extension}_filtered_SD.vcf"

    # Construct the BCFtools command with the appropriate file ID
    bctools_cmd_SD="tabix -p vcf ${name}; bcftools view -i 'FILTER==\"PASS\"  | FILTER==\"lc\"' ${name} -o ${test_vcf_filtered_SD}; bgzip ${test_vcf_filtered_SD}"

    # Run the swiss_army_knife app with the constructed command
    dx run app-swiss-army-knife \
    -i in="${PROJECT_ID}":"${file_id}" \
    -i cmd="$bctools_cmd_SD" \
    --destination="${folder}" -y

    # Print a message indicating the completion of processing for the current file
    echo "Processed file ${name}"
done < "$TSV_FILE"
