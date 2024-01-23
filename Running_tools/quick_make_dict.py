# Read the download links from the text file
with open("download_links.txt", "r") as file:
    download_links = file.read().splitlines()

# Initialize an empty dictionary
download_dict = {}

# Iterate through the download links and extract the ID
for link in download_links:
    parts = link.split('/')
    id = parts[-3]

    # Check if the link ends with ".bam" to ensure it's a BAM file
    if link.endswith(".bam"):
        download_dict[id] = link

# Print the resulting dictionary
print(download_dict)
