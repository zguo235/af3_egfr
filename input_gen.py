import json

# Define the protein sequences
proteins = [
    {
        "id": ["A"],
        "sequence": (
            "SVEMTQSPSSFSVSLGDRVTITCKASEDIYNRLAWYQQKPGNAPRLLISGATSLETEVPSRFSGSGSGKDYTLSITSLQTEDVATYYCQQYWSTWTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC"
        )
    },
    {
        "id": ["B"],
        "sequence": (
            "AVKLQESGPGILKPSQTLSLTCSFSGFSLTTYGMGVGWIRQSSGKGLEWLAHIWWDDDKYYNPSLKSRLTISKDTSRNQVFLKITSVATADTATYYCARRAPFYGNHAMDYWGQGTTVTVSSAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSPRPSETVTCNVAHPASSTKVDKKI"
        )
    },
    {
        "id": ["C"],
        "sequence": (
            "EMGTADLGPSSVPTPTNVTIESYNMNPIVYWEYQIMPQVPVFTVEVKNYGVKNSEWIDACINISHHYCNISDHVGDPSNSLWVRVKARVGQKESAYAKSEEFAVSRDG"
        )
    }
]

# Generate a list of 100 unique random seeds
model_seeds = list(range(1, 101))

# Construct the input dictionary
input_data = {
    "name": "Protein_Complex_1JRH",
    "sequences": [{"protein": protein} for protein in proteins],
    "modelSeeds": model_seeds,
    "dialect": "alphafold3",
    "version": 1
}

# Define the output file path
output_file = "input.json"

# Write the JSON data to a file
with open(output_file, "w") as f:
    json.dump(input_data, f, indent=2)

print(f"JSON input file '{output_file}' has been generated successfully.")
