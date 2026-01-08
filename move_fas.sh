#!/bin/bash

# Base directories
FASTA_DIR="./fastas"
BASE_DIR="./Nudibranch_Files"

# Loop through each species in base directory
for fasta in "$BASE_DIR"/*; do
    # Remove trailing slash and extract the species name
    # Original naming pattern: "Gen_sps_sc" or "Gen_sps_chr"
    species_name=$(basename "$fasta" .fna)
    echo "Processing species: $species_name"

    # Define the new structured directory path
    base_path="${BASE_DIR}/${species_name}"
    input_dir="${base_path}/input"
    alignments_dir="${base_path}/alignments"
    bins_dir="${base_path}/bins"
    coverage_dir="${base_path}/coverage"
    logs_dir="${base_path}/logs"

    # Try creating the directories
    mkdir "$base_path"
    for dir in "$input_dir" "$alignments_dir" "$bins_dir" "$coverage_dir" "$logs_dir"; do
        if mkdir -p "$dir"; then
            echo "    Created directory: $dir"
        else
            echo "    ERROR: Failed to create $dir"
        fi
    done

    # Find matching fasta file (either _sc.fna or _chr.fna)
    fasta_file=$(find "$FASTA_DIR" -type f -name "${species_name}*.fna" | head -n 1)

    if [[ -f "$fasta_file" ]]; then
        echo "    Found FASTA: $fasta_file"
        # Move into input/
        cp "$fasta_file" "$input_dir/"
    else
        echo "    WARNING: No matching FASTA found for $species_name"
    fi
done

done
