#!/bin/bash
## Version: 1.2

module load anaconda
conda activate nudis_cont

# Get Species from file
FAS_NAMES_TXT="/storage/gcim/oleanna/nudi_genomes/fasta_names.txt"
while read -r line; do
    SPECIES=$line

    # Directories to be used
    BASE_DIR="/storage/gcim/oleanna"
    SPS_DIR="${BASE_DIR}/nudi_genomes/Alignments/${SPECIES}"
    INPUT_DIR="${SPS_DIR}/input"
	BIN_DIR="${SPS_DIR}/bins"
    FASTA="${INPUT_DIR}/${SPECIES}.fna"
    SP=$(head -1 $FASTA | awk '{print $2"_"$3}')
    COV_DIR="${SPS_DIR}/coverage"
    BARR_DIR="${BASE_DIR}/Barrnap_dir/Gff_files"
    K2_DIR="${BASE_DIR}/Kraken2"
    CONTIG_DIR="${BASE_DIR}/nudi_genomes/contig_tables"

    # files to be used (FASTA described previously)
    K2_FILE="${K2_DIR}/Std_out/${SP}_kraken2.out"
    BARR_FILE="${BARR_DIR}/${SP}.gff"
    COV_FILE="${COV_DIR}/${SPECIES}.coverage.txt"
    TSV_IN="${CONTIG_DIR}/${SPECIES}.fna.tsv"
    TSV_OUT="${CONTIG_DIR}/${SPECIES}_merged_results.tsv"
	BIN_FILE="${BIN_DIR}/contig_bin.tsv"

    # --- Barrnap results --- #
    if [[ -f "${BARR_FILE}" ]]; then
        echo "Adding information from ${BARR_FILE}"
        python3 results_merger.py -i $TSV_IN -r $BARR_FILE -p b -o $TSV_OUT -d $BASE_DIR
    else
        echo "File not found: ${BARR_FILE}"
    fi

    # --- Kraken2 results --- #
    if [[ -f "${K2_FILE}" ]]; then
        echo "Adding information from ${K2_FILE}"
        python3 results_merger.py -i $TSV_IN -r $K2_FILE -p k2 -o $TSV_OUT -d $BASE_DIR
    else
        echo "File not found: ${K2_FILE}"
    fi

    # --- Coverage results --- #
    if [[ -f "${COV_FILE}" ]]; then
        echo "Adding information from ${COV_FILE}"
        python3 results_merger.py -i $TSV_IN -r $COV_FILE -p cov -o $TSV_OUT -d $BASE_DIR
    else
        echo "File not found: ${COV_FILE}"
    fi

	# --- Metabat results --- #
	if [[ -f "${BIN_FILE}" ]]; then
		echo "Ading information from ${BIN_FILE}"
		python3 results_merger.py -i $TSV_IN -r $BIN_FILE -p m -o $TSV_OUT -d $BASE_DIR
	else
		echo "File not found: ${BIN_FILE}"
	fi

done < "$FAS_NAMES_TXT"