#!/bin/bash
#SBATCH --job-name=bwasw_BergSte
#SBATCH --qos=short
#SBATCH --array=1-6
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/alignment_%A_%a.out
#SBATCH --error=logs/alignment_%A_%a.err

module load anaconda
conda activate nudis_cont

# Separate tasks by SPECIES
SPECIES=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sps_multiple_fq.txt)

# Directories to be used
BASE_DIR="/storage/gcim/oleanna/nudi_genomes/Alignments/${SPECIES}"
INPUT_DIR="${BASE_DIR}/input"
ALIGN_DIR="${BASE_DIR}/alignments"
COV_DIR="${BASE_DIR}/coverage"
LOG_DIR="${BASE_DIR}/logs"

# universal files
ACC_FILE="${BASE_DIR}/run_acc.txt"
FASTA="${INPUT_DIR}/${SPECIES}.fna"
MERGED_BAM="${ALIGN_DIR}/${SPECIES}.merged.bam"
COV_OUT="${COV_DIR}/${SPECIES}.coverage.txt"

# Separate by accession number (one per line in $ACC_FILE)
cat $ACC_FILE | while read ACCESSION
do
        # files to be used
        FQ1="${INPUT_DIR}/${ACCESSION}_1.fastq.gz"
        FQ2="${INPUT_DIR}/${ACCESSION}_2.fastq.gz"
        FQ="${INPUT_DIR}/${ACCESSION}*.fastq.gz"
	    SAM="${ALIGN_DIR}/${ACCESSION}.sam"
        BAM="${ALIGN_DIR}/${ACCESSION}.bam"
        SORTED_SAM="${ALIGN_DIR}/${ACCESSION}.sorted.sam"
        SORTED_BAM="${ALIGN_DIR}/${ACCESSION}.sorted.bam"

        # --- Index genome --- #
        if [[ ! -f "${FASTA}.bwt" ]]; then
                echo "Indexing genome: ${SPECIES}"
                bwa index "$FASTA"
        else
                echo "Already indexed: ${SPECIES}"
        fi

        # --- Align reads --- #

        # Check if reads are paired and create sam
        if [[ ! -f "$FQ2" ]]; then
		echo "Aligning unpaired reads: ${ACCESSION}"
		bwa bwasw -t 8 "$FASTA" $FQ > "$SAM"
        else
                echo "Aligning paired reads: ${ACCESSION}"
                bwa bwasw -t 8 "$FASTA" $FQ1 $FQ2 > "$SAM"
        fi

        # Get sorted bam file
        if [[ ! -f "$SORTED_BAM" ]]; then
                echo "Creating and sorting bam: ${SPECIES}"
                samtools view -@ 8 -bS "$SAM" > "$BAM"
                samtools sort -@ 8 "$BAM" -o "$SORTED_BAM"
        else
                echo "${SORTED_BAM} already exists"
        fi
done

# --- Merge Bam files --- #
if [[ ! -f "$MERGED_BAM" ]]; then
        echo "Merging all bam files: ${MERGED_BAM}"
        samtools merge -@ 8 -o $MERGED_BAM $ALIGN_DIR/*.sorted.bam
        samtools index $MERGED_BAM
else
        echo "${MERGED_BAM} already exists"
fi

# --- Get coverage --- #
if [[ ! -f "$COV_OUT" ]]; then 
	echo "Calculating coverage: ${SPECIES}"
	# coverage with ascii art
	samtools coverage "$MERGED_BAM" > "$COV_OUT"
	
	# only coverage
	samtools depth -a "$MERGED_BAM" > "${COV_DIR}/${SPECIES}.only_cover.txt"
else
	echo "$COV_OUT" already exists
fi

# --- Delete intermediate files --- #
echo "intermediate files removed: ${ALIGN_DIR}/*.sam"
rm "${ALIGN_DIR}/*.sam"