#!/bin/bash
#SBATCH --job-name=bwa_mem_alignment
#SBATCH --qos=short
#SBATCH --array=1-26%6
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/alignment_%A_%a.out
#SBATCH --error=logs/alignment_%A_%a.err

module load anaconda
conda activate nudis_cont

# Separate tasks by species
SPECIES=$(sed -n "${SLURM_ARRAY_TASK_ID}p" short_paired.txt)

# Directories to be used
BASE_DIR="Alignments/${SPECIES}"
INPUT_DIR="${BASE_DIR}/input"
ALIGN_DIR="${BASE_DIR}/alignments"
COV_DIR="${BASE_DIR}/coverage"
LOG_DIR="${BASE_DIR}/logs"

# files to be used
FASTA="${INPUT_DIR}/${SPECIES}.fna"
FQ1="${INPUT_DIR}/*_1.fastq.gz"
FQ2="${INPUT_DIR}/*_2.fastq.gz"
SAM="${ALIGN_DIR}/${SPECIES}.sam"
BAM="${ALIGN_DIR}/${SPECIES}.bam"
SORTED_SAM="${ALIGN_DIR}/${SPECIES}.sorted.sam"
SORTED_BAM="${ALIGN_DIR}/${SPECIES}.sorted.bam"
COV_OUT="${COV_DIR}/${SPECIES}.coverage.txt"

# --- Index genome --- #
if [[ ! -f "${FASTA}.bwt" ]]; then
	echo "Indexing genome: ${SPECIES}"
	bwa index "$FASTA"
else
	echo "Already indexed: ${SPECIES}"
fi

# --- Align reads --- #

# Get sam file
if [[ ! -f "$SAM" ]]; then
	echo "Aligning reads: ${SPECIES}"
	bwa mem -t 8 "$FASTA" $FQ1 $FQ2 > "$SAM"
else
	echo "${SAM} already exists"
fi

# Get sorted bam file
if [[ ! -f "$SORTED_BAM" ]]; then
	echo "Creating and sorting bam: ${SPECIES}"
	samtools view -@ 8 -bS "$SAM" > "$BAM"
	samtools sort -@ 8 "$BAM" -o "$SORTED_BAM"
	samtools index "$SORTED_BAM"
else
	echo "${SORTED_BAM} already exists"
fi

# --- Get coverage --- #
if [[ ! -f "$COV_OUT" ]]; then 
	echo "Calculating coverage: ${SPECIES}"
	# coverage by contig
	samtools coverage "$SORTED_BAM" > "$COV_OUT"
	
	# coverage/depth of each
	samtools depth -a "$SORTED_BAM" > "${COV_DIR}/${SPECIES}.only_cover.txt"
else
	echo "$COV_OUT" already exists
fi