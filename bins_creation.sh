#!/bin/bash
#SBATCH --job-name=bins_creation
#SBATCH --qos=short
#SBATCH --array=1-31%6
#SBATCH --cpus-per-task=32
#SBATCH --mem=100gb
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/bins_%A_%a.out
#SBATCH -e logs/bins_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=oleanna.guerra@uv.es

module load anaconda
mamba activate nudis_cont

# Separate tasks by SPECIES
SPECIES=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sps_FEMS.txt)


# Directories to be used
BASE_DIR="/storage/gcim/oleanna/nudi_genomes/Alignments/${SPECIES}"
INPUT_DIR="${BASE_DIR}/input"
ALIGN_DIR="${BASE_DIR}/alignments"
COV_DIR="${BASE_DIR}/coverage"
LOG_DIR="${BASE_DIR}/logs"
BIN_DIR="${BASE_DIR}/bins"
GTDBTK_DIR="${BIN_DIR}/gtdbtk"
CHECKM_DIR="${BIN_DIR}/Checkm2"
WHOKARYOTE_DIR="${BASE_DIR}/whokaryote"

# Files
FASTA="${INPUT_DIR}/${SPECIES}.fna"
MERGED_BAM="${ALIGN_DIR}/${SPECIES}.merged.bam"
SORTED_BAM="${ALIGN_DIR}/${SPECIES}.sorted.bam"
ABDFILE="${COV_DIR}/${SPECIES}_abdFile_coverage.tsv"

# --- Metabat2 --- #

# Get abd file for Metabat
if [[ ! -f $MERGED_BAM ]]; then
	echo "Summarizing bam contig depths from ${SORTED_BAM}"
	jgi_summarize_bam_contig_depths --includeEdgeBases --referenceFasta $FASTA --outputDepth $ABDFILE $SORTED_BAM 
else
	echo "Summarizing bam contig depths from ${MERGED_BAM}"
	jgi_summarize_bam_contig_depths --includeEdgeBases --referenceFasta $FASTA --outputDepth $ABDFILE $MERGED_BAM
fi

# Run metabat2
echo "Running Metabat2 with minContig 1500, minCV 1.0, minCVSum 1.0, maxP 95%, minS 60, maxEdges 200 and minClsSize 150000"
metabat2 -i $FASTA -o "${BIN_DIR}/${SPECIES}" -a $ABDFILE --unbinned -t 32 -m 1500 -s 150000

# --- CheckM2 --- #

mamba deactivate
mamba activate checkm2

echo "Checking bins with checkm2 predict: ${CHECKM_DIR}"
checkm2 predict -t 32 --input $BIN_DIR/*.fa --output-directory $CHECKM_DIR --database_path /home/oguefont/CheckM2_database/uniref100.KO.1.dmnd

# --- gtdbtk --- #

mamba deactivate
mamba activate gtdbtk

echo "Classifying bins with gtdbtk classify_wf: ${GTDBTK_DIR}"
gtdbtk classify_wf --genome_dir $BIN_DIR --out_dir $GTDBTK_DIR --cpus 32 --extension fa --skip_ani_screen

# --- whokaryote --- #

mamba deactivate
mamba activate whokaryote 

whokaryote.py --contigs $FASTA --outdir $WHOKARYOTE_DIR
