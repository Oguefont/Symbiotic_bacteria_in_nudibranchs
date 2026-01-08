#!/bin/bash
#SBATCH --job-name=gtdbtk
#SBATCH --qos=short
#SBATCH --array=1-31%6
#SBATCH --cpus-per-task=32
#SBATCH --mem=150gb
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/gtdbtk_%A_%a.out
#SBATCH -e logs/gtdbtk_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=oleanna.guerra@uv.es

module load anaconda
mamba activate gtdbtk

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

# --- gtdbtk --- #

echo "Classifying bins with gtdbtk classify_wf: ${GTDBTK_DIR}"
gtdbtk classify_wf --genome_dir $BIN_DIR --out_dir $GTDBTK_DIR --cpus 32 --extension fa --skip_ani_screen

