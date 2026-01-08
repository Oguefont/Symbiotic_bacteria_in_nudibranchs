#!/bin/bash
#SBATCH --job-name=gtdbtk_aegi
#SBATCH --qos=short
#SBATCH --cpus-per-task=64
#SBATCH --mem=150gb
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/gtdbtk_%j.out
#SBATCH -e logs/gtdbtk_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=oleanna.guerra@uv.es


module load anaconda
mamba activate gtdbtk

gtdbtk classify_wf --genome_dir 30_junio/C --out_dir 30_junio/C/gtdbtk --cpus 64 --extension fa --skip_ani_screen

