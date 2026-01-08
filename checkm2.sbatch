#!/bin/bash
#SBATCH --job-name=checkm_aegi
#SBATCH --qos=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=50gb
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/checkm2_%j.out
#SBATCH -e logs/checkm2_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=oleanna.guerra@uv.es

checkm2 predict -t 16 --input *.fa --output-directory Checkm2 --database_path /home/oguefont/CheckM2_database/uniref100.KO.1.dmnd
