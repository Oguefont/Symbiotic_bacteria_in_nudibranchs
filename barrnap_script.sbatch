#!/bin/bash
#SBATCH --job-name=Barrnap_all_nudis
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --time=1-00:00:00
#SBATCH --output=Barrnap_dir/Barrnap_%j.out
#SBATCH -e Barrnap_dir/Barrnap_%j.err

module load anaconda
conda activate nudis_cont

for f in nudi_genomes/fastas/*.fna
do
	if [[ ! -f Barrnap_dir2/${sp}.gff ]]; then 
       	sp=$(head -1 $f | awk '{print $2"_"$3}')	
		barrnap --outseq Barrnap_dir2/Barr_hits__${sp}.fas  ${f} > Barrnap_dir2/${sp}.gff --threads 16
	else
		echo "${f} already analized with barrnap"
	fi
done