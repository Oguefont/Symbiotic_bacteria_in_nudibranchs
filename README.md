**Master's Thesis:**

------------------------------------------------------------------
# Revealing uncharted symbiotic bacteria in nudibranch genome assemblies
- Author: Oleanna Guerra-Font
- Tutors: Mária Džunková and Vicente Arnau
------------------------------------------------------------------

Associated programs mentioned in the study and created by the author. A through explanation on how each program was used can be found in the body of the original work, section 7 of the introduction, titled "programs created for this study". What follows is a brief introduction of each program.

## assembly2reads.py

Divides the original assembly into 300 bp pseudoreads, each separated by 150 bp.

## barrnap_script.sbatch

Runs Barrnap v0.9 in a HPC Cluster using slurm.

## bins_creation.sh

Launches the entire metagenomic workflow as well as whokaryote on a HPC Cluster using slurm.

Metabat2 --> bins --> CheckM2 & GTDB_Tk

## bwamem_align.sh & bwasw_align.sh

Aligns the original reads to the final assembly by using the Burrows-Wheeler Aligner (BWA) software package for short and long reads, respectivelly.

## gtdb.sh

Launches GTDB-Tk on a HPC cluster with the results from Metabat2 and CheckM2.

## merger_launcher.sh

Launches the results_merger.py program for each species.

## move_fas.sh

Organises original fasta files into coherent directories and subdirectories, following:

<img width="1902" height="1275" alt="files" src="https://github.com/user-attachments/assets/e70adbb9-45ae-41e5-8c83-704567443713" />

## results_merger.py

Merges the results obtained by each of the tools used in the proposed workflow.
