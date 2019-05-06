#!/bin/bash

# run_blast2rma.sh
# 11.30.17
# This script builds a series of file paths and sends a command  to /opt/megan/tools/blast2rma
#    which is a program that computes MEGAN .rma files from BLAST files

# Set the parent path so we can use absolute paths later
parent_path=$(cd "$(dirname "$0")/.."; pwd)

# Set the paths to the directories using the parent_path variable
BLASTX_DIR=${parent_path}/data/MEGAN_input/BLASTX/
RMA_DIR=${parent_path}/analysis/results/RMA/
FASTA_DIR=${parent_path}/data/MEGAN_input/FASTA/

# Grab each BLASTX file, use the basename to generate the pattern
#    to find each sample's associate fasta file and to generate the
#    name of each sample's .rma file from MEGAN.

blast_files=($(find ${BLASTX_DIR} -name "*out*"))

for i in ${blast_files[@]}
do
	sample=`echo $i | awk -F "." '{print $1}'` # sample name with path but w/o the .out
	sample_basename=$(basename "${sample}") # sample name w/o path
	fasta=($(find ${FASTA_DIR} -name "${sample_basename}*.fa")) # includes path
	rma=${sample_basename}.rma # doesn't inclue a path

	/opt/megan/tools/blast2rma -i $i -f BlastText -bm BlastX -o $RMA_DIR$rma -r $fasta -ms 40 -me 0.001 -mpi 0 -top 10 -supp 0 -mrc 0 -v
done
