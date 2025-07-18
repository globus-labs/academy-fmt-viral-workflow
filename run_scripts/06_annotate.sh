#!/bin/bash

pwd; hostname; date
source ./config_py.sh

BLAST_HITS="/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline/results_testing/05D_mergeblast/AVrC_allrepresentatives.fasta/clusterRes_rep_seq.fasta.txt"
ANNOTATIONS="/xdisk/bhurwitz/databases/AVrC/database_csv/"
OUTPUT="/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline/results_testing/06_annotate"

# Debug: Print all files in the annotations directory

echo "Files in annotation directory:"
echo ${ANNOTATIONS}*

mkdir -p /xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/06_annotate

# Loop over all files in the annotation directory

for file in $(ls ${ANNOTATIONS}); do

    ${SCRIPT_DIR}/solution1_manual.py -b ${BLAST_HITS} -a ${ANNOTATIONS}"$file" -o ${OUTPUT}/annotated_${file} -p ${PCTID} -l ${LENGTH}

done


