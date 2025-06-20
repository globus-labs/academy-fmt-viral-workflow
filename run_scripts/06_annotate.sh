#pwd; hostname; date
#source ./config.sh

#BLAST_HITS=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/05D_mergeblast/AVrC_allrepresentatives.fasta/clusterRes_rep_seq.fasta.txt
#ANNOTATIONS="xdisk/bhurwitz/databases/AVrC/annotation/database_csv/"
##OUTPUT=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/06_annotate/annotated_cluster.csv
#PCTID=85
#LENGTH=1000

#echo ${ANNOTATIONS}*

#for file in ${ANNOTATIONS}*; do
 #  ${SCRIPT_DIR}/solution1_manual.py -b ${BLAST_HITS} -a "$file" -o ${OUTPUT} -p ${PCTID} -l ${LENGTH}
#   done

#!/bin/bash



pwd; hostname; date
source ./config.sh

BLAST_HITS="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/05D_mergeblast/AVrC_allrepresentatives.fasta/clusterRes_rep_seq.fasta.txt"
ANNOTATIONS="/xdisk/bhurwitz/databases/AVrC/annotation/database_csv/"
OUTPUT="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/06_annotate"

# Debug: Print all files in the annotations directory

echo "Files in annotation directory:"
echo ${ANNOTATIONS}*

mkdir -p /xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/06_annotate

# Loop over all files in the annotation directory

for file in $(ls ${ANNOTATIONS}); do

    ${SCRIPT_DIR}/solution1_manual.py -b ${BLAST_HITS} -a ${ANNOTATIONS}"$file" -o ${OUTPUT}/annotated_${file} -p ${PCTID} -l ${LENGTH}

done


