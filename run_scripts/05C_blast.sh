#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --partition=standard
#SBATCH --account=bhurwitz
#SBATCH --output=./logs/05C_blast.out
#SBATCH --error=./logs/05C_blast.err
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=5G

#
# This script runs blast jobs for each of the split fasta files
#

pwd; hostname; date
source $WORK_DIR/config.sh

# get the split file, that we are currently blasting
names=($(cat ${SPLIT_FILES_LIST}))
SPLIT_FILE=${names[${SLURM_ARRAY_TASK_ID}]}

# set up the results, stderr and stdout directories for this script
PROG="05C_blast"
RESULTS_DIR="$WORK_DIR/results/$PROG"

# create dir if it doesn't exist
create_dir "$RESULTS_DIR"

# we want to blast each array value for each split file against each database
i=0
while read DB; do
    let i++

    #NAME=`basename $SPLIT_FILE` 
    printf "%5d: %s\n" $i "$DB"
    RESULTS_BY_DB="$RESULTS_DIR/$DB/$FILE_NAME"
    create_dir "$RESULTS_BY_DB"
    BLAST_OUT="$RESULTS_BY_DB/$SPLIT_FILE"
    BLAST_DB="$DB_DIR/$DB"

    # run blast against each split file and database
    apptainer run ${BLAST} $BLAST_TYPE -num_threads 48 -db $BLAST_DB -query $SPLIT_DIR/$SPLIT_FILE -out $BLAST_OUT -evalue $EVAL -outfmt $OUT_FMT -max_target_seqs $MAX_TARGET_SEQS
done < "$DB_DIR/db-list"

echo "Finished `date`"
