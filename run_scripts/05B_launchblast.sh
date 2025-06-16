#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=standard
#SBATCH --account=bhurwitz
#SBATCH --output=./logs/05B_launchblast.out
#SBATCH --error=./logs/05B_launchblast.err
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5G

#
# This script splits the input files, blasts them against the databases, and collates results
# Note that this script runs additional jobs on each of the files as a job array for the split files
#

pwd; hostname; date
source ./config.sh

# set up the results, stderr and stdout directories for this script
PROG="05B_launchblast"
RESULTS_DIR="$WORK_DIR/results/$PROG"

# initialize directories, this will remove prior runs and create new directories
init_dir "$RESULTS_DIR"

cd "$FASTA_DIR"

export FILES_LIST="fasta-files"

find . -type f -name \*.fasta | sed "s/^\.\///" > $FILES_LIST

NUM_FILES=$(lc $FILES_LIST)

echo Found \"$NUM_FILES\" files in \"$FASTA_DIR\"

if [ $NUM_FILES -gt 0 ]; then
    i=0
    while read FILE; do
        let i++
    
        export FILE_NAME=`basename $FILE`
    
        printf "%5d: %s\n" $i "$FILE_NAME"

        OUT_DIR="$RESULTS_DIR/$FILE_NAME"
        export SPLIT_DIR="$OUT_DIR/fa_split"
        init_dir "$OUT_DIR" "$SPLIT_DIR" 

        # first we need to split up the fasta files to run quickly
        apptainer run ${FASPLIT} faSplit about "$FILE" "$FA_SPLIT_FILE_SIZE" "$SPLIT_DIR/"

        # now launch the blast for the split files against all of the blast databases
        cd "$SPLIT_DIR"
        export SPLIT_FILES_LIST="$SPLIT_DIR/split-fasta-files"
        find . -type f -name \*.fa | sed "s/^\.\///" > $SPLIT_FILES_LIST

        NUM_SPLIT_FILES=$(lc $SPLIT_FILES_LIST)
        ARRAY_NUM=$((NUM_SPLIT_FILES - 1))
        echo Found \"$NUM_SPLIT_FILES\" files in \"$SPLIT_DIR\"

        cd $WORK_DIR

        # run all blast jobs for each split file
        job1=$(sbatch --array=0-${ARRAY_NUM} ${SCRIPT_DIR}/05C_blast.sh)
        jid1=$(echo $job1 | sed 's/^Submitted batch job //')
        echo $jid1 $FILE_NAME $SPLIT_DIR

        # collate the results for each split file
        job2=$(sbatch --dependency=afterok:$jid1 ${SCRIPT_DIR}/05D_mergeblast.sh)
        jid2=$(echo $job2 | sed 's/^Submitted batch job //')
        echo $jid2 $FILE_NAME $SPLIT_DIR

        cd "$FASTA_DIR"
    done < "$FILES_LIST"
else
    echo No input fasta files.
fi

echo "Finished `date`"
