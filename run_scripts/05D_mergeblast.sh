#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=standard
#SBATCH --account=bhurwitz
#SBATCH --output=./logs/05D_mergeblast.out
#SBATCH --error=./logs/05D_mergeblast.err
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5G

#
# This script collates the blast results for each file against each database
#

pwd; hostname; date
source $WORK_DIR/config.sh

# set up the results, stderr and stdout directories for this script
PROG="05D_mergeblast"
RESULTS_DIR="$WORK_DIR/results/$PROG"

# initialize directories, this will remove prior runs and create new directories
create_dir "$RESULTS_DIR"

i=0
while read DB; do
    let i++
    printf "%5d: %s\n" $i "$DB_NAME"
    RESULTS_BY_DB="$RESULTS_DIR/$DB"
    create_dir "$RESULTS_BY_DB"
    BLAST_RESULTS="$RESULTS_BY_DB/${FILE_NAME}.txt"
    BLAST_GFF="$RESULTS_BY_DB/${FILE_NAME}.gff"

    BLAST_OUT="$PWD/results/05C_blast/$DB/$FILE_NAME"

    cat $BLAST_OUT/* > $BLAST_RESULTS

    # convert to GFF format
    awk '{print $1"\tblast\tgene\t"$7"\t"$8"\t.\t.\t.\tID=Gene"$7";Name="$2}' $BLAST_RESULTS > $BLAST_GFF

done < "$DB_DIR/db-list"

echo "Finished `date`"
