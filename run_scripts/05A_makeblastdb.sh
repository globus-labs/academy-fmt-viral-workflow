#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --partition=standard
#SBATCH --account=bhurwitz
#SBATCH --output=./logs/05A_makeblastdb.out
#SBATCH --error=./logs/05A_makeblastdb.err
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G

#
# This script will create blast databases with all *fa files in the DB_DIR 
#

pwd; hostname; date

# get the configurations
source $WORK_DIR/config.sh

cd "$DB_DIR"

export DB_LIST="db-list"

find . -type f -name \*.fasta | sed "s/^\.\///" > $DB_LIST

if [[ ! -e "$DB_LIST" ]]; then
    echo Cannot find database list \"$DB_LIST\"
    exit 1
fi

i=0
while read DB; do
    let i++

    DB_NAME=`basename $FILE`

    printf "%5d: %s\n" $i "$DB_NAME"

    #
    # create blast database for fasta file 
    #
    apptainer run ${BLAST} makeblastdb -title ${DB_NAME} -out ${DB_NAME} -in ${DB} -dbtype nucl -max_file_sz ${MAX_DB_SIZE} 
done < "$DB_LIST"

echo Finished `date`
