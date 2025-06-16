#!/bin/bash
#SBATCH --output=./logs/unzip/Job-%a.out
#SBATCH --account=bhurwitz
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000
#SBATCH --array=0-52

pwd; hostname; date

source ./config.sh
names=($(cat $XFILE_DIR/$XFILE))
SAMPLE_ID=${names[${SLURM_ARRAY_TASK_ID}]}

SPADES=${SPADES_DIR}/${SAMPLE_ID}/contigs.fasta.gz
UNZIPPED_SPADES=${SPADES_DIR}/${SAMPLE_ID}/contigs.fasta

gzip -d ${SPADES}
