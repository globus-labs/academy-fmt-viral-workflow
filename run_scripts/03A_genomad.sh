#!/bin/bash
#SBATCH --output=./logs/03A_genomad/Job-%a.out
#SBATCH --account=bhurwitz
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5000
#SBATCH --array=0-52

pwd; hostname; date

source ./config.sh
names=($(cat $XFILE_DIR/$XFILE))
SAMPLE_ID=${names[${SLURM_ARRAY_TASK_ID}]}

SPADES=${SPADES_DIR}/${SAMPLE_ID}/contigs.fasta.gz
UNZIPPED_SPADES=${SPADES_DIR}/${SAMPLE_ID}/contigs.fasta

#load environment
CONDA="/groups/bhurwitz/miniconda3"
source $CONDA/etc/profile.d/conda.sh
conda activate genomad_env

genomad end-to-end --cleanup  ${UNZIPPED_SPADES} ${OUT_GENOMAD}/${names[${SLURM_ARRAY_TASK_ID}]} ${GENOMAD_DB}

