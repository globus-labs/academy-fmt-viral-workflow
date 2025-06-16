#!/bin/bash
#SBATCH --output=./logs/01B_checkv_virsorter2/Job-%a.out
#SBATCH --account=bhurwitz
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --array=0-52

pwd; hostname; date

source ./config.sh
names=($(cat $XFILE_DIR/$XFILE))
SAMPLE_ID=${names[${SLURM_ARRAY_TASK_ID}]}

VIRSORT=${OUT_VIRSORT}/${SAMPLE_ID}/final-viral-combined.fa

#load environment
CONDA="/groups/bhurwitz/miniconda3"
source $CONDA/etc/profile.d/conda.sh
conda activate checkv_env

checkv end_to_end ${VIRSORT} ${OUT_CHECKV_VIRSORT}/${names[${SLURM_ARRAY_TASK_ID}]} -t 4
