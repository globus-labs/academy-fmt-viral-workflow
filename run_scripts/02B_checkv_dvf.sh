#!/bin/bash
#SBATCH --output=./logs/02B_checkv_dvf/Job-%a.out
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

DVF=${OUT_DVF}/${SAMPLE_ID}/contigs.fasta_gt1500bp_dvfpred.txt

SPADES=${SPADES_DIR}/${SAMPLE_ID}/contigs.fasta.gz 
UNZIPPED_SPADES=${SPADES_DIR}/${SAMPLE_ID}/contigs.fasta

DVF_FASTA=${OUT_DVF}/${SAMPLE_ID}/dvf.fasta

#load environment
CONDA="/groups/bhurwitz/miniconda3"
source $CONDA/etc/profile.d/conda.sh
conda activate checkv_env

${SCRIPT_DIR}/id_from_fasta.py -c ${UNZIPPED_SPADES} -d ${DVF} -o ${DVF_FASTA}

checkv end_to_end ${DVF_FASTA} ${OUT_CHECKV_DVF}/${names[${SLURM_ARRAY_TASK_ID}]} -t 16

