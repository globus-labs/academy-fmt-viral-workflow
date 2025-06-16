#!/bin/bash
#SBATCH --output=./logs/02A_dvf/Job-%a.out
#SBATCH --account=bhurwitz
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
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
conda activate dvf

cd ${DVF_DB}

python dvf.py -i ${UNZIPPED_SPADES} -o ${OUT_DVF}/${names[${SLURM_ARRAY_TASK_ID}]} -l 1500

cd ${WORK_DIR}

cd ${SCRIPT_DIR}

./id_from_fasta.py -c ${UNZIPPED_SPADES} -d ${OUT_DVF}/${names[${SLURM_ARRAY_TASK_ID}]}/contigs.fasta_gt1500bp_dvfpred.txt -o ${OUT_DVF}/${names[${SLURM_ARRAY_TASK_ID}]}/dvf.fasta

cd ${WORK_DIR}
