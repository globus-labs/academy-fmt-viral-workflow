#!/bin/bash
#SBATCH --output=./logs/03B_checkv_genomad/Job-%a.out
#SBATCH --account=bhurwitz
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --array=0-2

pwd; hostname; date
source ./config_py.sh

names=($(cat $XFILE_DIR/$XFILE))
SAMPLE_ID=${names[${SLURM_ARRAY_TASK_ID}]}

#SPADES=${SPADES_DIR}/${SAMPLE_ID}/contigs.fasta.gz
UNZIPPED_SPADES=${SPADES_DIR}/${SAMPLE_ID}/contigs.fasta

GENOMAD=${OUT_GENOMAD}/${SAMPLE_ID}/contigs_summary/contigs_virus.fna
PARSE_INPUT=${OUT_CHECKV_GENOMAD}/${SAMPLE_ID}/contamination.tsv

#load environment
CONDA="/xdisk/bhurwitz/miniconda3"
source $CONDA/etc/profile.d/conda.sh
conda activate checkv_env

checkv end_to_end ${GENOMAD} ${OUT_CHECKV_GENOMAD}/${SAMPLE_ID} -t 4
conda deactivate

cd ${OUT_CHECKV_GENOMAD}/${SAMPLE_ID}/

conda activate r_env
Rscript ${CHECKV_PARSER} -i ${PARSE_INPUT} -l ${PARSE_LENGTH} -o selection2_viral.csv
conda deactivate

conda activate seqtk_env
seqtk subseq ${UNZIPPED_SPADES} selection2_viral.csv > subset_spades.fasta
cd ${WORK_DIR}
