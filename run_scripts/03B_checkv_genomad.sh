#!/bin/bash
#SBATCH --output=./logs/03B_checkv_genomad/Job-%a.out
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

#SPADES=${SPADES_DIR}/${SAMPLE_ID}/contigs.fasta.gz
UNZIPPED_SPADES=${SPADES_DIR}/${SAMPLE_ID}/contigs.fasta

GENOMAD=${OUT_GENOMAD}/${SAMPLE_ID}/contigs_summary/contigs_virus.fna
PARSE_INPUT=${OUT_CHECKV_GENOMAD}/${SAMPLE_ID}/contamination.tsv

#load environment
CONDA="/groups/bhurwitz/miniconda3"
source $CONDA/etc/profile.d/conda.sh
conda activate checkv_env

checkv end_to_end ${GENOMAD} ${OUT_CHECKV_GENOMAD}/${SAMPLE_ID} -t 4
conda deactivate
conda activate seqtk_env

cd ${OUT_CHECKV_GENOMAD}/${SAMPLE_ID}/
${RSCRIPT_DIR} ${CHECKV_PARSER} -i ${PARSE_INPUT} -l ${PARSE_LENGTH}
seqtk subseq ${UNZIPPED_SPADES} selection2_viral.csv > subset_spades.fasta
cd ${WORK_DIR}
