#!/bin/bash
#SBATCH --output=./logs/04A_dereplicate/Job-%a.out
#SBATCH --account=bhurwitz
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5000
#SBATCH --array=0-52

pwd; hostname; date
source ./config.sh

names=($(cat $XFILE_DIR/$XFILE))
SAMPLE_ID=${names[${SLURM_ARRAY_TASK_ID}]}

SUBSET_SPADES=${OUT_CHECKV_GENOMAD}/${SAMPLE_ID}/subset_spades.fasta
CLUSTER_RES=${OUT_DEREP}/${SAMPLE_ID}/clusterRes
TMP=${OUT_DEREP}/${SAMPLE_ID}/tmp

#load environment
CONDA="/groups/bhurwitz/miniconda3"
source $CONDA/etc/profile.d/conda.sh
conda activate mmseqs2_env

mkdir ${OUT_DEREP}/${SAMPLE_ID}

mmseqs easy-cluster ${SUBSET_SPADES} ${CLUSTER_RES} ${TMP} --min-seq-id 0.99 -c 0.90 --cov-mode 1

awk '/^>/{if($0!=prev){print; prev=$0}} !/^>/' ${OUT_DEREP}/${SAMPLE_ID}/clusterRes_all_seqs.fasta > ${OUT_DEREP}/${SAMPLE_ID}/cleaned_clusterRes_all_seqs.fasta
