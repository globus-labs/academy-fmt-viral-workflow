#!/bin/bash
#SBATCH --output=./logs/04B_cluster/Job-%a.out
#SBATCH --account=bhurwitz
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5000

pwd; hostname; date
source ./config.sh

CLUSTER_RES=${OUT_CLUSTER}/clusterRes
TMP=${OUT_CLUSTER}/tmp

#load environment
CONDA="/groups/bhurwitz/miniconda3"
source $CONDA/etc/profile.d/conda.sh
conda activate mmseqs2_env

cat ${OUT_DEREP}/*/cleaned_clusterRes_all_seqs.fasta > ${OUT_DEREP}/dereplicated.fasta

mmseqs easy-cluster ${OUT_DEREP}/dereplicated.fasta ${CLUSTER_RES} ${TMP} --min-seq-id 0.95 -c 0.75 --cov-mode 1

cp ${OUT_CLUSTER}/clusterRes_rep_seq.fasta ${WORK_DIR}/query/
