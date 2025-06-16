#defining the log and scripts directories
export WORK_DIR=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline
export SCRIPT_DIR=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/run_scripts
export LOG_DIR=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/logs

# defining input from assembly 
export XFILE=xaa
export XFILE_DIR=/xdisk/bhurwitz/virus_hunting/data/all_assemblies
export SPADES_DIR=/xdisk/bhurwitz/virus_hunting/data/all_assemblies/out_spades

# variables used for 3 viromics tools as well as checkv 
#export OUT_VIRSORT=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/01B_checkv_virsorter
export CHECKVDB=/groups/bhurwitz/databases/checkv-db-v1.5
#export OUT_DVF=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/02A_dvf
#export DVF_DB=/groups/bhurwitz/databases/DeepVirFinder
#export OUT_CHECKV_DVF=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/02B_checkv_dvf
export OUT_GENOMAD=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/03A_genomad
export GENOMAD_DB=/groups/bhurwitz/databases/genomad_db
export OUT_CHECKV_GENOMAD=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/03B_checkv_genomad
export CHECKV_PARSER=/groups/bhurwitz/tools/R_scripts/CheckV_parser.R
export PARSE_LENGTH=5000
export RSCRIPT_DIR=/groups/bhurwitz/miniconda3/bin/Rscript

# dereplication and clustering 
export OUT_DEREP=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/04A_dereplicate
export OUT_CLUSTER=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/04B_cluster

# step 1 create blastdb
export DB_DIR=/xdisk/bhurwitz/databases/AVrC
export MAX_DB_SIZE="0.5GB" 

# step 2 : blast query against blast db
export FASTA_DIR=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/query
export FA_SPLIT_FILE_SIZE=5000000 # in bytes, 5000 in KB

# containers for tools
export FASPLIT=/contrib/singularity/shared/bhurwitz/ucsc-fasplit:469--h9b8f530_0.sif
export BLAST=/contrib/singularity/shared/bhurwitz/blast:2.16.0--hc155240_2.sif

# BLAST parameters
export BLAST_TYPE=blastn
export MAX_TARGET_SEQS=1
export EVAL=1e-3
export OUT_FMT=6 # tabular format with no headings

#
# Some custom functions for our scripts
#
# --------------------------------------------------
function init_dir {
    for dir in $*; do
        if [ -d "$dir" ]; then
            rm -rf $dir/*
        else
            mkdir -p "$dir"
        fi
    done
}

# --------------------------------------------------
function create_dir {
    for dir in $*; do
        if [[ ! -d "$dir" ]]; then
          echo "$dir does not exist. Directory created"
          mkdir -p $dir
        fi
    done
}

# --------------------------------------------------
function lc() {
    wc -l $1 | cut -d ' ' -f 1
}
