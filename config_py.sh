#defining the log and scripts directories
export WORK_DIR=/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline
export SCRIPT_DIR=/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline/run_scripts
export LOG_DIR=/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline/logs

# defining input from assembly 
export XFILE=xac
export XFILE_DIR=/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses
export SPADES_DIR=/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/out_spades

# variables used for 3 viromics tools as well as checkv 
#export OUT_VIRSORT=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/01B_checkv_virsorter
export CHECKVDB=/xdisk/bhurwitz/databases/checkv-db-v1.5  
#export OUT_DVF=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/02A_dvf
#export DVF_DB=/groups/bhurwitz/databases/DeepVirFinder
#export OUT_CHECKV_DVF=/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/02B_checkv_dvf
export OUT_GENOMAD=/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline/results_testing/03A_genomad
export GENOMAD_DB=/xdisk/bhurwitz/databases/genomad_db          
export OUT_CHECKV_GENOMAD=/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline/results_testing/03B_checkv_genomad
export CHECKV_PARSER=/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline/run_scripts/CheckV_parser.R 
export PARSE_LENGTH=5000
#export RSCRIPT_DIR=/groups/bhurwitz/miniconda3/bin/Rscript  '''redownload'''

# dereplication and clustering 
export OUT_DEREP=/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline/results_testing/04A_dereplicate
export OUT_CLUSTER=/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline/results_testing/04B_cluster

# step 1 create blastdb
export DB_DIR=/xdisk/bhurwitz/databases/AVrC    
export MAX_DB_SIZE="0.5GB" 

# step 2 : blast query against blast db
export FASTA_DIR=/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline/query
export FA_SPLIT_FILE_SIZE=5000000 # in bytes, 5000 in KB

# BLAST parameters
export BLAST_TYPE=blastn
export MAX_TARGET_SEQS=1
export EVAL=1e-3
export OUT_FMT=6 # tabular format with no headings

# Annotation parameters
export PCTID=85
export LENGTH=1000
export BLAST_HITS="/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline/results_testing/05D_mergeblast/AVrC_allrepresentatives.fasta/clusterRes_rep_seq.fasta.txt"
export ANNOTATIONS="/xdisk/bhurwitz/databases/AVrC/database_csv/"
export OUTPUT="/xdisk/bhurwitz/virus_hunting/kolodisner/fmt_viruses/viral_detection_pipeline/results_testing/06_annotate"


#
# Some custom functions for our scripts
#
''' --------------------------------------------------
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
}'''
