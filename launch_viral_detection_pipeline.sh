#! /bin/bash
# script dir
export SCRIPT_DIR="./run_scripts"

# Unzip - no dependencies
job1=$(sbatch $SCRIPT_DIR/unzip.sh)
$jid1=$(echo $job1 | sed 's/^Submitted batch job //')
echo $jid1 >> job_ids

# Virsorter2 - 1 dependency
#job2=$(sbatch  --dependency=afterany:$jid1 $SCRIPT_DIR/01A_virsorter2.sh)
#jid2=$(echo $job2 | sed 's/^Submitted batch job //')
#echo $jid2 >> job_ids

# DeepVirFinder - 1 dependency
#job3=$(sbatch  --dependency=afterany:$jid1 $SCRIPT_DIR/02A_dvf.sh)
#jid3=$(echo $job3 | sed 's/^Submitted batch job //')
#echo $jid3 >> job_ids

# Genomad - 1 dependency
job4=$(sbatch  --dependency=afterany:$jid1 $SCRIPT_DIR/03A_genomad.sh)
jid4=$(echo $job4 | sed 's/^Submitted batch job //')
echo $jid4 >> job_ids

# CheckV Virsorter2 - 1 dependency
#job5=$(sbatch  --dependency=afterok:$jid2 $SCRIPT_DIR/01B_checkv_virsorter2.sh)
#jid5=$(echo $job5 | sed 's/^Submitted batch job //')
#echo $jid5 >> job_ids

# CheckV DeepVirFinder - 1 dependency
#job6=$(sbatch  --dependency=afterok:$jid3 $SCRIPT_DIR/02B_checkv_dvf.sh)
#jid6=$(echo $job6 | sed 's/^Submitted batch job //')
#echo $jid6 >> job_ids

# CheckV Genomad - 1 dependency
job7=$(sbatch  --dependency=afterok:$jid4 $SCRIPT_DIR/03B_checkv_genomad.sh)
jid7=$(echo $job7 | sed 's/^Submitted batch job //')
echo $jid7 >> job_ids

# Zip - 1 dependency
job8=$(sbatch  --dependency=afterok:$jid7 $SCRIPT_DIR/zip.sh)
jid8=$(echo $job8 | sed 's/^Submitted batch job //')
echo $jid8 >> job_ids

# Dereplicate - 1 dependency
job9=$(sbatch  --dependency=afterok:$jid7 $SCRIPT_DIR/04A_dereplicate.sh)
jid9=$(echo $job9 | sed 's/^Submitted batch job //')
echo $jid9 >> job_ids

# Cluster - 1 dependency
job10=$(sbatch  --dependency=afterok:$jid9 $SCRIPT_DIR/04B_cluster.sh)
jid10=$(echo $job10 | sed 's/^Submitted batch job //')
echo $jid10 >> job_ids

# 05A_makeblastdb.sh - create the blast databases - no dependencies
#job11=$(sbatch --dependency=afterok:$jid10 $SCRIPTS_DIR/05A_makeblastdb.sh)
#jid11=$(echo $job11 | sed 's/^Submitted batch job //')
#echo $jid11

# 05B_launchblast.sh - jid12 depends on jid11
# This script: 
# 1. splits the query files (into small chunks)
# 2. runs job 05C_blast.sh to blast each chunk vs the databases
# 3. runs jub 05D_mergeblast.sh to collate results by input file
#job12=$(sbatch --dependency=afterok:$jid11 $SCRIPTS_DIR/05B_launchblast.sh)
job12=$(sbatch --dependency=afterok:$jid10 $SCRIPTS_DIR/05B_launchblast.sh)
jid12=$(echo $job12 | sed 's/^Submitted batch job //')
echo $jid12

# show dependencies in squeue output:
squeue -u $USER -o "%.8A %.4C %.10m %.20E"
