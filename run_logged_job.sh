#!/bin/bash

# Usage: ./run_logged_job.sh <script_to_run.py> [optional args]

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <python_script> [script_args...]"
    exit 1
fi

SCRIPT="$1"
shift
SCRIPT_ARGS="$@"
CSV_FILE="job_log.csv"

START_TIME=$(date +"%Y-%m-%d %H:%M:%S")

# === Run the script and capture return code ===
echo "Running: python $SCRIPT $SCRIPT_ARGS"
python "$SCRIPT" $SCRIPT_ARGS

RET_CODE=$?
END_TIME=$(date +"%Y-%m-%d %H:%M:%S")
# === Extract metadata ===

# Default values
NODES="N/A"
CPUS="N/A"
MEM="N/A"
CONFIG_TYPE="unknown"

# Check if Parsl config is embedded in script (no runtime capture needed)
if grep -q "SlurmProvider" "$SCRIPT"; then
    CONFIG_TYPE="parsl-embedded"

    NODES=$(grep -Po "nodes_per_block\s*=\s*\K\d+" "$SCRIPT" | head -1)
    CPUS=$(grep -Po "cores_per_node\s*=\s*\K\d+" "$SCRIPT" | head -1)
    MEM=$(grep -Po "mem_per_node\s*=\s*\K\d+" "$SCRIPT" | head -1)

    # Fallbacks if any value was missed
    NODES=${NODES:-"N/A"}
    CPUS=${CPUS:-"N/A"}
    MEM=${MEM:-"N/A"}
elif [ ! -z "$SLURM_JOB_ID" ]; then
    CONFIG_TYPE="slurm"
    NODES=$SLURM_JOB_NUM_NODES
    CPUS=$SLURM_CPUS_ON_NODE
    MEM=$((SLURM_MEM_PER_CPU * SLURM_CPUS_ON_NODE / 1000))  # in GB
fi

# === Create CSV header if not present ===
if [ ! -f "$CSV_FILE" ]; then
    echo "start_time,end_time,script,config_type,nodes,cpus,mem_gb,exit_code" > "$CSV_FILE"
fi

# === Append to log ===
echo "$START_TIME,$END_TIME,$SCRIPT,$CONFIG_TYPE,$NODES,$CPUS,$MEM,$RET_CODE" >> "$CSV_FILE"

echo "Logged run to $CSV_FILE"


