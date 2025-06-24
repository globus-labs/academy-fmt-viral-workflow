#!/bin/bash
#SBATCH --output=./logs/pipeline/Job-%j.out
#SBATCH --error=./logs/pipeline/Job-%j.err
#SBATCH --account=bhurwitz
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5000


python viral_detection_pipeline.py
