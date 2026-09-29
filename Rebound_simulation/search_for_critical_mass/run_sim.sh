#!/bin/bash
#SBATCH --job-name=toi178_moon
#SBATCH --partition=compute
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --array=0-199
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

module load devel/python/3.11

python main_sim.py "${SLURM_ARRAY_TASK_ID}"