#!/bin/bash
#SBATCH --job-name=toi178_moon_around_d
#SBATCH --partition=compute
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --requeue

cd "$SLURM_SUBMIT_DIR"

module load devel/miniforge
conda activate main

export OMP_NUM_THREADS=1
python -u sim.py