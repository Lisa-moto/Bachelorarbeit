#!/bin/bash
#SBATCH --job-name=sim-with-moon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=14-00:00:00
#SBATCH --partition=compute
#SBATCH --mem=8G

set -euo pipefail

cd /home/tu/tu_tu/tu_zxown31/Rebound_simulation

# Use the Python interpreter from the existing conda environment directly.
/home/tu/tu_tu/tu_zxown31/.conda/envs/main/bin/python3 sim_with_moon.py >> output.txt