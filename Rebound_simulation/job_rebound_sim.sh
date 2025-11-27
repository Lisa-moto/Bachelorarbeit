#!/bin/bash
#SBATCH --job-name=rebound-sim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --mem=4G
#SBATCH --output=rebound_sim.%j.out
#SBATCH --error=rebound_sim.%j.err

set -euo pipefail
cd /home/tu/tu_tu/tu_zxown31/Rebound_simulation

# direct python from your env (replace with the real path)
# this avoids activating conda/venv in the job
/home/tu/tu_tu/tu_zxown31/.conda/envs/main/bin/python3 rebound_sim.py