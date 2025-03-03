#!/bin/bash
#SBATCH --job-name=malathion_analysis
#SBATCH --output=malathion_%j.out
#SBATCH --error=malathion_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=40

# Load R module (replace ?? with the correct version, e.g., 4.2.1)
module load R/4.2.1

# Run the scripts
Rscript problem1.R
Rscript problem2.R
Rscript problem3.R

# Check efficiency (optional)
seff $SLURM_JOB_ID
