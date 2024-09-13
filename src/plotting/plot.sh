#!/bin/bash
#SBATCH --job-name=mcmc
#SBATCH --partition=private-dpt-cpu
#SBATCH --time=00-01:00:00
#SBATCH --output=/home/users/r/radhakrb/cardy/data/slurm-output/slurm-outfiles/%j.out
#SBATCH --error=/home/users/r/radhakrb/cardy/data/slurm-output/slurm-errors/%j.err
#SBATCH -N 1 # total number of nodes
#SBATCH --mem=30G

srun python3 /home/users/r/radhakrb/cardy/src/plotting/evdensity.py $@