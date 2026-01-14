#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --output=result.out
#SBATCH --nodes=3
#SBATCH --ntasks=6
#SBATCH --ntasks-per-node=2
#
sbcast -f test.py /tmp/test.py
srun python3 /tmp/test.py