#!/bin/bash
#
#SBATCH --job-name=test_kep
#SBATCH --output=res.txt
#
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=4096

srun python AutoTest.py V-PT-50-0.8-3-1-Said-8.txt
srun ./RobustKEPSolver 0Config.txt

