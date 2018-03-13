#!/bin/bash
#
#SBATCH --job-name=kep-V-PT-50-0.8-3-1-500-Said-6
#SBATCH --output=res-V-PT-50-0.8-3-1-500-Said-6.txt
#
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=4096

srun ./RobustKEPSolver 0Config-V-PT-50-0.8-3-1-500-Said-6.txt
srun ./RobustKEPSolverOld 0Config-V-PT-50-0.8-3-1-500-Said-6.txt