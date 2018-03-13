#!/bin/bash
#
#SBATCH --job-name=kep-V-PT-50-0.2-3-1-Said-4
#SBATCH --output=res-V-PT-50-0.2-3-1-Said-4.txt
#
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=4096

srun ./RobustKEPSolver 0Config-V-PT-50-0.2-3-1-Said-4.txt
srun ./RobustKEPSolverOld 0Config-V-PT-50-0.2-3-1-Said-4.txt