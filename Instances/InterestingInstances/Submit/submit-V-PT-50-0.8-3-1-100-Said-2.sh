#!/bin/bash
#
#SBATCH --job-name=kep-V-PT-50-0.8-3-1-100-Said-2
#SBATCH --output=res-V-PT-50-0.8-3-1-100-Said-2.txt
#
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=4096

srun ./RobustKEPSolver 0Config-V-PT-50-0.8-3-1-100-Said-2.txt
srun ./RobustKEPSolverOld 0Config-V-PT-50-0.8-3-1-100-Said-2.txt