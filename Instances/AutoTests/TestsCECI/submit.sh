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
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.8-3-1-Said-3.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.8-3-1-500-Said-6.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.8-3-1-100-Said-8.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.8-3-1-100-Said-3.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.8-3-1-100-Said-2.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.6-3-1-100-Said-8.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.4-3-1-Said-4.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.4-3-1-100-Said-5.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.4-3-1-100-Said-4.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.2-3-1-Said-4.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.2-3-1-Said-2.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.2-3-1-100-Said-6.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt
srun python AutoTest.py V-PT-50-0.2-3-1-100-Said-3.txt
srun ./RobustKEPSolver 0Config.txt
srun ./RobustKEPSolverOld 0Config.txt

