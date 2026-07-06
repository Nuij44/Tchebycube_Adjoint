#!/bin/bash
#SBATCH --account=project_465002977
#SBATCH --partition=small
#SBATCH --job-name=Q1_Lesur
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --output=Q1_Lesur%j.out
#SBATCH --error=error_Q1_Lesur%j.out
#SBATCH --exclusive


cd ~/Tchebycube_Adjoint/turbulent_TC/Q1_Lesur_Longaretti
source ~/./start.sh

set -x
 
srun -n 128 ../build/Q1_kerswell.x --input run_cart.in
