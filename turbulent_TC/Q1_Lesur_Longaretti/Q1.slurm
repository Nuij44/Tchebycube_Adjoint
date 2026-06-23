#!/bin/bash
#SBATCH --account=c2017511
#SBATCH --constraint=GENOA
#SBATCH --job-name=Q1_Lesur
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --time=1-00:00:00
#SBATCH --output=Q1_Lesur%j.out
#SBATCH --error=error_Q1_Lesur%j.out
#SBATCH --exclusive


cd ~/Tchebycube_Adjoint/turbulent_TC/Q1_Lesur_Longaretti
source ~/./start.sh

set -x
 
srun -n 96 ../build/Q1_kerswell.x --input run_cart.in
