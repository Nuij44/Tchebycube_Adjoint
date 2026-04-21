#!/bin/bash
#SBATCH --account=circe
#SBATCH --constraint=amd9124
#SBATCH --job-name=Opti_TC
#SBATCH --nodes=10
#SBATCH --ntasks=320
#SBATCH --ntasks-per-node=32
#SBATCH --time=3-0:00:00
#SBATCH --output=opti_SSP%j.out
#SBATCH --error=error_opti_SSPt%j.out
#SBATCH --exclusive

source ~/./start.sh
cd ~/tchebycube/Taylor_Couette

set -x
module load miniconda/25.1.1

source /softs/miniconda3/etc/profile.d/conda.sh
conda activate Opti_1.8

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

python3 V2optimize_TC.py --np $SLURM_NTASKS --input_jac run_TC.in --input_fun run_TC_fun.in --iter 10 --output Result --exe_dir build --E0 1E-4 --tol 0.01 -v --armijo #--multi 
