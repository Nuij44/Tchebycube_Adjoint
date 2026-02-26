#!/bin/bash
#SBATCH --account=cad17050
#SBATCH --constraint=GENOA
#SBATCH --job-name=2Z_3000
#SBATCH --nodes=4
#SBATCH --ntasks=384
#SBATCH --time=1-00:00:00
#SBATCH --output=opti_SSP%j.out
#SBATCH --error=error_opti_SSPt%j.out
#SBATCH --exclusive


cd ~/Tchebycube_Adjoint/Taylor_Couette
source ~/./start.sh
module load anaconda3
conda activate Opti_1.8
set -x
 
python3 optimize_TC.py --np $SLURM_NTASKS --input_jac run_TC.in --input_fun run_TC_fun.in --iter 10 --output Result --exe_dir build --E0 1E-5 --tol 0.01 -v
