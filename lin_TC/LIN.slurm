#!/bin/bash
#SBATCH --account=cad16836
#SBATCH --constraint=GENOA
#SBATCH --job-name=Opti_lin_TC
#SBATCH --nodes=4
#SBATCH --ntasks=384
#SBATCH --time=1-00:00:00
#SBATCH --output=opti_SSP%j.out
#SBATCH --error=error_opti_SSPt%j.out
#SBATCH --exclusive


cd ~/tchebycube/lin_TC
source ~/./start.sh
module load anaconda3
conda activate Opti_1.8
set -x
 
python3 optimize_LIN.py --np $SLURM_NTASKS --input_jac run_TC.in --input_fun run_TC_fun.in --iter 10 --output Result --exe_dir build --E0 1E-3 --tol 0.01 -v
