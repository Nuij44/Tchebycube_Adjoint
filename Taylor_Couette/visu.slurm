#!/bin/bash
#SBATCH --account=circe
#SBATCH --constraint=amd9124
#SBATCH --job-name=visu
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=32
#SBATCH --time=5:00:00
#SBATCH --output=visu%j.out
#SBATCH --error=error_visu%j.out
#SBATCH --exclusive

cd ~/tchebycube/Taylor_Couette
source ~/./start.sh

set -x

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

rm -rf outdir/ && srun -n $SLURM_NTASKS ./build/visu_TC.x --input ./run_visu.in

