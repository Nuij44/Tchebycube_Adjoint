#!/bin/bash
#SBATCH --account=circe
#SBATCH --constraint=amd9124
#SBATCH --job-name=visu
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --time=5:00:00
#SBATCH --output=opti_SSP%j.out
#SBATCH --error=error_visu%j.out
#SBATCH --exclusive

cd ~/tchebycube/lin_TC
source ~/./start.sh

set -x

rm -rf outdir/ && srun -n 32 --cpu-bind=threads -- build/visu_TC.x --input ./run_visu.in

