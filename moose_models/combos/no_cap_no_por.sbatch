#!/bin/bash -l
#SBATCH --job-name=no_cap_no_por
#SBATCH --time=0:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#### Following is 100GB/node, which is optimal for pearcey
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrew.wilkins@csiro.au
##SBATCH --export=NONE

source ~/.profile
ulimit -s 102400

export PORFLOW_DIR=/projects/cmrp/projects_wil04q/moose/modules/porous_flow

cd $SLURM_SUBMIT_DIR
mpirun -n 20 ${PORFLOW_DIR}/porous_flow-opt -i no_cap_no_por.i > ${SLURM_SUBMIT_DIR}/no_cap_no_por.log 2>&1
