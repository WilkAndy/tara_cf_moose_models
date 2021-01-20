#!/bin/bash -l
#SBATCH --job-name=rfp
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#### Following is 100GB/node, which is optimal for pearcey
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrew.wilkins@csiro.au
##SBATCH --export=NONE

source ~/.profile

export PORFLOW_DIR=/projects/cmrp/projects_wil04q/moose/modules/porous_flow
export SIMULATE_DIR=/home/wil04q/cemvmc_home/papers/tara_cf/moose_models/rfp

ulimit -s 102400

cd $SIMULATE_DIR
mpirun -n 20 ${PORFLOW_DIR}/porous_flow-opt -i rfp.i > ${SIMULATE_DIR}/rfp.log 2>&1
