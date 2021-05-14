#!/bin/bash -l
#SBATCH --job-name=radial_thermal
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#### Following is 100GB/node, which is optimal for pearcey
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrew.wilkins@csiro.au
##SBATCH --export=NONE

source ~/.profile

export CSIRO_DIR=/datasets/work/en-cmrp/work/apps/projects/projects_wil04q/moose/modules/porous_flow
export SIMULATE_DIR=/datasets/work/en-cmrp/work/apps/projects/projects_wil04q/tara_cf_moose_models/moose_models/radial_thermal

ulimit -s 102400

cd $SIMULATE_DIR
mpirun -n 20 ${CSIRO_DIR}/porous_flow-opt -i radial_thermal.i > ${SIMULATE_DIR}/radial_thermal.log 2>&1
