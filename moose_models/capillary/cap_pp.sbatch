#!/bin/bash -l
#SBATCH --job-name=cap_pp
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#### Following is 100GB/node, which is optimal for pearcey
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrew.wilkins@csiro.au
##SBATCH --export=NONE

source ~/.profile

export CSIRO_DIR=/projects/moose/projects_wil04q/csiro
export SIMULATE_DIR=/home/wil04q/anlec/capillary_rerun_radio

ulimit -s 102400

cd $SIMULATE_DIR
mpirun -n 20 ${CSIRO_DIR}/csiro-opt -i cap_pp.i > ${SIMULATE_DIR}/cap_pp.log 2>&1
