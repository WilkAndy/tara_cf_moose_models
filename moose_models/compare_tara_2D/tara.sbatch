#!/bin/bash -l
#SBATCH --job-name=tara2D
#SBATCH --time=00:40:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#### Following is 100GB/node, which is optimal for pearcey
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrew.wilkins@csiro.au
##SBATCH --export=NONE

source ~/.profile

export CSIRO_DIR=/projects/moose/projects_wil04q/csiro
export SIMULATE_DIR=/home/wil04q/anlec/compare_tara_2D_rerun_radio

ulimit -s 102400

cd $SIMULATE_DIR
mpirun -n 20 ${CSIRO_DIR}/csiro-opt -i tara.i > ${SIMULATE_DIR}/tara.log 2>&1
