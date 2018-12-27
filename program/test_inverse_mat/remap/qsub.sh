#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1
#PBS -q normal
#PBS -N remap
#PBS -j oe
#PBS -V


export OMP_NUM_THREADS=1
export KMP_AFFINITY=disabled
export OMP_STACKSIZE=2048M

cd /home/jhkim/work/program/test_inverse_mat/remap

./run.exe > exp.log 2>&1
