#!/bin/bash
#PBS -l select=5:ncpus=16:mpiprocs=16:ompthreads=1
#PBS -N test
#PBS -q normal
#PBS -V
#PBS -j oe

#100#PBS -l select=7:ncpus=16:mpiprocs=16:ompthreads=1
#075#PBS -l select=5:ncpus=16:mpiprocs=16:ompthreads=1
#051#PBS -l select=4:ncpus=16:mpiprocs=16:ompthreads=1
#050#PBS -l select=4:ncpus=16:mpiprocs=16:ompthreads=1
#021#PBS -l select=2:ncpus=16:mpiprocs=16:ompthreads=1
#020#PBS -l select=2:ncpus=16:mpiprocs=16:ompthreads=1

export OMP_NUM_THREADS=1
export KMP_AFFINITY=disabled
export OMP_STACKSIZE=2048M
EXPDIR=/home/jhkim/work/program/task_control
cd $EXPDIR
rm -rf ./exp.log
touch ./exp.log

mpirun -np 75 -f $PBS_NODEFILE ./run.exe >> ./exp_075.log 2>&1
