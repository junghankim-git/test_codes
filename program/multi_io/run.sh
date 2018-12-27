#!/bin/bash
#PBS -l select=7:ncpus=16:mpiprocs=16:ompthreads=1
#PBS -N log_101
#PBS -q normal
#PBS -W block=true
#PBS -V
#PBS -j oe

cd /home/jhkim/work/program/multi_io

mpirun -np 101 -f $PBS_NODEFILE ./run.exe | tee exp.txt
