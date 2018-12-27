#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1
#PBS -N remap
##PBS -l place=scatter:excl
#PBS -q serialq
#PBS -V
#PBS -j oe


WORKDIR=/home/jhkim/work/program/remapper
cd $WORKDIR

./remap_all.py > ./exp.log 2>&1
