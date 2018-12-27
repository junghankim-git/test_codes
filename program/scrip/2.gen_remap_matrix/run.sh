#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1
#PBS -N SCRIP
#PBS -q serialq
#PBS -V
#PBS -j oe


EXPDIR=/home/jhkim/Study/Library/Main/SCRIP/2.gen_remap_matrix

cd $EXPDIR

./scrip > log.txt
