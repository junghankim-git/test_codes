#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1
#PBS -N SCRIP
#PBS -q serialq
#PBS -V
#PBS -j oe


EXPDIR=/home/jhkim/work/program/GenRemapMat_withLatLon

cd $EXPDIR

./run.exe > log.txt 2>&1

./scrip  >> log.txt 2>&1
