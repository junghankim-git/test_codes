#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1
#PBS -N main_qsub
#PBS -q serialq
#PBS -V
#PBS -j oe

exp_dir="/home/jhkim/work/program/multi_io"
out_dir="/scratch/jhkim/data/multi_io"

cd $exp_dir

./master.sh $exp_dir $out_dir
