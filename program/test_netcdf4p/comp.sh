#!/bin/bash

#VER=4.6.2
VER=4.7.1

#NETCDF_PATH=/app/compilers/gcc/4.6.2/mpi/mvapich2/1.8.1/applib2/NETCDF4/4.1.3
NETCDF_PATH=/app/compilers/gcc/$VER/mpi/mvapich2/1.8.1/applib2/NETCDF4/4.1.3
#NETCDF_PATH=/app/compilers/gcc/$VER/applib1/NETCDF4/4.1.3

HDF5_PATH=/app/compilers/gcc/4.7.1/mpi/mvapich2/1.8.1/applib2/HDF5/1.8.9

echo mpif90 -O2 -ffree-line-length-none -o run.exe main_program.F90 kinds.o mod_netcdf.o \
            -L$NETCDF_PATH/lib -lnetcdff -lnetcdf -I$NETCDF_PATH/include            \
            -L$HDF5_PATH/lib -lhdf_hl -lhdf       -I$HDF5_PATH/include
mpif90 -O2 -ffree-line-length-none -o run.exe main_program.F90 kinds.o mod_netcdf.o \
            -L$NETCDF_PATH/lib -lnetcdff -lnetcdf -I$NETCDF_PATH/include            \
            -L$HDF5_PATH/lib -lhdf_hl -lhdf       -I$HDF5_PATH/include





#mpif90 -O2 -ffree-line-length-none -o run.exe main_program.F90 kinds.o mod_netcdf.o \
#            $NETCDF_PATH/lib/libnetcdff.a   \
#            $NETCDF_PATH/lib/libnetcdf.so   \
#          -I$NETCDF_PATH/include \
#            $NETCDF_PATH/lib/libnetcdff.a   \
#            $NETCDF_PATH/lib/libnetcdf.so   \
#          -I$NETCDF_PATH/include
