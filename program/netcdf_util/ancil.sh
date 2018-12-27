#!/bin/bash


INPUT_PATH=/data/KIM3.3/inputdata/ne768np3_rotated
OUTPUT_PATH=out_ancil

IFILES=()
IFILES[0]="clim_aerosol_BC.nc.NETCDF4"
IFILES[1]="clim_aerosol_DUST.nc.NETCDF4"
IFILES[2]="clim_aerosol_OC.nc.NETCDF4"
IFILES[3]="clim_aerosol_SEASALT.nc.NETCDF4"
IFILES[4]="clim_aerosol_SO4.nc.NETCDF4"
IFILES[5]="clim_ozone.nc.NETCDF4"
OFILES=()
OFILES[0]="clim_aerosol_BC.nc"
OFILES[1]="clim_aerosol_DUST.nc"
OFILES[2]="clim_aerosol_OC.nc"
OFILES[3]="clim_aerosol_SEASALT.nc"
OFILES[4]="clim_aerosol_SO4.nc"
OFILES[5]="clim_ozone.nc"


for ((i=1; i<6; i++))
do
echo " "
echo "nc_util -m 2 $INPUT_PATH/${IFILES[$i]} -o $OUTPUT_PATH/${OFILES[$i]}"
nc_util -m 2 $INPUT_PATH/${IFILES[$i]} -o $OUTPUT_PATH/${OFILES[$i]}
done

#> nc_util -m 2 $infile -o $oufile

#nc_util -m 2 /data/KIM3.3/inputdata/ne768np3_rotated/clim_aerosol_BC.nc.NETCDF4 -o out_ancil/clim_aerosol_BC.nc
