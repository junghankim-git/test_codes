#!/bin/bash

#draw_cs -v ps -v u10m /scratch/jhkim/TestBed/KIM/Output/3.0/ne30/gnu/10h/SW_G/101/UP-20110725220000-000010.nc -f org.png

#draw_cs -v ps -v u10m -l ../../scrip/3.remap_matrix/CS_HYCOM/ne030_rotated_to_hycom.nc -x dst_grid_center_lon -y dst_grid_center_lat result.nc


# BACK
draw_cs -v sst -v u10m ./result_cs.nc
