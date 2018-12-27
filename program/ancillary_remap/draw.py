#!/usr/bin/env python

import sys

# load My Class
sys.path.append('/home/jhkim/work/share/python')
from class_draw_remap import *


# remap file info
src_sizename = 'src_grid_size'
src_latname  = 'src_grid_center_lat'
src_lonname  = 'src_grid_center_lon'
dst_sizename = 'dst_grid_size'
dst_latname  = 'dst_grid_center_lat'
dst_lonname  = 'dst_grid_center_lon'


'''
#=============== test
# source file and target file(remaped file)
remapfile = '/home/jhkim/Study/Library/Main/GenRemapMat_withLatLon/remap_matrix/ne030rotated_ne030unrotated/remap_12.nc'
title     = 'test'
srcfile   = '/scratch/jhkim/TestBed/KIM/Output/2.4/ne30/gnu/10h/SW_G/101/UP-20110725220000-000011.nc'
dstfile   = '/home/jhkim/Study/Library/Main/Remapper/result.nc'
# variable informations
varnames  = ['tsfc', 'T']
ndims     = [1, 2]
#=============== test
'''


#remapfile = '/home/jhkim/Study/Library/Main/GenRemapMat_withLatLon/remap_matrix/ne120rotated_ne030variable/remap_12.nc'
remapfile = '/home/jhkim/work/program/gen_remap_mat_with_latlon/remap_matrix/ne030rotated_ne030variable/remap_12.nc'
#remapfile = '/home/jhkim/Study/Library/Main/GenRemapMat_withLatLon/remap_matrix/ne030rotated_ne030unrotated/remap_12.nc'


'''
#=============== gfrac.nc
# source file and target file(remaped file)
title     = 'clim_aerosol' # just title
srcfile   = '/data/KIM2.4/inputdata/ne120np4_rotated/gfrac.nc'  # IN
dstfile   = '/home/jhkim/Study/Library/Main/Remapper/outputs/gfrac.nc'  # OUT
varnames  = ['gfrac']
ndims     = [1]
#=============== gfrac.nc
'''


'''
#=============== clim_aerosol.nc
# source file and target file(remaped file)
title     = 'clim_aerosol' # just title
srcfile   = '/data/KIM2.4/inputdata/ne030np4_rotated/clim_aerosol.nc'  # IN
dstfile   = './outputs/clim_aerosol.nc'  # OUT
varnames  = ['BC', 'OC']
ndims     = [3, 3]
#varnames  = ['BC', 'OC', 'SO4', 'SEASALT', 'DUST']
#ndims     = [3, 3, 3, 3, 3]
#=============== clim_aerosol.nc
'''



remap = class_draw_remap(remapfile, src_sizename, src_latname, src_lonname, \
                         remapfile, dst_sizename, dst_latname, dst_lonname)
remap.set_variables(srcfile, dstfile, varnames, ndims)
remap.draw(title, 0)
