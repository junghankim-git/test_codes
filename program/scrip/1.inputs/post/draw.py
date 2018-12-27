#!/usr/bin/env python

# load My Class
#sys.path.append('./')
from SCRIP_Input import *

#os.system('rm -rf ./FigCube')
#os.system('mkdir ./FigCube')



#infilename = '../CS/khkim_ne030_rotated.nc'
infilename = '../HYCOM/hycom_0240x0193.nc'
#infilename = '../HYCOM/hycom_0240x0193_wo_pole.nc'
#infilename = '/home/jhkim/work/program/gen_remap_mat_with_latlon/scrip_1.nc'

cs = SCRIP_Input(infilename)


cs.Draw('HYCOM grid (240x193)', False)


cs.Show()
