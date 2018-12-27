#!/usr/bin/env python

# load My Class
#sys.path.append('./')
from class_scrip_in import *


infilename = './HYCOM/hycom_0240x0193.nc'
title      = 'hycom (240x193)'
infilename = './CS/khkim_ne030_rotated.nc'
title      = 'cubed-sphere (ne030_r)'
remap = class_remap_in(infilename)


# hycom
# type: 0(oth), 1(sphere), 2(local)
#remap.draw('ne030(rotated) to hycom(240x193)')
#remap.draw('ne030(rotated) to hycom(240x193)',ix=42001)

# cubed-sphere
remap.draw_corners(title,175,175)

remap.show()
