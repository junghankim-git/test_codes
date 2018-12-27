#!/usr/bin/env python

# load My Class
#sys.path.append('./')
from Checker import *


kimfilenames =  [\
   './output/UP-20110725120000-000001.nc', \
   './output/UP-20110725130000-000002.nc', \
   './output/UP-20110725140000-000003.nc', \
   './output/UP-20110725150000-000004.nc', \
   './output/UP-20110725160000-000005.nc', \
   './output/UP-20110725170000-000006.nc', \
   './output/UP-20110725180000-000007.nc', \
   './output/UP-20110725190000-000008.nc', \
   './output/UP-20110725200000-000009.nc', \
   './output/UP-20110725210000-000010.nc', \
   './output/UP-20110725220000-000011.nc' ]

#wavfilename   = '/home/jhkim/Study/Library/Main/SCRIP/4.check_wav_kim/output/file_wav.nc'
wavfilename   = '/home/jhkim/TestBed/KIM/2.3.micros/2.3.02/01.dev/Exp_10h/SW_G/file_wav.nc'
remapfilename = '/home/jhkim/Study/Library/Main/SCRIP/3.remap_matrix/CS_WW3/ne030_rotated_to_WW3.nc'

check = Checker(kimfilenames, wavfilename, remapfilename)

#check.Draw('u10m', type=0)
#check.Draw_grid_last('u10m', type=0)
check.Draw_flat('u10m', type=0)

check.Show()
