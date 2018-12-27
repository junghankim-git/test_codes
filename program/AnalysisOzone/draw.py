#!/usr/bin/env python

# load My Class
#sys.path.append('./')
from Checker import *

varname = 'oz'

# org
#kimfilepath  = '/scratch/jhkim/TestBed/Data/2.3.03_2/10h/SW_G/101_org'
#kimfilepath  = '/scratch/jhkim/TestBed/Data/2.3.02_03/10h/SW_G/101'
#kimfilepath  = '/scratch/jhkim/TestBed/Data/2.3.03_2/10h/SW_G/101_org'
#kimfilepath  = '/scratch/jhkim/TestBed/Data/2.3.09_01/10h/SW_G/101'
#kimfilepath  = '/scratch/jhkim/TestBed/Data/3.0.00.01/10h/SW_G/101'
kimfilepath  = '/scratch/jhkim/TestBed/Data/3.0.01/10h/SW_G/101'
kimfilenames =  [\
   'UP-20110725120000-000001.nc', \
   'UP-20110725130000-000002.nc', \
   'UP-20110725140000-000003.nc', \
   'UP-20110725150000-000004.nc', \
   'UP-20110725160000-000005.nc', \
   'UP-20110725170000-000006.nc', \
   'UP-20110725180000-000007.nc', \
   'UP-20110725190000-000008.nc', \
   'UP-20110725200000-000009.nc', \
   'UP-20110725210000-000010.nc', \
   'UP-20110725220000-000011.nc' ]


check = Checker(kimfilepath, kimfilenames, varname)

check.Draw(ilev=1,  type=0)
check.Draw(ilev=25, type=0)
check.Draw(ilev=50, type=0)

check.Show()
