#!/usr/bin/env python

# load My Class
import sys
sys.path.append('./post')
from class_lonslats import *

infilename = './lonslats.nc'
cs = lonslats_plot(infilename)

#infilename = '/scratch/jhkim/TestBed/KIM/Output/2.5.16/ne30/gnu/10h/SW_G/101/UP-20110725130000-000001.nc'
#cs = lonslats_plot(infilename, False)

cs.draw('test', True)


cs.show()
