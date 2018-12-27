#!/usr/bin/env python

# load My Class
import sys
#sys.path.append('./post')
from class_wav_grid import *

os.system('rm -rf ./FigCube')
os.system('mkdir ./FigCube')

infilename = './cs.nc'

cs = wav_grid()

#cs.draw_grid_structure('aa')
cs.draw_all()

cs.show()
