#!/usr/bin/env python

# load My Class
import sys
SHARE_DIR='/home/jhkim/work/share/python'
sys.path.append(SHARE_DIR)
from class_cubed_sphere import *


infilename = './cs.nc'

cs = cubed_sphere(infilename)

cs.draw_one_elem('elem')

cs.show()
