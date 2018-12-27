#!/usr/bin/env python

# load My Class
import sys
SHARE_DIR='../../share/python'
sys.path.append(SHARE_DIR)
from class_cubed_sphere import *

#os.system('rm -rf ./fig')
#os.system('mkdir ./fig')


infilename = './cs.nc'

cs = cubed_sphere(infilename)

cs.draw_cube_structure('test1')       # 6 faces & face number, quadrature points
cs.draw_cube_cart_glb_index('test2')  # catecian coordinates of cube, element numbers, 
cs.draw_cube_glb_sfc_index('test3')   # element numbers, sfc connection lines
cs.draw_sfc_rank_dist('test4')        # sfc number, rank number
cs.draw_domain('test5')               # local index, neiborhoods at specific rank
cs.draw_local_domain('test6')               # local index, neiborhoods at specific rank
#cs.draw_cubed_sphere('figure_1_1',True)         # element line at sphere view

cs.show()
