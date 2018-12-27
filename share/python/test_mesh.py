
import os
import sys
import numpy as np
import traceback
sys.path.append('/home/jhkim/work/share/Python')
from control_mesh import *


mesh = [[1,3],[2,4]]
print_mesh_2d(2,2,mesh,'original')

cmesh = symmetry_mesh(2,0,mesh)
print_mesh_2d(2,2,cmesh,'x-axis')

cmesh = symmetry_mesh(2,1,mesh)
print_mesh_2d(2,2,cmesh,'y-axis')

cmesh = symmetry_mesh(2,2,mesh)
print_mesh_2d(2,2,cmesh,'xy-axis')

cmesh = symmetry_mesh(2,3,mesh)
print_mesh_2d(2,2,cmesh,'mxy-axis')

cmesh = symmetry_mesh(2,4,mesh)
print_mesh_2d(2,2,cmesh,'0-axis')





mesh = [[0,2],[1,3]]
print_mesh_2d(2,2,mesh,'original',dir=True)

cmesh = symmetry_dir_mesh(2,0,mesh)
print_mesh_2d(2,2,cmesh,'x-axis',dir=True)

cmesh = symmetry_dir_mesh(2,1,mesh)
print_mesh_2d(2,2,cmesh,'y-axis',dir=True)

cmesh = symmetry_dir_mesh(2,2,mesh)
print_mesh_2d(2,2,cmesh,'xy-axis',dir=True)

cmesh = symmetry_dir_mesh(2,3,mesh)
print_mesh_2d(2,2,cmesh,'mxy-axis',dir=True)

cmesh = symmetry_dir_mesh(2,4,mesh)
print_mesh_2d(2,2,cmesh,'0-axis',dir=True)

