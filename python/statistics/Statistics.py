#!/usr/bin/env python

import os
import sys
import math
SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
sys.path.append(SHARE_DIR)
from class_statistics import *


var1 = [1,2,3,4,5,6,7,8]
var2 = [[1,2,3,4],[5,6,7,8]]
var3 = [[[1,2],[3,4]],[[5,6],[7,8]]]
print 'L1 1D = ', get_lnorm(var1, 'L1')
print 'L1 2D = ', get_lnorm(var2, 'L1')
print 'L1 3D = ', get_lnorm(var3, 'L1')
print 'L2 1D = ', get_lnorm(var1, 'L2')
print 'L2 2D = ', get_lnorm(var2, 'L2')
print 'L2 3D = ', get_lnorm(var3, 'L2')
print 'Li 1D = ', get_lnorm(var1, 'Linf')
print 'Li 2D = ', get_lnorm(var2, 'Linf')
print 'Li 3D = ', get_lnorm(var3, 'Linf')

tvar1 = [0,1,2,3,4,5,6,7]
tvar2 = [[0,1,2,3],[4,5,6,7]]
tvar3 = [[[0,1],[2,3]],[[4,5],[6,7]]]

print 'L1E 1D = ', get_lerror(var1, tvar1, 'L1')
print 'L1E 2D = ', get_lerror(var2, tvar2, 'L1')
print 'L1E 3D = ', get_lerror(var3, tvar3, 'L1')
print 'L2E 1D = ', get_lerror(var1, tvar1, 'L2')
print 'L2E 2D = ', get_lerror(var2, tvar2, 'L2')
print 'L2E 3D = ', get_lerror(var3, tvar3, 'L2')
print 'LiE 1D = ', get_lerror(var1, tvar1, 'Linf')
print 'LiE 2D = ', get_lerror(var2, tvar2, 'Linf')
print 'LiE 3D = ', get_lerror(var3, tvar3, 'Linf')

a = 0

for i in range(8):
  a = a + var1[i]**2.0

a = a**0.5
print a
