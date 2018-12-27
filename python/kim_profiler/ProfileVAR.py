#!/usr/bin/env python
#-------------------------------------------------------------------------------
#
#> @brief
#>  - 
#>
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import cm
import os
import sys
import random
from optparse import OptionParser
#sys.path.append('/home/jhkim/Study/Library/Shared/Python/ProfileKIM')
from LibProfileVAR import *


opt = OptionParser()
#
## action: 'store', 'store_const', 'append', 'count', 'callback'
#opt.add_option('-v', '--vertical',dest='vertical',  default='Eulerian', action='store',  help='Dynamics (vertically Eulerian, Lagrangian)')
#opt.add_option('-i', '--input',   dest='isInput',   default=False, action='store_true',  help='Input')
#opt.add_option('-d', '--diff',    dest='isDiff',    default=False, action='store_true',  help='Difference')
#opt.add_option('-k', '--keys',    dest='isKeys',    default=False, action='store_true',  help='Print Keys')
#opt.add_option('-s', '--save',    dest='isSave',    default=False, action='store_true',  help='Save Statistics')
#opt.add_option('-p', '--plot',    dest='isPlot',    default=True,  action='store_false', help='Plotting')
#opt.add_option('-w', '--wiki',    dest='isWiki',    default=False, action='store_true',  help='Print Wiki version')
#opt.add_option('-r', '--read',    dest='isRead',    default=False, action='store_true',  help='Read Statistics')
#
(options, args) = opt.parse_args()


nprocs = [320] 
exps   = ['./DataVAR/000320']
# Make a profiler of KIM
prof = ProfileVAR()

#prof.PrintKeysInfo()
prof.Initialize(nprocs, exps)
#prof.PrintWiki()
prof.DoPrediction()
prof.DoPlot()



print ' '
print ' '
print '------------------------------------------------------------------------------'
print '# Contact: Junghan Kim (jh.kim@kiaps.org)'
print '------------------------------------------------------------------------------'
