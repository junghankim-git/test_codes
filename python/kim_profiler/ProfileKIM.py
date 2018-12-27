#!/usr/bin/env python
#-------------------------------------------------------------------------------
#
#> @brief
#>  - 
#>
#> @date 13MAR2015
#>  - JH KIM(jh.kim@kiaps.org) : first 
#>
#> @date 13MAR2015
#>  - JH KIM(jh.kim@kiaps.org) : Make the structures
#>
#> @date 03APR2015
#>  - JH KIM(jh.kim@kiaps.org) : Initailize with input path
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import cm
import os
import sys
import random
from optparse import OptionParser
#SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
#sys.path.append(SHARE_DIR)
from LibProfileKIM import *


opt = OptionParser()
#
## action: 'store', 'store_const', 'append', 'count', 'callback'
opt.add_option('-v', '--vertical',dest='vertical',  default='Eulerian', action='store',  help='Dynamics (vertically Eulerian, Lagrangian)')
opt.add_option('-i', '--input',   dest='isInput',   default=False, action='store_true',  help='Input')
opt.add_option('-d', '--diff',    dest='isDiff',    default=False, action='store_true',  help='Difference')
opt.add_option('-k', '--keys',    dest='isKeys',    default=False, action='store_true',  help='Print Keys')
opt.add_option('-s', '--save',    dest='isSave',    default=False, action='store_true',  help='Save Statistics')
opt.add_option('-p', '--plot',    dest='isPlot',    default=True,  action='store_false', help='Plotting')
opt.add_option('-w', '--wiki',    dest='isWiki',    default=False, action='store_true',  help='Print Wiki version')
opt.add_option('-r', '--read',    dest='isRead',    default=False, action='store_true',  help='Read Statistics')
#
(options, args) = opt.parse_args()


def CheckDiff(exps1, exps2):
   nexps1 = len(exps1)
   nexps2 = len(exps2)
   nmatch = 0
   for i in range(nexps2):
      for j in range(nexps1):
         if exps2[i] == exps1[j]:
            nmatch = nmatch + 1
            break
   if nexps2 != nmatch:
      print 'Check number of diff. exps....'
      print '  - exps1 = ', exps1
      print '  - exps2 = ', exps2
      quit()


if options.vertical == 'Eulerian':
   onLagrangian = False
   vert    = 'Eulerian'
elif options.vertical == 'Lagrangian':
   onLagrangian = True
   vert    = 'Lagrangian'
else:
   print 'check vertical... (-v option)'
   quit()



if not onLagrangian:
   # GAON2
   nprocs  = [5, 16, 32, 64, 80, 160, 320, 400, 640, 800, 960]
   baseDIR = '/home/jhkim/TestBed/KIM/0.21.micros/0.21.01/04.perf/exps'
   nexps   = len(nprocs)
   exps    = ['' for i in range(nexps)]
   for i in range(nexps):
      exps[i] = baseDIR+'/'+'%05d'%(nprocs[i])+'/profile'

   tag = 'gaon2'
   
   
   if options.isDiff:
       # - diff data
      diff_vert    = 'Eulerian'
      diff_nprocs  = [ 400 ]
      diff_baseDIR = '/home/jhkim/TestBed/KIM//0.21.micros/0.21.01/04.perf/checkMVAPICH2'
      diff_nexps   = len(diff_nprocs)
      diff_exps    = ['' for i in range(diff_nexps)]
      #diff_exps    = [diff_baseDIR+'/profile']
      for i in range(diff_nexps):
         diff_exps[i] = diff_baseDIR+'/'+'%05d'%(diff_nprocs[i])+'/profile'
   
      CheckDiff(nprocs, diff_nprocs)

else:

   # Haedam
   nprocs  = [1080, 2400, 3000, 3360, 4008, 4680, 5400, 6000, 6600, 7200]
   baseDIR = '/scratch/yslee/0.21.02/profiling'
   nexps   = len(nprocs)
   exps    = ['' for i in range(nexps)]
   for i in range(nexps):
      exps[i] = baseDIR+'/'+'%04d'%(nprocs[i])

   tag = 'haedam'


# Check arguments
isManual = False
needArgs = 0
if options.isInput: needArgs = needArgs + 1
if options.isDiff: needArgs = needArgs + 1

if len(args) < needArgs or len(args) == 0:
   isManual = True



# Make a profiler of KIM
prof = ProfileKIM(vert)

if options.isKeys:
   prof.PrintKeysInfo()

# Initialize of a profiler
if options.isRead:
   prof.Initialize_DB()
else:
   if isManual:
      prof.Initialize(nprocs, exps)
   else:
      prof.Initialize_DIR(args[0])
   if options.isSave:
      prof.WriteFile()

if not options.isDiff:
   if options.isWiki:
      prof.PrintWiki()
   if options.isPlot:
      prof.DoPlot(tag)
else:
   if isManual:
      prof.Initialize(diff_nprocs, diff_exps, isDiff=True)
   else:
      prof.Initialize_DIR(args[1], isDiff=True)
   if options.isWiki:
      prof.PrintWiki_Diff()
   if options.isPlot:
      prof.DoPlot_Diff()



print ' '
print ' '
print '------------------------------------------------------------------------------'
print '# Contact: Junghan Kim (jh.kim@kiaps.org)'
print '------------------------------------------------------------------------------'
