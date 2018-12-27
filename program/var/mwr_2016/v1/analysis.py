#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import scipy as sp
import netCDF4 as nc
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_toolkits.basemap import Basemap
#from matplotlib import cm



print '# 1. start'
print ' '
#====================================================================
# Variables
#====================================================================
# experiment(0, ne60np4, T179) and expectation(1, ne120, T320)
ndata    = 2
wns      = [179+1, 320+1]
nes      = [60,    120]
np       = 4
nUPs     = [6*((np-1)**2)*(nes[i]**2)+2 for i in range(ndata)]
size_mat = [wns[i]**2*nUPs[i] for i in range(ndata)]
size_com = [wns[i]**2 for i in range(ndata)]
frac_mat = float(size_mat[1])/float(size_mat[0])
frac_com = float(size_com[1])/float(size_com[0])
print '   - nups     = ', nUPs
print '   - frac_mat = ', frac_mat
print '   - frac_com = ', frac_com
print '   - wns      = ', float(wns[1])/float(wns[0])
print '   - nUPs     = ', float(nUPs[1])/float(nUPs[0])
# walltimes (calculation, communication, others, total)
nwall    = 4

nps      = 0
nprocs   = []
walltime = [[[] for i in range(nwall)] for j in range(ndata)]
speedup  = [[[] for i in range(nwall)] for j in range(ndata)]
basewall = [[0.0 for i in range(nwall)] for j in range(ndata)]



print ' '
print '# 2. read walltime and expectation'
#====================================================================
# Read data file (index: nprocs = 0, 6, 7, 11)
#====================================================================
ip = 0
im = 6
ic = 7
it = 11
filename = './walltime.dat'
infile   = open(filename, 'r')

il       = 0
for line in infile:
  if il != 0:
    sline = line.split()
    nprocs.append(int(sline[ip]))
    walltime[0][0].append(float(sline[im]))  # caculation
    walltime[0][1].append(float(sline[ic]))  # communication
    walltime[0][2].append(float(sline[it]))  # total
    others = walltime[0][2][il-1] - walltime[0][1][il-1] - walltime[0][0][il-1]
    walltime[0][3].append(others)    # others
    print nprocs[il-1], walltime[0][0][il-1], walltime[0][1][il-1], walltime[0][2][il-1], walltime[0][3][il-1]
  il = il + 1
infile.close()
nps = len(nprocs)



print ' '
print '# 3. expectaion'
#====================================================================
# Expectation (ne120, T320)
#====================================================================
for ip in range(nps):
  walltime[1][0].append(walltime[0][0][ip]*frac_mat)
  walltime[1][1].append(walltime[0][1][ip]*frac_com)
  walltime[1][3].append(walltime[0][3][ip])
  walltime[1][2].append(walltime[1][0][ip] + walltime[1][1][ip] + walltime[1][3][ip])
  print nprocs[ip], walltime[1][0][ip], walltime[1][1][ip], walltime[1][2][ip], walltime[1][3][ip]



print ' '
print '# 4. speedup'
#====================================================================
# Speed-up
#====================================================================
for id in range(ndata):
  for iw in range(nwall):
    basewall[id][iw] = walltime[id][iw][0]*nprocs[0]
    for ip in range(nps):
      speedup[id][iw].append(basewall[id][iw]/walltime[id][iw][ip])






print ' '
print '# 5. draw'
#====================================================================
# Draw
#====================================================================

def IniAxis(pltfig_in, xmin_in, xmax_in, ymin_in, ymax_in, xlabel_in, ylabel_in, title_in, nposx_in, nposy_in, pos_in):
   axis_out = pltfig_in.add_subplot(nposx_in,nposy_in,pos_in)
   axis_out.set_xlim(xmin_in,xmax_in)
   axis_out.set_ylim(ymin_in,ymax_in)
   axis_out.set_title(title_in)
   axis_out.set_xlabel(xlabel_in)
   axis_out.set_ylabel(ylabel_in)
   return axis_out

# canvas
margin = 0.04
#pltfig = plt.figure(figsize=(18,8))
pltfig = plt.figure(figsize=(18,16))
#pltfig.subplots_adjust(left=margin, bottom=2.0*margin, right=1.0-margin,top=1.0-margin-0.06, wspace=0.12, hspace=0.1)
pltfig.subplots_adjust(left=margin, bottom=1.5*margin, right=1.0-margin,top=1.0-margin-0.02, wspace=0.13, hspace=0.18)
pltfig.suptitle('Analysis 3DVar', size='18')

# default draw setting
colors  = ['k', 'k', 'k']
formats = ['o--', 'o-.', 'o-']
labels  = ['Matrix calculation', 'Communication', 'Total']
xmin   = 0
xmax   = nprocs[-1]*1.05

idx = 0
# plot 1-1
ymin   = 0
ymax   = max(walltime[idx][2])*1.05
xlabel = 'Number of CPUs [#]'
ylabel = 'Wall-clock time [s]'
title  = 'Wall-clock time (experiment: ne60np4, T179)'
axis11 = IniAxis(pltfig, xmin, xmax, ymin, ymax, xlabel, ylabel, title, 2, 2, 1)
for iw in range(nwall-1):
    axis11.plot(nprocs, walltime[idx][iw], colors[iw]+formats[iw], label=labels[iw])
axis11.grid(True)
axis11.legend(loc=1,shadow=True)

# plot 1-2
ymin   = 0
#ymax   = max(speedup[0][0])*1.05
ymax   = xmax
xlabel = 'Number of CPUs [#]'
ylabel = 'Speed-up'
title  = 'Speed-up (experiment: ne60np4, T179)'
axis12 = IniAxis(pltfig, xmin, xmax, ymin, ymax, xlabel, ylabel, title, 2, 2, 2)
axis12.plot(nprocs, nprocs, 'r-', label='ideal')
for iw in range(nwall-1):
    axis12.plot(nprocs, speedup[idx][iw], colors[iw]+formats[iw], label=labels[iw])
axis12.grid(True)
axis12.legend(loc=2,shadow=True)

idx = 1
# plot 2-1
ymin   = 0
ymax   = max(walltime[idx][2])*1.05
xlabel = 'Number of CPUs [#]'
ylabel = 'Wall-clock time [s]'
title  = 'Wall-clock time (expectation: ne120np4, T320)'
axis21 = IniAxis(pltfig, xmin, xmax, ymin, ymax, xlabel, ylabel, title, 2, 2, 3)
for iw in range(nwall-1):
    axis21.plot(nprocs, walltime[idx][iw], colors[iw]+formats[iw], label=labels[iw])
axis21.grid(True)
axis21.legend(loc=1,shadow=True)

# plot 2-2
ymin   = 0
#ymax   = max(speedup[0][0])*1.05
ymax   = xmax
xlabel = 'Number of CPUs [#]'
ylabel = 'Speed-up'
title  = 'Speed-up (expectation: ne120np4, T320)'
axis22 = IniAxis(pltfig, xmin, xmax, ymin, ymax, xlabel, ylabel, title, 2, 2, 4)
axis22.plot(nprocs, nprocs, 'r-', label='ideal')
for iw in range(nwall-1):
    axis22.plot(nprocs, speedup[idx][iw], colors[iw]+formats[iw], label=labels[iw])
axis22.grid(True)
axis22.legend(loc=2,shadow=True)




pltfig.savefig('perf_anal.png')


plt.show()

