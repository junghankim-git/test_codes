#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt



print '# 1. start'
print ' '
#====================================================================
# Variables
#====================================================================
# experiment(0, ne120np4, T323) and expectation(1, ne120, T323)
wns      = [323+1, 323+1]
nes      = [120,   120]

#=== comm .vs. mm
#np       = 4
#nUPs     = [6*((np-1)**2)*(nes[i]**2)+2 for i in range(ndata)]
#size_mat = [wns[i]**2*nUPs[i] for i in range(ndata)]
#size_com = [wns[i]**2 for i in range(ndata)]
#frac_mat = float(size_mat[1])/float(size_mat[0])
#frac_com = float(size_com[1])/float(size_com[0])
#print '   - nups     = ', nUPs
#print '   - frac_mat = ', frac_mat
#print '   - frac_com = ', frac_com
#print '   - wns      = ', float(wns[1])/float(wns[0])
#print '   - nUPs     = ', float(nUPs[1])/float(nUPs[0])
#=== comm .vs. mm

# walltimes (calculation, communication, total)
nwall    = 3

max_nprocs = 0
nps        = 0
nprocs     = []
walltime   = [[] for i in range(nwall)]
speedup    = [[] for i in range(nwall)]
basewall   = [0.0 for i in range(nwall)]

nps_exp        = 0
nprocs_exp     = []
walltime_exp   = [[] for i in range(nwall)]
speedup_exp    = [[] for i in range(nwall)]
basewall_exp   = [0.0 for i in range(nwall)]



print ' '
print '# 2. read walltime and expectation'
#====================================================================
# Read data file (index: nprocs = 0, 6, 7, 11)
#====================================================================
ip = 0
im = 4
ic = 5
it = 6
filename = './walltime.dat'
infile   = open(filename, 'r')

il       = 0
for line in infile:
  if il != 0:
    sline = line.split()
    nprocs.append(int(sline[ip]))
    walltime[0].append(float(sline[im]))  # caculation
    walltime[1].append(float(sline[ic]))  # communication
    walltime[2].append(float(sline[it]))  # total
  il = il + 1
infile.close()
nps = len(nprocs)

max_nprocs = nprocs[-1]
max_nprocs_exp = max_nprocs*2
#max_nprocs_exp = max_nprocs*4
#max_nprocs_exp = max_nprocs*8

print 'max nproc = ', max_nprocs_exp


print ' '
print '# 3. expectaion'
#====================================================================
# Expectation (nprocs: max_nprocs_exp)
#====================================================================
for ip in range(nps):
  nprocs_exp.append(nprocs[ip])
  walltime_exp[0].append(walltime[0][ip])
  walltime_exp[1].append(walltime[1][ip])
  walltime_exp[2].append(walltime[2][ip])

ii = nps-1
numprocs = nprocs_exp[-1]
while True:
  numprocs = numprocs + 432
  if numprocs > max_nprocs_exp: break

  nprocs_exp.append(numprocs)
  walltime_exp[0].append(walltime_exp[0][-1]*float(numprocs-432)/float(numprocs))
  walltime_exp[1].append(walltime_exp[1][-1])
  walltime_exp[2].append(walltime_exp[0][-1]+walltime_exp[1][-1])

  nps_exp = len(nprocs_exp)

print '  - comm:', walltime[1][:]
print '  - experiment'
print '    walltime = ', walltime[2][-1]
print '  - expecation'
print '    walltime = ', walltime_exp[2][-1]

print ' '
print '# 4. speedup'
#====================================================================
# Speed-up
#====================================================================
for iw in range(nwall):
  basewall[iw] = walltime[iw][0]*nprocs[0]
  for ip in range(nps):
    speedup[iw].append(basewall[iw]/walltime[iw][ip])

for iw in range(nwall):
  basewall_exp[iw] = walltime_exp[iw][0]*nprocs_exp[0]
  for ip in range(nps_exp):
    speedup_exp[iw].append(basewall_exp[iw]/walltime_exp[iw][ip])

def speedup_func(n, p):
        w=1.0/(p+(1.0+p)/n)
        return w

speedup_fit  = [0.0 for i in range(nps_exp)]
popt, pcov = curve_fit(speedup_func, nprocs, speedup[-1])
speedup_fit  = speedup_func(nprocs_exp, popt)


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
pltfig = plt.figure(figsize=(14,10))
#pltfig.subplots_adjust(left=margin, bottom=2.0*margin, right=1.0-margin,top=1.0-margin-0.06, wspace=0.12, hspace=0.1)
pltfig.subplots_adjust(left=margin, bottom=1.5*margin, right=1.0-margin,top=1.0-margin-0.02, wspace=0.13, hspace=0.18)
pltfig.suptitle('Analysis 3DVar', size='18')

# default draw setting
colors  = ['k', 'k', 'k']
formats = ['o--', 'o-.', 'o-']
labels  = ['Matrix calculation', 'Communication', 'Total']
xmin   = 0
xmax   = nprocs[-1]*1.05

# plot 1-1
ymin   = 0
ymax   = max(walltime[2])*1.05
xlabel = 'Number of CPUs [#]'
ylabel = 'Wall-clock time [s]'
title  = 'Wall-clock time (experiment: ne120np4, T323, max. proc='+str(max_nprocs)+')'
axis11 = IniAxis(pltfig, xmin, xmax, ymin, ymax, xlabel, ylabel, title, 2, 2, 1)
for iw in range(nwall):
    axis11.plot(nprocs, walltime[iw], colors[iw]+formats[iw], label=labels[iw])
axis11.grid(True)
axis11.legend(loc=1,shadow=True)
#text
axis11.arrow(nprocs[-1], walltime[2][-1], 0.0, 1.5, head_width=150.0, head_length=0.7, fc='k', ec='k')
axis11.text(nprocs[-2], walltime[2][-1]+3.0, 'walltime = '+str(walltime[2][-1])+'s', color='k', size=15.0, ha='center', va='center', style='italic')

# plot 1-2
ymin   = 0
ymax   = xmax
xlabel = 'Number of CPUs [#]'
ylabel = 'Speed-up'
title  = 'Speed-up (experiment: ne120np4, T323, max. proc='+str(max_nprocs)+')'
axis12 = IniAxis(pltfig, xmin, xmax, ymin, ymax, xlabel, ylabel, title, 2, 2, 2)
axis12.plot(nprocs, nprocs, 'r-', label='Ideal')
for iw in range(nwall):
    axis12.plot(nprocs, speedup[iw], colors[iw]+formats[iw], label=labels[iw])
axis12.grid(True)
axis12.legend(loc=2,shadow=True)



xmin   = 0
xmax   = nprocs_exp[-1]*1.05
# plot 2-1
ymin   = 0
ymax   = max(walltime_exp[2])*1.05
xlabel = 'Number of CPUs [#]'
ylabel = 'Wall-clock time [s]'
title  = 'Wall-clock time (expectation: ne120np4, T323, max. proc='+str(max_nprocs_exp)+')'
axis21 = IniAxis(pltfig, xmin, xmax, ymin, ymax, xlabel, ylabel, title, 2, 2, 3)
for iw in range(nwall):
    axis21.plot(nprocs_exp, walltime_exp[iw], colors[iw]+formats[iw], label=labels[iw]+' (expectation)')
#for iw in range(nwall):
#    axis21.plot(nprocs, walltime[iw], colors[iw]+formats[iw], label=labels[iw]+'(experiment)')
axis21.grid(True)
axis21.legend(loc=1,shadow=True)
#text
axis21.arrow(nprocs[-1], walltime[2][-1], 0.0, 1.5, head_width=150.0, head_length=0.7, fc='k', ec='k')
axis21.text(nprocs[-2], walltime[2][-1]+3.0, 'walltime = '+str(walltime[2][-1])+'s', color='k', size=15.0, ha='center', va='center', style='italic')
axis21.arrow(nprocs_exp[-1], walltime_exp[2][-1], 0.0, 1.5, head_width=150.0, head_length=0.7, fc='k', ec='k')
axis21.text(nprocs_exp[-2], walltime_exp[2][-1]+3.0, 'walltime = '+str(walltime_exp[2][-1])+'s', color='k', size=15.0, ha='center', va='center', style='italic')


# plot 2-2
ymin   = 0
ymax   = xmax
xlabel = 'Number of CPUs [#]'
ylabel = 'Speed-up'
title  = 'Speed-up (experiment & expectation: ne120np4, T323, max. proc='+str(max_nprocs_exp)+')'
axis22 = IniAxis(pltfig, xmin, xmax, ymin, ymax, xlabel, ylabel, title, 2, 2, 4)
axis22.plot(nprocs_exp, nprocs_exp, 'r-', label='Ideal')
#for iw in range(nwall):
#    axis22.plot(nprocs_exp, speedup_exp[iw], colors[iw]+formats[iw], label=labels[iw]+'(expectation)')
#for iw in range(nwall):
#    axis22.plot(nprocs, speedup[iw], colors[iw]+formats[iw], label=labels[iw]+'(experiment)')
axis22.plot(nprocs, speedup[iw], 'ko-', label='Speed-up'+' (experiment)')
axis22.plot(nprocs_exp, speedup_exp[iw], 'k-.', label='Speed-up'+' (expectation)')
#axis22.plot(nprocs_exp, speedup_fit, 'r-.', label='Speed-up'+' (expectation)')
axis22.grid(True)
axis22.legend(loc=2,shadow=True)



pltfig.savefig('perf_anal.png')


plt.show()

