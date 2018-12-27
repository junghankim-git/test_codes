#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


datafile = './walltime.dat'


def nlen_file(filename):
    lines = 0
    for line in open(filename):
        lines = lines + 1
    return lines

print '# 1. start'
print ' '
#====================================================================
# Variables
#====================================================================

# walltimes (total, comm, io)
nps        = nlen_file(datafile)
nws        = 3

nprocs     = []
walltime   = [[] for i in range(nws)]
speedup    = [[] for i in range(nws)]



print ' '
print '# 2. read walltime and expectation'
#====================================================================
# Read data file 
#====================================================================
infile   = open(datafile, 'r')

for line in infile:
    sline = line.split()
    nprocs.append(int(sline[0]))
    walltime[2].append(float(sline[1]))                  # i/o
    walltime[0].append(float(sline[2]))                  # total
    walltime[1].append(float(sline[2])-float(sline[1]))  # computation
infile.close()
nps = len(nprocs)

for i in range(nps):
    print nprocs[i], walltime[0][i], walltime[1][i], walltime[2][i]

nprocs_u = 7680
total_u  = 1306.
model_u  = 1012.
io_u     =  293.

print ' '
print '# 4. speedup'
#====================================================================
# Speed-up
#====================================================================
basewall = [0.0 for i in range(nws)]
for iw in range(nws):
    basewall[iw] = walltime[iw][0]*nprocs[0]
    for ip in range(nps):
        speedup[iw].append(basewall[iw]/walltime[iw][ip])

def speedup_func(n, p):
    w=1.0/(p+(1.0+p)/n)
    return w

speedup_id = [0.0 for i in range(nps)]
for i in range(nps):
    speedup_id[i] = nprocs[i]



print ' '
print '# 5. draw'
#====================================================================
# Draw
#====================================================================

def IniAxis(pltfig_in,xmin_in,xmax_in,ymin_in,ymax_in,xlabel_in,ylabel_in,title_in,nposx_in,nposy_in,pos_in):
   fsize = 11
   axis_out = pltfig_in.add_subplot(nposx_in,nposy_in,pos_in)
   axis_out.set_xlim(xmin_in,xmax_in)
   axis_out.set_ylim(ymin_in,ymax_in)
   axis_out.set_title(title_in)
   axis_out.set_xlabel(xlabel_in,fontsize=fsize)
   axis_out.set_ylabel(ylabel_in,fontsize=fsize)
   axis_out.tick_params(axis='both',labelsize=fsize)
   #axis_out.set_xticklabels(labelsize=10)
   #axis_out.set_yticklabels(labelsize=10)
   return axis_out

# canvas
margin = 0.08
#pltfig = plt.figure(figsize=(18,8))
pltfig = plt.figure(figsize=(7.5,5))
pltfig.subplots_adjust(left=margin,bottom=margin+0.02,right=1.0-margin,top=1.0-margin-0.02,wspace=0.0,hspace=0.0)
#pltfig.suptitle('KIM: ne240(12km)',size='18')

# default draw setting
colors  = ['k', 'r', 'b']
formats = ['o-', 'o-', 'o-']
labels  = ['Total', 'Computation', 'I/O']
xmin    = 0
xmax    = nprocs[-1]*1.05

# plot 1
ymin   = 0
ymax   = max(walltime[0])*1.05
xlabel = 'number of CPUs [#]'
ylabel = 'wall-clock time [s]'
title  = 'KIM scalability (v2.5)'
axis1 = IniAxis(pltfig,xmin,xmax,ymin,ymax,xlabel,ylabel,title,1,1,1)
for iw in range(nws):
    axis1.plot(nprocs,walltime[iw],colors[iw]+formats[iw],label=labels[iw])
y_labels = axis1.get_yticklabels()
for label in y_labels:
    label.set_rotation(90)
axis1.grid(True)
axis1.legend(loc=1,fontsize=11,shadow=True)

pltfig.savefig('perf_anal_1.png')

axis1.plot([nprocs_u],[total_u],'k^',label=labels[0]+' (update)')
axis1.plot([nprocs_u],[model_u],'r^',label=labels[1]+' (update)')
axis1.plot([nprocs_u],[io_u],'b^',label=labels[2]+' (update)')
axis1.arrow(nprocs_u,walltime[0][2],0,(total_u-walltime[0][2])*0.8,head_width=350.0,head_length=70.,fc='m',ec='m')
axis1.arrow(nprocs_u,walltime[1][2],0,(model_u-walltime[1][2])*0.8,head_width=350.0,head_length=70.,fc='m',ec='m')
axis1.arrow(nprocs_u,walltime[2][2],0,(io_u-walltime[2][2])*0.8,head_width=350.0,head_length=70.,fc='m',ec='m')
axis1.legend(loc=1,fontsize=11,shadow=True)

pltfig.savefig('perf_anal_2.png')

plt.show()

