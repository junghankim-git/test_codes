#!/usr/bin/env python
from numpy import *
from scipy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time
from matplotlib import cm
import os
import sys



####  Read Meta data            ####
metafilename = 'info.dat'
metafile = open(metafilename, 'r')
testcase = int(metafile.readline())
npoints  = int(metafile.readline())
zmin     = float(metafile.readline())
zmax     = float(metafile.readline())
metafile.close()

nplot = 4
npart = 2


####  Initialization variables  ####
x = array([[0.0 for i in arange(npoints)] for j in range(nplot)])
datatmp = array([[[0.0]*npoints for i in arange(npart)] for j in arange(nplot)])


####  Read input file           ####
for i in arange(nplot):
  for j in arange(npoints):
    if i == 1:
      x[i][j] = float(j)
    else:
      x[i][j] = (zmax-zmin)*float(j)/npoints + zmin

filename = 'Result.dat'
infile = open(filename, 'r')

for i in arange(nplot):
  for j in arange(npoints):
    line = infile.readline()
    splitline = str(line).split()
    datatmp[i][0][j] = float(splitline[0])
    datatmp[i][1][j] = float(splitline[1])

infile.close()




# - scale
xmin   = array([0.0 for i in arange(nplot)])
xmax   = array([0.0 for i in arange(nplot)])
ymin   = array([[0.0]*npart for i in arange(nplot)])
ymax   = array([[0.0]*npart for i in arange(nplot)])
title  = array([[' '*35]*npart for i in arange(nplot)])
color  = ['r-', 'g-']


for i in range(nplot):
  if i == 1:
    xmin[i] = 0
    xmax[i] = npoints
  else:
    xmin[i] = zmin
    xmax[i] = zmax
  for j in range(npart):
    ymin[i][j] = min(datatmp[i][j])
    ymax[i][j] = max(datatmp[i][j])
    margin = (ymax[i][j]-ymin[i][j])*0.1
    ymin[i][j] = ymin[i][j]-margin
    ymax[i][j] = ymax[i][j]+margin
    if ymin[i][j] == ymax[i][j]:
      ymin[i][j] = -0.01
      ymax[i][j] =  0.01
    if i == 0:
      title[i][j] = 'Input ['
    elif i == 1:
      title[i][j] = 'Frequency ['
    elif i == 2:
      title[i][j] = 'Differential ['
    elif i == 3:
      title[i][j] = 'Integration ['
    if j == 0:
      title[i][j] = title[i][j] + 'real part]'
    elif j == 1:
      title[i][j] = title[i][j] + 'imaginary part]'



####  Plotting                  ####
pltfig   = plt.figure(figsize=(10,12))
pltfig.suptitle('FFT', fontsize=18)

pltgrid = gridspec.GridSpec(nplot, npart)
axisdic = dict()

for i in range(nplot):
  for j in range(npart):
    axis = plt.subplot(pltgrid[i,j])
    axis.set_title(title[i][j])
    axis.set_xlim(xmin[i],xmax[i])
    axis.set_ylim(ymin[i][j],ymax[i][j])
    axis.plot(x[i], datatmp[i][j], color[j])
    axisdic[(i,j)] = axis


#axis.set_xlabel(r'$M/M_H$')
#axis.set_ylabel(r'$\frac{Nz}{U}$')
#axis.set_xlim(xmin,xmax)
#axM.plot(x, input[0], '-', label=legend_labels[i])
#axM.legend(shadow=True, loc=2)



outfilename = 'FFT'+str(testcase)+'.png'
pltfig.savefig(outfilename)

plt.show(True)
