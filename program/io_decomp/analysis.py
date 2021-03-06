#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


nplots = 2

print ' '
print '# 1. start'
#====================================================================
# Variables
#====================================================================
nios       = 4
nps        = 100
max_nprocs = 400
dproc      = max_nprocs/nps

max_iproc  = 20

nprocs = [dproc*i for i in range(nps)]
speed  = [[0.0 for j in range(nps)] for i in range(nios)]




print ' '
print '# 2. gen base function'

xs = [[0.0 for i in range(nps)] for j in range(nios)]
for i in range(nios):
    amp = 2.1-(0.05*i)
    for j in range(nps):
        nproc = float(nprocs[j])/float(i+1)
        if nproc<16:
            xx = float(nproc-max_iproc)/20
        else:
            xx = float(nproc-max_iproc)/80
        speed[i][j] = amp*np.exp(-0.5*(xx+np.exp(-xx)))




print ' '
print '# 3. experiment'
#====================================================================
# Expectation (nprocs: max_nprocs_exp)
#====================================================================

speed_016 = 0


print ' '
print '# 4. draw'
#====================================================================
# Draw
#====================================================================

if nplots==1:
  ncols = 1
  nrows = 1
else:
  ncols = 2
  nrows = 1

def axis_initialize(pltfig_in,xmin_in,xmax_in,ymin_in,ymax_in,xlabel_in,ylabel_in,title_in,nposx_in,nposy_in,pos_in):
   axis_out = pltfig_in.add_subplot(nposx_in,nposy_in,pos_in)
   axis_out.set_xlim(xmin_in,xmax_in)
   axis_out.set_ylim(ymin_in,ymax_in)
   axis_out.set_title(title_in)
   axis_out.set_xlabel(xlabel_in)
   axis_out.set_ylabel(ylabel_in)
   return axis_out

# canvas
margin = 0.08
pltfig = plt.figure(figsize=(8*ncols,6))
#pltfig.subplots_adjust(left=margin,bottom=2.0*margin,right=1.0-margin,top=1.0-margin-0.06,wspace=0.12,hspace=0.1)
pltfig.subplots_adjust(left=margin,bottom=1.5*margin,right=1.0-margin,top=1.0-margin-0.02,wspace=0.13,hspace=0.18)
pltfig.suptitle('I/O performance (write)',size='16')

# default draw setting
colors  = ['k','b','g','c']
#formats = ['-','--','-.',':']
formats = ['-','-','-','-']
labels  = ['nios=ncoms','nios=ncoms/2','nios=ncoms/3','nios=ncoms/4']

xmin   = 0
xmax   = nprocs[-1]*1.02
ymin   = 0
ymax   = max(speed[0])*1.3

# plot 1-1
xlabel = 'number of procs [#]'
ylabel = 'writing speed [GB/s]'
title  = 'expectaion'
axis = axis_initialize(pltfig,xmin,xmax,ymin,ymax,xlabel,ylabel,title,nrows,ncols,1)
for io in range(nios):
    axis.plot(nprocs,speed[io],colors[io]+formats[io],label=labels[io])
axis.grid(True)
axis.legend(loc=1,shadow=True)


if nplots>1:
    title = 'experiment'
    width = 16
    axis2 = axis_initialize(pltfig,xmin,xmax,ymin,ymax,xlabel,ylabel,title,nrows,ncols,2)
    #axis2.bar( 16-width/2.0,1.473,width=width/1.0,bottom=0.0,color='k',edgecolor='k',label='016p-016io')
    #axis2.bar(160-width/1.0,0.636,width=width/2.0,bottom=0.0,color='k',edgecolor='k',label='160p-160io')
    #axis2.bar(160-width/2.0,1.146,width=width/2.0,bottom=0.0,color='r',edgecolor='r',label='160p-016io')
    axis2.plot([16,160],[1.473,0.636],'ko',label='nios=ncoms')
    axis2.plot(   [160],      [1.146],'ro',label='nios=ncoms/10')
    axis2.grid(True)
    axis2.legend(loc=1,shadow=True)

pltfig.savefig('test.png')
print ' - save fig'
print ' '


plt.show()

