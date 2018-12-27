#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


ntasks = 100

procs  = [i+1 for i in range(ntasks+1)]

nprocs = len(procs)

wtime_task   = 1.0
wtime_all    = float(ntasks)*wtime_task
#wtime      = [[0.0 for i in range(nprocs)] for j in range(2)]
wtime_ideal0 = [0.0 for i in range(nprocs)]
wtime_ideal1 = [0.0 for i in range(nprocs)]
wtime_ideal2 = [0.0 for i in range(nprocs)]
speed_ideal0 = [0.0 for i in range(nprocs)]
speed_ideal1 = [0.0 for i in range(nprocs)]
speed_ideal2 = [0.0 for i in range(nprocs)]


procs_t  = [1, 11, 20, 21, 50, 51, 75, 101]
nprocs_t = len(procs_t)
#wtime_t  = [0.0 for i in range(nprocs_t)]
wtime_t  = [100.009, 10.002, 6.008, 5.003, 3.008, 2.003, 2.002, 1.008]
speed_t  = [0.0 for i in range(nprocs_t)]
for ip in range(nprocs_t):
    speed_t[ip] = wtime_all/wtime_t[ip]


for ip in range(nprocs):
    wtime_ideal0[ip] = wtime_all/float(procs[ip])

    wtime_ideal1[ip] = float(int(wtime_all)/procs[ip])
    if ntasks%procs[ip]!=0:
        wtime_ideal1[ip] = wtime_ideal1[ip]+wtime_task

    if procs[ip]==1:
        wtime_ideal2[ip] = wtime_all
    else:
        wtime_ideal2[ip] = float(int(wtime_all)/(procs[ip]-1))
        if ntasks%(procs[ip]-1)!=0:
            wtime_ideal2[ip] = wtime_ideal2[ip]+wtime_task

    speed_ideal0[ip] = wtime_all/wtime_ideal0[ip]
    speed_ideal1[ip] = wtime_all/wtime_ideal1[ip]
    speed_ideal2[ip] = wtime_all/wtime_ideal2[ip]


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
pltfig = plt.figure(figsize=(7.5,6))
pltfig.subplots_adjust(left=margin,bottom=margin+0.02,right=1.0-margin,top=1.0-margin-0.02,wspace=0.0,hspace=0.0)

xmin   = 0
xmax   = procs[-1]*1.05
ymin   = 0
ymax   = xmax
xlabel = 'number of cores [#]'
ylabel = 'speed-up'
title  = 'master-slave scalability, {} tasks'.format(ntasks)
axis1  = IniAxis(pltfig,xmin,xmax,ymin,ymax,xlabel,ylabel,title,1,1,1)
axis1.plot(procs,speed_ideal0,'k--',label='ideal')
axis1.plot(procs,speed_ideal1,'b--',label='ideal (integer)')
axis1.plot(procs,speed_ideal2,'r-' ,label='ideal (master-slave)')
axis1.plot(procs_t,speed_t,   'g*' ,label='experiment')
axis1.legend(loc=2,fontsize=11,shadow=True)
axis1.grid(True)


pltfig.savefig('perf.png')
plt.show()





