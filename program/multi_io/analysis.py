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

#==== ideal cases ====#
wtime_task   = 1.0
wtime_all    = float(ntasks)*wtime_task
wtime_ideal0 = [0.0 for i in range(nprocs)]
wtime_ideal1 = [0.0 for i in range(nprocs)]
wtime_ideal2 = [0.0 for i in range(nprocs)]
speed_ideal0 = [0.0 for i in range(nprocs)]
speed_ideal1 = [0.0 for i in range(nprocs)]
speed_ideal2 = [0.0 for i in range(nprocs)]


#==== test cases ====#
procs_t  = [1, 11, 20, 21, 50, 51, 75, 101]
nprocs_t = len(procs_t)
#wtime_t  = [0.0 for i in range(nprocs_t)]
wtime_t = [100.009, 10.002,  6.008,  5.003,  3.008,  2.003,  2.002,  1.008]
wtime_w = [ 25.906, 12.749,  7.412,  7.570,  5.689,  5.585,  5.182,  5.784]
wtime_r = [ 49.012, 12.549,  8.881, 10.502,  5.307,  5.262,  4.947,  4.665]
wtime_c = [ 68.941, 20.951, 12.197, 12.966,  7.791,  7.701,  7.335,  6.694]
speed_t = [0.0 for i in range(nprocs_t)]
speed_w = [0.0 for i in range(nprocs_t)]
speed_r = [0.0 for i in range(nprocs_t)]
speed_c = [0.0 for i in range(nprocs_t)]
wtime_w_all = wtime_w[0]
wtime_r_all = wtime_r[0]
wtime_c_all = wtime_c[0]
for ip in range(nprocs_t):
    speed_t[ip] = wtime_all/wtime_t[ip]
    speed_w[ip] = wtime_w_all/wtime_w[ip]
    speed_r[ip] = wtime_r_all/wtime_r[ip]
    speed_c[ip] = wtime_c_all/wtime_c[ip]


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


def axis_initialize(pltfig_in,xmin_in,xmax_in,ymin_in,ymax_in,xlabel_in,ylabel_in,title_in,nposx_in,nposy_in,pos_in):
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
axis1  = axis_initialize(pltfig,xmin,xmax,ymin,ymax,xlabel,ylabel,title,1,1,1)
axis1.plot(procs,speed_ideal0,'k--',label='ideal')
axis1.plot(procs,speed_ideal1,'b--',label='ideal (integer)')
axis1.plot(procs,speed_ideal2,'r-' ,label='ideal (master-slave)')
axis1.plot(procs_t,speed_t,   'g*' ,label='experiment (sleep)')
axis1.plot(procs_t,speed_w,   'bo' ,label='experiment (write)')
axis1.plot(procs_t,speed_r,   'ro' ,label='experiment (read)')
#axis1.plot(procs_t,speed_c,   'go' ,label='experiment (copy)')
axis1.legend(loc=2,fontsize=11,shadow=True)
axis1.grid(True)


pltfig.savefig('perf.png')
plt.show()





