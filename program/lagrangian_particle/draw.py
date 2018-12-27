#!/usr/bin/env python

import os
import sys
import random
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


# setting plot
nplot = 50
# view on-off
onoffview = True

def load_file():
   infile = nc.Dataset('result.nc',mode='r')
   #
   ivis   = infile.ivis
   uptype = infile.uptype
   bndry  = infile.bndry
   nstep  = infile.nstep
   dt     = infile.dt
   #
   ndims  = len(infile.dimensions['ndims'])
   npart  = len(infile.dimensions['npart'])
   time   = len(infile.dimensions['time'])
   #
   xmins = infile.variables['xmin'][:].tolist()
   xmaxs = infile.variables['xmax'][:].tolist()
   xtmp  = infile.variables['x'][:].tolist()
   vtmp  = infile.variables['v'][:].tolist()
   ax    = infile.variables['ax'][:].tolist()
   av    = infile.variables['av'][:].tolist()
   infile.close()

   x = [[[0.0 for i in range(npart)] for j in range(ndims)] for k in range(time)]
   v = [[[0.0 for i in range(npart)] for j in range(ndims)] for k in range(time)]
   for it in range(time):
      for ip in range(npart):
         for id in range(ndims):
            x[it][id][ip] = xtmp[it][ip][id]
            v[it][id][ip] = vtmp[it][ip][id]

   ex = [0.0 for i in range(time)]
   ev = [0.0 for i in range(time)]
   for it in range(time):
     dx ,dy  = ax[0][it]-x[it][0][0],ax[1][it]-x[it][1][0]
     dvx,dvy = av[0][it]-v[it][0][0],av[1][it]-v[it][1][0]
     ex[it]  = np.sqrt(dx*dx+dy*dy)
     ev[it]  = np.sqrt(dvx*dvx+dvy*dvy)

   return ivis,uptype,bndry,nstep,dt,ndims,npart,time,xmins,xmaxs,ax,av,x,v,ex,ev
   

ivis,uptype,bndry,nstep,dt,ndims,npart,time,xmins,xmaxs,ax,av,x,v,ex,ev = load_file()
print 'ivis    = {}'.format(ivis)
print 'uptype  = {}'.format(uptype)
print 'bndry   = {}'.format(bndry)
print 'nstep   = {}'.format(nstep)
print 'dt      = {}'.format(dt)
print 'ndims   = {}'.format(ndims)
print 'npart   = {}'.format(npart)
print 'time    = {}'.format(time)
x0 = [0.0 for i in range(ndims)]
v0 = [0.0 for i in range(ndims)]
for i in range(ndims):
   x0[i] = x[0][i][0]
   v0[i] = v[0][i][0]

# make directory for figures
namedir = './testcase{0:}'.format(npart)
command = 'rm -rf '+namedir
os.system(command)
command = 'mkdir '+namedir
os.system(command)
xmin = xmins[0]
xmax = xmaxs[0]
ymin = xmins[1]
ymax = xmaxs[1]


# plot
px = 0.1*(xmax-xmin)
py = 0.08*(ymax-ymin)
pltfig = plt.figure(figsize=(px,py))
ma = 0.07
pltfig.subplots_adjust(left=ma,right=1.0-ma,bottom=ma,top=1.0-ma,wspace=0.0,hspace=0.0)
axis = pltfig.add_subplot(1,1,1)
axis.clear()

it = 0
title = 't = {:5.3f}sec, x = {:6.2f}m, y = {:6.2f}m'.format(dt*it,x[it][0][0],x[it][1][0])
axis.set_title(title)
axis.set_xlabel('x ['+r'$m$'+']')
axis.set_ylabel('y ['+r'$m$'+']')
axis.set_xlim(xmin,xmax)
axis.set_ylim(ymin,ymax)
plt.ion()

if npart==1:
   axis.axhline(y=x[0][1][0],xmin=xmin,xmax=xmax)
   ana, = axis.plot(ax[0],ax[1],'k--',label='analytic') 
   pana, = axis.plot(ax[0][0],ax[1][0],'ko',label='analytic')
pexp, = axis.plot(x[0][0][:],x[0][1][:],'ro',label='numerical')

arrsize = 100.0/(xmax-xmin)
arrmag  = 40.0/(xmax-xmin)
arr = []
for ip in range(npart):
   arr.append(axis.arrow(x[0][0][ip],x[0][1][ip],arrmag*v[0][0][ip],arrmag*v[0][1][ip],head_width=arrsize,head_length=arrsize,fc='k',ec='k'))

if onoffview: plt.show()

# time loop
for it in range(time):
   if it%nplot == 0:
      print '# step = {0:00005d}'.format(it)
      if npart==1:
         title = 't = {:5.3f}sec, x = {:6.2f}m, y = {:6.2f}m, ex = {:6.2f}m, ev = {:6.2f}m/s'.format(dt*it,x[it][0][0],x[it][1][0],ex[it],ev[it])
         pana.set_data(ax[0][it],ax[1][it])
      else:
         title = 't = {:5.3f}sec, x = {:6.2f}m, y = {:6.2f}m'.format(dt*it,x[it][0][0],x[it][1][0])
      axis.set_title(title)
      pexp.set_data(x[it][0][:],x[it][1][:])
      for ip in range(npart):
         arr[ip].remove()
         arr[ip] = axis.arrow(x[it][0][ip],x[it][1][ip],arrmag*v[it][0][ip],arrmag*v[it][1][ip],head_width=arrsize,head_length=arrsize,fc='k',ec='k')
      filename = '{0:}/{1:000005d}.png'.format(namedir,int(it/nplot))
      pltfig.savefig(filename)
      if onoffview:
         if it == time-1:
            plt.show(True)
         else:
            plt.draw()




