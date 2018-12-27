#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import cm
import os
import sys
import random
SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
sys.path.append(SHARE_DIR)
from LibBndryCondition import *
from LibTimeIntegration import *
from LibInteractionParticles import *


fgravity = -9.80665

def FreeFall(ndims_in, t_in, x0, v0):
   x1 = [0.0 for i in range(ndims_in)]
   v1 = [0.0 for i in range(ndims_in)]
   a  = [0.0 for i in range(ndims_in)]
   a[ndims_in-1] = fgravity

   for i in range(ndims_in):
      x1[i] = x0[i] + t_in*v0[i] + 0.5*a[i]*t_in*t_in
      v1[i] = v0[i] + t_in*a[i]
   return x1, v1

######################
###  Test Case
######################

#testcase = 1  # 1 particle with a gravity force
#testcase = 11 # 1 particle with a gravity and viscosity
testcase = 2  # two particles with repulsive force
#testcase = 3  # three particles with repulsive force
#testcase = 4  # 40 particles with repulsive force


######################
###  Make Directory for Figures
######################
namedir = './testcase'+str(testcase)
command = 'rm -rf '+namedir
os.system(command)
command = 'mkdir '+namedir
os.system(command)

######################
### Dimensions
######################
ndims   = 2

######################
### Viscosity
######################
lvis = False


######################
### Domain Size
######################
xmin   = -50.0
xmax   =  50.0
ymin   = -50.0
ymax   =  50.0

######################
### Initial State
######################
if testcase == 1:
   nparts  = 1
elif testcase == 2:
   nparts  = 2
elif testcase == 3:
   nparts  = 3
elif testcase == 4:
   nparts  = 40
elif testcase == 11:
   nparts  = 1
else:
   print 'Check1 [testcase]...'
   quit()

x0      = np.array([[0.0]*ndims for i in np.arange(nparts)])
v0      = np.array([[0.0]*ndims for i in np.arange(nparts)])
gravity = np.array([[0.0]*ndims for i in np.arange(nparts)])


if testcase == 1 or testcase == 11:
   vmag    = 10.0
   vang    = 0.0
   x0[0,:] = [-50.0, 30.0]
   v0[0,:] = [vmag*np.cos(vang), vmag*np.sin(vang)]
elif testcase == 2:
   x0[0,0] =-30.0
   x0[0,1] = -1.0
   x0[1,0] =+30.0
   x0[1,1] =  1.0

   vmag    = 15.0
   v0[0,0] = vmag
   v0[0,1] = 00.0
   v0[1,0] = -vmag
   v0[1,1] = 00.0
elif testcase == 3:
   vmag    = 15.0
   x0[0,0] =-30.0
   x0[0,1] =-11.0
   v0[0,0] = vmag
   v0[0,1] =  0.0

   x0[1,0] = 30.0
   x0[1,1] =-09.0
   v0[1,0] =-vmag
   v0[1,1] =  0.0

   x0[2,0] = 00.0
   x0[2,1] = 10.0
   v0[2,0] =  0.0
   v0[2,1] =-vmag
elif testcase == 4:
   vmag    = 20.0
   for ipart in range(nparts):
      xrand   = (xmax - xmin)*random.random() + xmin
      yrand   = (ymax - ymin)*random.random() + ymin
      vang    = 2.0*3.141592*random.random()
      x0[ipart,:] = [xrand, yrand]
      v0[ipart,:] = [vmag*np.cos(vang), vmag*np.sin(vang)]
else:
   print 'Check2 [testcase]...'
   quit()

######################
### Time Integration
######################
dt     = 0.001
nstep  = 10000
uptype = 'ForwardEuler'
#uptype = 'Leapfrog'  # not operated boundary condition
#uptype = 'RK2'       # not completed

if testcase == 2 or testcase == 3:
   nstep = 5000

######################
### Analytic Solution
######################
if testcase == 1 or testcase == 11:
   tana   = np.sqrt(2.0*(x0[0,1]-ymin)/-fgravity)
   xana   = xmin+v0[0,0]*tana
   print 'analytic solution: t = ', tana, ', x = ', xana

   # theory
   x_t     = np.array([[0.0 for j in range(nstep)] for i in np.arange(ndims)])

   for istep in range(nstep):
      t_t = dt*istep
      if t_t <= tana:
         [x_t[0][istep], x_t[1][istep]], tmp = FreeFall(ndims, t_t, x0[0], v0[0])
      else:
         [x_t[0][istep], x_t[1][istep]], tmp = FreeFall(ndims, tana, x0[0], v0[0])
         tr = t_t - tana
         xr = [x_t[0][istep], x_t[1][istep]]
         vr = [tmp[0], -tmp[1]]
         [x_t[0][istep], x_t[1][istep]], tmp = FreeFall(ndims, tr, xr, vr)



######################
### Objects: Time Integration, Boundary Condition, Solver
######################
tl     = TimeIntegration('Lagrangian',dt,nstep,uptype,ndims,nparts)
if testcase == 1:
   xbndry = BndryCondition('Lagrangian', 'Periodic', nparts, ndims, 0, xmin, xmax)
   ybndry = BndryCondition('Lagrangian', 'Reflection', nparts, ndims, 1, ymin, ymax)
   solver = InteractionParticles(nparts, ndims, 'Gravity', lvis=lvis)
elif testcase == 11:
   xbndry = BndryCondition('Lagrangian', 'Periodic', nparts, ndims, 0, xmin, xmax)
   ybndry = BndryCondition('Lagrangian', 'Reflection', nparts, ndims, 1, ymin, ymax)
   solver = InteractionParticles(nparts, ndims, 'Gravity', lvis=True)
else:
   xbndry = BndryCondition('Lagrangian', 'Reflection', nparts, ndims, 0, xmin, xmax)
   ybndry = BndryCondition('Lagrangian', 'Reflection', nparts, ndims, 1, ymin, ymax)
   solver = InteractionParticles(nparts, ndims, 'BetweenParticles',lvis=lvis)

######################
### Variables
######################
x      = np.array([[[0.0]*ndims for i in np.arange(nparts)] for k in np.arange(tl.ntl)])
v      = np.array([[[0.0]*ndims for i in np.arange(nparts)] for k in np.arange(tl.ntl)])
acc    = np.array([[[0.0]*ndims for i in np.arange(nparts)] for k in np.arange(tl.ntl)])
xRHS   = np.array([[0.0]*ndims for i in np.arange(nparts)])
vRHS   = np.array([[0.0]*ndims for i in np.arange(nparts)])

######################
### Set Initial Values
######################
x[tl.n0,:,0] = x0[:,0]
x[tl.n0,:,1] = x0[:,1]
v[tl.n0,:,0] = v0[:,0]
v[tl.n0,:,1] = v0[:,1]



######################
### Setting Plot
######################
# periodic
#nplot = 1000
nplot = 50
# view on-off
onoffview = True


pltfig = plt.figure(figsize=(0.1*(xmax-xmin),0.08*(ymax-ymin)))
axis = pltfig.add_subplot(1,1,1)
axis.clear()
title = 't = %000005.3f'%(tl.dt*tl.istep) + 'sec, ' + 'x = %00007.2f'%(x[tl.n0,0,0])+'m, ' + 'y = %00007.2f'%(x[tl.n0,0,1]) + 'm'
axis.set_title(title)
axis.set_xlabel('x ['+r'$m$'+']')
axis.set_ylabel('y ['+r'$m$'+']')
axis.set_xlim(xmin,xmax)
axis.set_ylim(ymin,ymax)
plt.ion()
pos, = axis.plot(x[tl.n0,:,0], x[tl.n0,:,1], 'ro') 
if testcase == 1 or testcase == 11:
   axis.axhline(y=x0[0,1], xmin=xmin, xmax=xmax)
   axis.axvline(x=xana,   ymin=ymin, ymax=ymax)
   pos2, = axis.plot(x_t[0], x_t[1], 'k--') 

harrsize = 100.0/(xmax-xmin)
arrmag   = 20.0/(xmax-xmin)
arr = []
for ipart in range(nparts):
   arr.append(axis.arrow(x[tl.n0,ipart,0], x[tl.n0,ipart,1], arrmag*v[tl.n0,ipart,0], arrmag*v[tl.n0,ipart,1], head_width=harrsize, head_length=harrsize, fc='k', ec='k'))

if onoffview: plt.show()


# Time Loop
for istep in np.arange(tl.nstep+1):

   if istep%nplot == 0:
      print '# istep = ', istep
      # Update Plot
      title = 't = %000005.3f'%(tl.dt*(tl.istep)) + 'sec, ' + 'x = %00007.2f'%(x[tl.n0,0,0])+'m, ' + 'y = %00007.2f'%(x[tl.n0,0,1]) + 'm'
      axis.set_title(title)
      pos.set_data(x[tl.n0,:,0], x[tl.n0,:,1])
      for ipart in range(nparts):
         arr[ipart].remove()
         arr[ipart] = axis.arrow(x[tl.n0,ipart,0], x[tl.n0,ipart,1], arrmag*v[tl.n0,ipart,0], arrmag*v[tl.n0,ipart,1], head_width=harrsize, head_length=harrsize, fc='k', ec='k')
      # START: saving figs
      txt_step = '%000005d' % int(istep/nplot)
      filename = namedir+'/'+txt_step+'.png'
      pltfig.savefig(filename)
      # END  : saving figs
      print ' '
      if onoffview:
         if istep == nstep:
            plt.show(True)
         else:
            plt.draw()

   # ComputeRHS
   solver.ComputeRHS(x[tl.n0,:,:],v[tl.n0,:,:],acc[tl.n0,:,:],xRHS,vRHS)

   # Apply the Time Integration
   tl.Update(x,xRHS)
   tl.Update(v,vRHS)
   tl.UpdateLevel()

   # Apply the Boundary Condition
   xbndry.CheckBoundary(tl,x,v,vRHS)
   ybndry.CheckBoundary(tl,x,v,vRHS)



