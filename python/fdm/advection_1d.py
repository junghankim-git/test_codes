#!/usr/bin/env python

import os
import sys
import random
import numpy as np
from math import *
import matplotlib.pyplot as plt
SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
sys.path.append(SHARE_DIR)
from LibBndryCondition  import *
from LibTimeIntegration import *
from LibDifferentialOperator import *


def IniScalar(ndims_in, dimsize_in, halo_in):
   #opt = 0 # increment
   opt = 1 # gaussian
   if ndims_in == 1:
      psi_o = [0.0 for i in range(dimsize_in+2*halo_in)]
      for i in range(halo_in,dimsize_in+halo_in):
         if opt == 0:
            psi_o[i] = float(i)
         elif opt == 1:
            mu  = 50.0
            sig = 5.0
            amp = 5
            psi_o[i] = amp*exp(-pow(float(i)-mu,2.0)/2.0/sig/sig)
   elif ndims_in == 2:
      psi_o = [[0.0 for j in range(dimsize_in[1]+2*halo_in)] for i in range(dimsize_in[0]+2*halo_in)]
      for i in range(halo_in,dimsize_in[0]+halo_in):
         for j in range(halo_in,dimsize_in[1]+halo_in):
            psi_o[i][j] = float(i+j)
   elif ndims_in == 3:
      psi_o = [[[0.0 for k in range(dimsize_in[2]+2*halo_in)] for j in range(dimsize_in[1]+2*halo_in)] for i in range(dimsize_in[0]+2*halo_in)]
      for i in range(halo_in,dimsize_in[0]+halo_in):
         for j in range(halo_in,dimsize_in[1]+halo_in):
            for k in range(halo_in,dimsize_in[2]+halo_in):
               psi_o[i][j][k] = float(i+j+k)
   else:
      print 'error...'
   return psi_o



def ComputeRHS(ndims_in, dimsize_in, halo_in, psi_in, oper_in):
   gpsi = oper_in.Gradient(psi_in)
   for i in range(halo,dimsize+halo):
      gpsi[i] = -1.0*gpsi[i]
   return gpsi


# Domain
ndims   = 1
dimsize = 200
xmin    = -100
xmax    = 100
method  = 'FDM'
order   = 2
halo    = 1
oper    = DifferentialOperator(ndims,1,dimsize,xmin,xmax,method,order,halo)
dt      = 0.2
nsteps  = 100000
uptype  = 'Leapfrog'
tl      = TimeIntegration('Eulerian', dt,nsteps,uptype,ndims,dimsize,halo)
ntls    = tl.ntls

x        = [float(i-100) for i in range(dimsize+2*halo)]
psi        = [[0.0 for i in range(dimsize+2*halo)] for itl in range(ntls)]
psi[tl.n0] = IniScalar(ndims, dimsize, halo)

pltfig = plt.figure(figsize=(18,10))
pltfig.suptitle('Advection',fontsize=18)
axis = pltfig.add_subplot(1,1,1)
axis.clear()
axis.set_title('1d advection')
axis.set_xlim(xmin,xmax)
axis.set_ylim(0.0,6.0)
axis.set_xlabel('x')
axis.set_ylabel('psi')

plt.ion()

for istep in range(nsteps):
   oper.BndryExchange(psi[tl.n0])
   if istep == 0:
      plot, = axis.plot(x,psi[tl.n0],'o-',lw=1.0)
   else:
      plot.set_data(x,psi[tl.n0])
   plt.show()
   if istep == nsteps:
      plt.show(True)
   else:
      plt.draw()

   RHS = ComputeRHS(ndims, dimsize, halo, psi[tl.n0], oper)
   tl.Update(psi, RHS)
   tl.UpdateLevel()



