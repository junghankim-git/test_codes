#!/usr/bin/env python

import os
import sys
import random
import numpy as np
import matplotlib.pyplot as plt
SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
sys.path.append(SHARE_DIR)
from control_mesh import *
from LibDomain import *
from LibFDM    import *
from LibDifferentialOperator import *


def IniScalar(ndims_in, ss_in, ee_in, psi_in):
   if ndims_in == 1:
      for i in range(ss_in,ee_in+1):
         psi_in[i] = float(i)
   elif ndims_in == 2:
      for i in range(ss_in[0],ee_in[0]+1):
         for j in range(ss_in[1],ee_in[1]+1):
            psi_in[i][j] = float(i+j)
   elif ndims_in == 3:
      for i in range(ss_in[0],ee_in[0]+1):
         for j in range(ss_in[1],ee_in[1]+1):
            for k in range(ss_in[2],ee_in[2]+1):
               psi_in[i][j][k] = float(i+j+k)
   else:
      print 'error...'


ndims   = 2
dimsize = [5,5]
order   = 2
xmin    = -50
xmax    =  50
ymin    = -50
ymax    =  50
dx      = [(xmax-xmin)/dimsize[0], (ymax-ymin)/dimsize[1]]

FDM     = FDM(ndims, dimsize, dx, order)
halo, ss, ee, ts, te = FDM.GetInfoFDM()
# halo   : size of halo regine
# ss, ee : first and last index of the computaional array
# ts, te : first and last index of the entire array

domain  = Domain(FDM)
psi     = domain.MakeDomain()

#IniScalar(ndims, ss, ee, psi)
IniScalar(ndims, ts, te, psi)
print_mesh_2d(dimsize[0]+2*halo, dimsize[1]+2*halo, psi, message='psi')
FDM.BndryExchange(psi)
print_mesh_2d(dimsize[0]+2*halo, dimsize[1]+2*halo, psi, message='Bndry Exchange psi')

gpsi   = FDM.Gradient(psi)
mgpsi  = vector_to_matrix_2d(dimsize[0]+2*halo, dimsize[1]+2*halo, gpsi)
print_mesh_2d(dimsize[0]+2*halo, dimsize[1]+2*halo, mgpsi[0], message='gradient 0')
print_mesh_2d(dimsize[0]+2*halo, dimsize[1]+2*halo, mgpsi[1], message='gradient 1')

#dpsi   = FDM.ApplyFDM(psi)
#print_mesh_2d(dimsize[0]+2*halo, dimsize[1]+2*halo, dpsi, message='dpsi')

print '####### Check Differential Operator ########'
diffoper = DifferentialOperator(ndims,1,dimsize,[xmin,ymin],[xmax,ymax],'FDM',order,halo)
gpsi   = diffoper.Gradient(psi)
mgpsi  = vector_to_matrix_2d(dimsize[0]+2*halo, dimsize[1]+2*halo, gpsi)
print_mesh_2d(dimsize[0]+2*halo, dimsize[1]+2*halo, mgpsi[0], message='gradient 0')
print_mesh_2d(dimsize[0]+2*halo, dimsize[1]+2*halo, mgpsi[1], message='gradient 1')




