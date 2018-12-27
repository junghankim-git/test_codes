import os
import sys
import numpy as np


class FDM:
   def __init__(self,ndims,dimsize,dx,order):
      self.method  = 'FDM'
      self.ndims   = ndims
      self.dimsize = dimsize
      self.dx      = dx
      if order%2 != 0:
         print 'FDM (__init__): check order...'
         quit()
      self.order  = order
      self.halo   = order/2
      self.StartStop()

   def StartStop(self):
      if self.ndims == 1:
         self.ss = self.halo
         self.ee = self.dimsize+self.halo-1
         self.ts = 0
         self.te = self.dimsize+2*self.halo-1
      elif self.ndims == 2:
         self.ss = [self.halo, self.halo]
         self.ee = [self.dimsize[0]+self.halo-1, self.dimsize[1]+self.halo-1]
         self.ts = [0, 0]
         self.te = [self.dimsize[0]+2*self.halo-1, self.dimsize[1]+2*self.halo-1]
      elif self.ndims == 3:
         self.ss = [self.halo, self.halo, self.halo]
         self.ee = [self.dimsize[0]+self.halo-1, self.dimsize[1]+self.halo-1, self.dimsize[2]+self.halo-1]
         self.ts = [0, 0, 0]
         self.te = [self.dimsize[0]+2*self.halo-1, self.dimsize[1]+2*self.halo-1, self.dimsize[2]+2*self.halo-1]

   def GetInfoFDM(self):
      return self.halo, self.ss, self.ee, self.ts, self.te


   def Gradient(self, psi_in):
      ndims    = self.ndims
      dimsize  = self.dimsize
      dx       = self.dx
      order    = self.order
      halo     = self.halo
      ss       = self.ss
      ee       = self.ee
      #self.BndryExchage(psi_in)

      if ndims == 1:
         rdx      = 1.0/dx
      elif ndims == 2:
         rdx      = [1.0/dx[0], 1.0/dx[1]]
      elif ndims == 3:
         rdx      = [1.0/dx[0], 1.0/dx[1], 1.0/dx[2]]

      psi_out = [[[0.0 for k in range(2)] for j in range(self.te[1]+1)] for i in range(self.te[0]+1)]
      for i in range(ss[0],ee[0]+1):
         for j in range(ss[1],ee[1]+1):
            psi_out[i][j][0] = 0.5*rdx[0]*(psi_in[i+1][j]-psi_in[i-1][j])
            psi_out[i][j][1] = 0.5*rdx[1]*(psi_in[i][j+1]-psi_in[i][j-1])

      return psi_out
      

   def BndryExchange(self, psi_in):
      ndims   = self.ndims
      dimsize = self.dimsize
      halo    = self.halo
      if ndims == 1:
         psi_in[0:halo-1] = psi_in[dimsize:dimsize+halo-1]
      elif ndims == 2:
         for i in range(halo):
            psi_in[i][:] = psi_in[i+dimsize[0]][:]
            psi_in[dimsize[0]+i+1][:] = psi_in[i+1][:]
         for i in range(self.te[0]):
            for j in range(halo):
               psi_in[i][j] = psi_in[i][j+dimsize[1]]
               psi_in[i][dimsize[1]+j+1] = psi_in[i][j+1]


