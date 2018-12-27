import os
import sys
import numpy as np
sys.path.append('/home/jhkim/Study/Library/Shared/Python')
from LibFDM import *


class Domain:
   def __init__(self,NumericalMethod=None):
      if NumericalMethod != None:
         self.method  = NumericalMethod.method
         self.ndims   = NumericalMethod.ndims
         self.dimsize = NumericalMethod.dimsize
         self.order   = NumericalMethod.order
         self.halo    = NumericalMethod.halo
      else:
         self.method  = None
         self.ndims   = None
         self.dimsize = None
         self.order   = None
         self.halo    = None

   def MakeDomain(self):
      ndims   = self.ndims
      dimsize = self.dimsize
      order   = self.order
      halo    = self.halo
      if ndims == 1:
         domain = [0.0 for i in range(dimsize+2*halo)]
      elif ndims == 2:
         domain = [[0.0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
      elif ndims == 3:
         domain = [[[0.0 for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
      return domain

   def MakeDomainValue(self,ndims_in,dimsize_in,order_in,halo_in):
      ndims   = ndims_in
      dimsize = dimsize_in
      order   = order_in
      halo    = halo_in
      if ndims == 1:
         domain = [0.0 for i in range(dimsize+2*halo)]
      elif ndims == 2:
         domain = [[0.0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
      elif ndims == 3:
         domain = [[[0.0 for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
      return domain
