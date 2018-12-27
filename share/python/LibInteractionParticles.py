import os
import sys
import numpy as np
sys.path.append('/home/jhkim/Study/Library/Shared/Python')


class InteractionParticles:
   def __init__(self,nparts,ndims,type,lvis=False):
      self.type   = type   # Gravity, BetweenParticles
      self.nparts = nparts
      self.ndims  = ndims
      #self.mag    = 1500.0
      self.mag    = 2000.0
      self.lvis   = lvis

   def ComputeRHS(self, position, velocity, acceleration, xRHS, vRHS):
      gravity = -9.80665
      ndims   = self.ndims
      nparts  = self.nparts
      acceleration = [[0.0 for j in range(ndims)] for i in range(nparts)]
      direction    = [[[0.0 for k in range(ndims)] for j in range(nparts)] for i in range(nparts)]
      distance     = [[0.0 for j in range(nparts)] for i in range(nparts)]
      if self.lvis:
         viscosity_l = self.GetViscosity(velocity)
      if self.type == 'Gravity':
         for ipart in range(nparts):
            acceleration[ipart][:]  = [0.0 for i in range(ndims)]
            acceleration[ipart][-1] = gravity
      elif self.type == 'BetweenParticles':
         for ipart in range(nparts):
            acceleration[ipart][:]  = [0.0 for i in range(ndims)]
            for jpart in range(nparts):
               if ipart == jpart:
                  direction[ipart][jpart][:] = [0.0 for i in range(ndims)]
                  distance[ipart][jpart]     = 0.0
               else:
                  direction[ipart][jpart][:] = [position[jpart][i]-position[ipart][i] for i in range(ndims)]
                  distance[ipart][jpart]     = self.GetDistance(direction[ipart][jpart][:])
#                  print distance[ipart][jpart]
                  for idim in range(ndims):
                     acceleration[ipart][idim] = acceleration[ipart][idim] - self.mag/pow(distance[ipart][jpart],3.0)*direction[ipart][jpart][idim]
      else:
         print '[InteractionParticles Class]: must be checked [type]...'

      if self.lvis:
         for ipart in range(nparts):
            for idim in range(ndims):
               acceleration[ipart][idim] = acceleration[ipart][idim] + viscosity_l[ipart][idim]

      for ipart in range(nparts):
         for idim in range(ndims):
            xRHS[ipart][idim] = velocity[ipart][idim]
            vRHS[ipart][idim] = acceleration[ipart][idim]

#   def ComputeAcc(self, position, velocity, acceleration):
#      if self.type == 'Gravity':


   def GetDistance(self, direction_in):
      distance_out = 0.0
      for idim in range(self.ndims):
         distance_out = distance_out + direction_in[idim]*direction_in[idim]
      return np.sqrt(distance_out)

   def GetViscosity(self, velocity_in):
      vis = 0.01
      viscosity_out = [[0.0 for j in range(self.ndims)] for i in range(self.nparts)]
      for ipart in range(self.nparts):
         distance_l = self.GetDistance(velocity_in[ipart][:])
         for idim in range(self.ndims):
            viscosity_out[ipart][idim] = -vis*distance_l*velocity_in[ipart][idim]
      return viscosity_out
