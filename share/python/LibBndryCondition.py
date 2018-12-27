import os
import sys



class BndryCondition:

   def __init__(self,method,type,nparts,ndims,idim,min,max):
      self.method = method # 'Lagrangian' or 'Eularian'
      self.type   = type   # 'Reflection' or 'Periodic'
      self.nparts = nparts
      self.ndims  = ndims
      self.idim   = idim
      self.min    = min
      self.max    = max

   # Check the Boundary Condition
   def CheckBoundary(self, tl, position, velocity, vRHS):
      nparts = self.nparts
      idim   = self.idim
      n0     = tl.n0
      nm1    = tl.nm1
      if self.method == 'Lagrangian':

         if self.type == 'Periodic':
            for ipart in range(nparts):
               if position[n0][ipart][idim] < self.min:
                  position[n0][ipart][idim]  = self.max - (self.min - position[n0][ipart][idim])
                  if tl.uptype == 'Leapfrog':
                     position[nm1][ipart][idim] = self.max - (self.min - position[nm1][ipart][idim])
               if position[n0][ipart][idim] > self.max:
                  position[n0][ipart][idim]  = self.min + (position[n0][ipart][idim] - self.max)
                  if tl.uptype == 'Leapfrog':
                     position[nm1][ipart][idim] = self.min + (position[nm1][ipart][idim] - self.max)

         elif self.type == 'Reflection':
            for ipart in range(nparts):
               if position[n0][ipart][idim] < self.min:
                  position[n0][ipart][idim] = 2.0*self.min - position[n0][ipart][idim]
                  velocity[n0][ipart][idim] = -velocity[n0][ipart][idim]
                  if tl.uptype == 'Leapfrog':
                     position[nm1][ipart][idim] = 2.0*self.min - position[nm1][ipart][idim]
                     velocity[nm1][ipart][idim] = velocity[n0][ipart][idim]-2.0*tl.dt*vRHS[ipart][idim]
               if position[n0][ipart][idim] > self.max:
                  position[n0][ipart][idim] = 2.0*self.max - position[n0][ipart][idim]
                  velocity[n0][ipart][idim] = -velocity[n0][ipart][idim]
                  if tl.uptype == 'Leapfrog':
                     position[nm1][ipart][idim] = 2.0*self.max - position[nm1][ipart][idim]
                     velocity[nm1][ipart][idim] = velocity[n0][ipart][idim]-2.0*tl.dt*vRHS[ipart][idim]

         else:
            print '[TimeIntegration Class]: must be checked [type] ...'

      elif self.method == 'Eularian':
         print '[TimeIntegration Class]: Eularian Scheme was not developed...'
      else:
         print '[TimeIntegration Class]: must be checked [method]...'
