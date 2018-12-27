import os
import sys
import numpy as np
sys.path.append('/home/jhkim/Study/Library/Shared/Python')


class DifferentialOperator:
   def __init__(self,ndims,nelems,dimsize,xmin,xmax,method,order,halo,spec=False):
      self.method  = method
      self.ndims   = ndims
      self.nelems  = nelems
      self.dimsize = dimsize
      self.xmin    = xmin
      self.xmax    = xmax
      self.order   = order
      self.halo    = halo
      if ndims == 1:
         self.emin   = -1.0
         self.emax   = 1.0
         self.dx     = 0.0
         self.D      = [0.0 for i in range(dimsize+2*halo)]
         self.DT     = [0.0 for i in range(dimsize+2*halo)]
         self.Dinv   = [0.0 for i in range(dimsize+2*halo)]
         self.DinvT  = [0.0 for i in range(dimsize+2*halo)]
         self.gij    = [0.0 for i in range(dimsize+2*halo)]
         self.gijinv = [0.0 for i in range(dimsize+2*halo)]
         self.J      = [0.0 for i in range(dimsize+2*halo)]
      elif ndims == 2:
         self.emin   = [-1.0, -1.0]
         self.emax   = [1.0, 1.0]
         self.dx     = [0.0, 0.0]
         mat0        = [[0.0,0.0],[0.0,0.0]]
         self.D      = [[mat0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         self.DT     = [[mat0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         self.Dinv   = [[mat0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         self.DinvT  = [[mat0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         self.gij    = [[mat0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         self.gijinv = [[mat0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         self.J      = [[0.0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
      elif ndims == 3:
         self.emin   = [-1.0, -1.0, -1.0]
         self.emax   = [1.0, 1.0, 1.0]
         self.dx     = [0.0, 0.0, 0.0]
         mat0        = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
         self.D      = [[[mat0 for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         self.DT     = [[[mat0 for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         self.Dinv   = [[[mat0 for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         self.DinvT  = [[[mat0 for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         self.gij    = [[[mat0 for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         self.gijinv = [[[mat0 for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         self.J      = [[[0.0 for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]

      emin = self.emin
      emax = self.emax
      # x = (xmax-xmin)/2(a+1)+xmin
      if method=='FDM':
         if spec == False:
            if ndims == 1:
               self.de  = (emax-emin)/(dimsize-1)
               self.rde = 1.0/self.de
               self.dx  = (xmax-xmin)/(dimsize-1)
               for i in range(halo,dimsize+halo):
                  self.D[i] = 0.5*(xmax-xmin)
                  self.DT[i],self.Dinv[i],self.DinvT[i],self.gij[i],self.gijinv[i],self.J[i]=self.GetMetricMatrix(self.D[i])
               self.BndryExchange(self.D)
               self.BndryExchange(self.DT)
               self.BndryExchange(self.Dinv)
               self.BndryExchange(self.DinvT)
               self.BndryExchange(self.gij)
               self.BndryExchange(self.gijinv)
               self.BndryExchange(self.J)
            elif ndims == 2:
               self.de  = [0.0 for idim in range(ndims)]
               self.rde = [0.0 for idim in range(ndims)]
               for idim in range(ndims):
                  self.de[idim]  = float(emax[idim]-emin[idim])/(dimsize[idim]-1)
                  self.rde[idim] = 1.0/self.de[idim]
                  self.dx[idim]  = float(xmax[idim]-xmin[idim])/(dimsize[idim]-1)
               for i in range(halo,dimsize[0]+halo):
                  for j in range(halo,dimsize[1]+halo):
                     self.D[i][j] = [[0.5*(self.xmax[0]-self.xmin[0]), 0.0],[0.0, 0.5*(self.xmax[1]-self.xmin[1])]]
                     self.DT[i][j],self.Dinv[i][j],self.DinvT[i][j],self.gij[i][j],self.gijinv[i][j],self.J[i][j]=self.GetMetricMatrix(self.D[i][j])
               self.BndryExchange(self.D)
               self.BndryExchange(self.DT)
               self.BndryExchange(self.Dinv)
               self.BndryExchange(self.DinvT)
               self.BndryExchange(self.gij)
               self.BndryExchange(self.gijinv)
               self.BndryExchange(self.J)
            elif ndims == 3:
               self.de = [0.0 for idim in range(ndims)]
               self.rde = [0.0 for idim in range(ndims)]
               for idim in range(ndims):
                  self.de[idim]  = float(emax[idim]-emin[idim])/(dimsize[idim]-1)
                  self.rde[idim] = 1.0/self.de[idim]
                  self.dx[idim]  = float(xmax[idim]-xmin[idim])/(dimsize[idim]-1)
               for i in range(halo,dimsize[0]+halo):
                  for j in range(halo,dimsize[1]+halo):
                     for k in range(halo,dimsize[2]+halo):
                        self.D[i][j][k] = [[0.5*(self.xmax[0]-self.xmin[0]), 0.0, 0.0],[0.0, 0.5*(self.xmax[1]-self.xmin[1]), 0.0], [0.0, 0.0, 0.5*(self.xmax[2]-self.xmin[2])]]
                        self.DT[i][j][k],self.Dinv[i][j][k],self.DinvT[i][j][k],self.gij[i][j][k],self.gijinv[i][j][k],self.J[i][j][k]=self.GetMetricMatrix(self.D[i][j][k])
               self.BndryExchange(self.D)
               self.BndryExchange(self.DT)
               self.BndryExchange(self.Dinv)
               self.BndryExchange(self.DinvT)
               self.BndryExchange(self.gij)
               self.BndryExchange(self.gijinv)
               self.BndryExchange(self.J)

#      self.D       = D_in
#      self.DT      = np.array(D_in).T.tolist()
#      self.Dinv    = np.linalg.inv(D_in).tolist()
#      self.DinvT   = np.linalg.inv(D_in).T.tolist()
#      self.gij     = np.dot(np.array(self.DT), np.array(self.D)).tolist()
#      self.gijinv  = np.linalg.inv(self.gij).tolist()
#      self.J       = np.abs(np.linalg.det(D_in)).tolist()


#   def GetMetricTensor(self):
#      return self.DT, self.Dinv, self.DinvT, self.gij, self.gijinv, self.J


#   def ComputaionalToPhysical(self,vec_in):
#      

   # private
   def GetMetricMatrix(self, D_in):
      if self.ndims == 1:
         D       = D_in
         DT      = D
         Dinv    = 1.0/D
         DinvT   = Dinv
         gij     = DT*D
         gijinv  = 1.0/D/D
         J       = np.sqrt(D*D)
      else:
         D       = D_in
         DT      = np.array(D_in).T.tolist()
         Dinv    = np.linalg.inv(D_in).tolist()
         DinvT   = np.linalg.inv(D_in).T.tolist()
         gij     = np.dot(np.array(DT), np.array(D)).tolist()
         gijinv  = np.linalg.inv(gij).tolist()
         J       = np.abs(np.linalg.det(D_in)).tolist()
      return DT, Dinv, DinvT, gij, gijinv, J


   def GetMatrix(self):
      return self.D, self.DT, self.Dinv, self.DinvT, self.gij, self.gijinv, self.J


   def Gradient(self, psi_in):
      ndims   = self.ndims
      dimsize = self.dimsize
      order   = self.order
      halo    = self.halo

      ref_vec = self.Gradient_ref(psi_in)

      if ndims == 1:
         vec_out = [0.0 for i in range(dimsize+2*halo)]
         for i in range(halo,dimsize+halo):
            vec_out[i] = self.DinvT[i]*ref_vec[i]
      elif ndims == 2:
         vec_out = [[[0.0,0.0] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         for i in range(halo,dimsize[0]+halo):
            for j in range(halo,dimsize[1]+halo):
               vec_out[i][j] = self.MatVecMul(self.DinvT[i][j], ref_vec[i][j])
      elif ndims == 3:
         vec_out = [[[[0.0,0.0,0.0] for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         for i in range(halo,dimsize[0]+halo):
            for j in range(halo,dimsize[1]+halo):
               for k in range(halo,dimsize[2]+halo):
                  vec_out[i][j][k] = self.MatVecMul(self.DinvT[i][j][k], ref_vec[i][j][k])
      return vec_out


   def Divergence(self, vec_in):
      ndims   = self.ndims
      dimsize = self.dimsize
      order   = self.order
      halo    = self.halo

      if ndims == 1:
         ref_vec = [0.0 for i in range(dimsize+2*halo)]
         #for i in range(halo,dimsize+halo):
         for i in range(dimsize+2*halo):
            ref_vec[i] = self.Dinv[i]*vec_in[i]
      elif ndims == 2:
         ref_vec = [[[0.0,0.0] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         #for i in range(halo,dimsize[0]+halo):
         #   for j in range(halo,dimsize[1]+halo):
         for i in range(dimsize[0]+2*halo):
            for j in range(dimsize[1]+2*halo):
               ref_vec[i][j] = self.MatVecMul(self.Dinv[i][j], vec_in[i][j])
      elif ndims == 3:
         ref_vec = [[[[0.0,0.0,0.0] for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         #for i in range(halo,dimsize[0]+halo):
         #   for j in range(halo,dimsize[1]+halo):
         #      for k in range(halo,dimsize[2]+halo):
         for i in range(dimsize[0]+2*halo):
            for j in range(dimsize[1]+2*halo):
               for k in range(dimsize[2]+2*halo):
                  ref_vec[i][j][k] = self.MatVecMul(self.Dinv[i][j][k], vec_in[i][j][k])


      psi = self.Divergence_ref(ref_vec)

      return psi

   def Gradient_ref(self, ref_psi_in):
      ndims   = self.ndims
      dimsize = self.dimsize
      order   = self.order
      halo    = self.halo
      rde     = self.rde
      if ndims == 1:
         ref_vec_out = [0.0 for i in range(dimsize+2*halo)]
      elif ndims == 2:
         ref_vec_out = [[[0.0,0.0] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
      elif ndims == 3:
         ref_vec_out = [[[[0.0,0.0,0.0] for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]


      if ndims == 1:
         for i in range(halo,dimsize+halo):
            ref_vec_out[i] = 0.5*rde*(ref_psi_in[i+1]-ref_psi_in[i-1])
      elif ndims == 2:
         for i in range(halo,dimsize[0]+halo):
            for j in range(halo,dimsize[1]+halo):
               ref_vec_out[i][j][0] = 0.5*rde[0]*(ref_psi_in[i+1][j]-ref_psi_in[i-1][j])
               ref_vec_out[i][j][1] = 0.5*rde[1]*(ref_psi_in[i][j+1]-ref_psi_in[i][j-1])
      elif ndims == 3:
         for i in range(halo,dimsize[0]+halo):
            for j in range(halo,dimsize[1]+halo):
               for k in range(halo,dimsize[2]+halo):
                  ref_vec_out[i][j][k][0] = 0.5*rde[0]*(ref_psi_in[i+1][j][k]-ref_psi_in[i-1][j][k])
                  ref_vec_out[i][j][k][1] = 0.5*rde[1]*(ref_psi_in[i][j+1][k]-ref_psi_in[i][j-1][k])
                  ref_vec_out[i][j][k][2] = 0.5*rde[2]*(ref_psi_in[i][j][k+1]-ref_psi_in[i][j][k-1])
      return ref_vec_out


   def Divergence_ref(self, ref_vec_in):
      ndims   = self.ndims
      dimsize = self.dimsize
      order   = self.order
      halo    = self.halo
      rde     = self.rde
      if ndims == 1:
         ref_psi_out = [0.0 for i in range(dimsize+2*halo)]
      elif ndims == 2:
         ref_psi_out = [[0.0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
      elif ndims == 3:
         ref_psi_out = [[[0.0 for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]


      if ndims == 1:
         for i in range(halo,dimsize+halo):
            ref_psi_out[i] = 1.0/self.J[i]*0.5*rde*((self.J[i+1]*ref_vec_in[i+1])-(self.J[i-1]*ref_vec_in[i-1]))
      elif ndims == 2:
         for i in range(halo,dimsize[0]+halo):
            for j in range(halo,dimsize[1]+halo):
               const = 1.0/self.J[i][j]
               xcomp = 0.5*rde[0]*(self.J[i+1][j]*ref_vec_in[i+1][j][0]-self.J[i-1][j]*ref_vec_in[i-1][j][0])
               ycomp = 0.5*rde[1]*(self.J[i][j+1]*ref_vec_in[i][j+1][1]-self.J[i][j-1]*ref_vec_in[i][j-1][1])
               ref_psi_out[i][j] = const*(xcomp+ycomp)
      elif ndims == 3:
         for i in range(halo,dimsize[0]+halo):
            for j in range(halo,dimsize[1]+halo):
               for k in range(halo,dimsize[2]+halo):
                  const = 1.0/self.J[i][j][k]
                  xcomp = 0.5*rde[0]*(self.J[i+1][j][k]*ref_vec_in[i+1][j][k][0]-self.J[i-1][j][k]*ref_vec_in[i-1][j][k][0])
                  ycomp = 0.5*rde[1]*(self.J[i][j+1][k]*ref_vec_in[i][j+1][k][1]-self.J[i][j-1][k]*ref_vec_in[i][j-1][k][1])
                  zcomp = 0.5*rde[2]*(self.J[i][j][k+1]*ref_vec_in[i][j][k+1][2]-self.J[i][j][k+1]*ref_vec_in[i][j][k+1][2])
                  ref_psi_out[i][j][k] = const*(xcomp+ycomp+zcomp)
      return ref_psi_out




   def MatVecMul(self, mat_in, vec_in):
      ndims   = self.ndims
      dimsize = self.dimsize
      halo    = self.halo
      if ndims == 1:
         vec_out = mat_in*vec_in
      else:
         vec_out = [0.0 for i in range(ndims)]
         for i in range(ndims):
            for j in range(ndims):
               vec_out[i] = vec_out[i] + mat_in[i][j]*vec_in[j]

      return vec_out




   def BndryExchange(self, psi_in):
      '''
        0        halo     dimsize
        |         |          |
      _________________________________________
      |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |
      |    |    |0000|0000|0000|0000|    |    |
      |                               |
      |                           dimsize+halo
      '''
      method  = self.method
      ndims   = self.ndims
      dimsize = self.dimsize
      halo    = self.halo
      if method=='FDM':
         if ndims == 1:
            for i in range(halo):
               psi_in[i]              = psi_in[dimsize+i]
               psi_in[dimsize+halo+i] = psi_in[halo+i]
         elif ndims == 2:
            for i in range(halo):
               psi_in[i][halo:halo+dimsize[1]]    = psi_in[dimsize[0]+i][halo:halo+dimsize[1]]
               psi_in[dimsize[0]+halo+i][halo:halo+dimsize[1]] = psi_in[halo+i][halo:halo+dimsize[1]]
            for jhalo in range(halo):
               for i in range(dimsize[0]):
                  psi_in[halo+i][dimsize[1]+halo+jhalo] = psi_in[halo+i][halo+jhalo]
                  psi_in[halo+i][jhalo] = psi_in[halo+i][dimsize[1]+jhalo]




