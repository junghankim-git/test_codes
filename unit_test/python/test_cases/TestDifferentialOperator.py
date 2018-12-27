#!/usr/bin/env python
import os
import sys
import random
import unittest
sys.path.append('/home/jhkim/work/share/python')
from LibTesting import *
from LibDifferentialOperator import *




class TestDifferentialOperator(unittest.TestCase):

   def setUp(self):
      ''' @brief this method is called before every test methods '''
      self.diff1D = self.VertifyDefault(1)
      self.diff2D = self.VertifyDefault(2)
      self.diff3D = self.VertifyDefault(3)

   def tearDown(self):
      ''' @brief this method is called after every test methods '''
      pass

   def test01_OnlyMetric1D(self):
      ndims = 1
      D_res, DT_res, Dinv_res, DinvT_res, gij_res, gijinv_res, J_res = self.VertifyMetric(ndims)
      DT, Dinv, DinvT, gij, gijinv, J = self.diff1D.GetMetricMatrix(D_res)
      isEqual(DT_res, DT)
      isAlmostEqual(Dinv_res, Dinv)
      isAlmostEqual(DinvT_res, DinvT)
      isEqual(gij_res, gij)
      isAlmostEqual(gijinv_res, gijinv, digit=14)
      isAlmostEqual(J_res, J)

   def test01_OnlyMetric2D(self):
      ndims = 2
      D_res, DT_res, Dinv_res, DinvT_res, gij_res, gijinv_res, J_res = self.VertifyMetric(ndims)
      DT, Dinv, DinvT, gij, gijinv, J = self.diff2D.GetMetricMatrix(D_res)
      isEqual(DT_res, DT)
      isAlmostEqual(Dinv_res, Dinv)
      isAlmostEqual(DinvT_res, DinvT)
      isEqual(gij_res, gij)
      isAlmostEqual(gijinv_res, gijinv, digit=13)
      isAlmostEqual(J_res, J)

   def test01_OnlyMetric3D(self):
      ndims = 3
      D_res, DT_res, Dinv_res, DinvT_res, gij_res, gijinv_res, J_res = self.VertifyMetric(ndims)
      DT, Dinv, DinvT, gij, gijinv, J = self.diff3D.GetMetricMatrix(D_res)
      isEqual(DT_res, DT)
      isAlmostEqual(Dinv_res, Dinv)
      isAlmostEqual(DinvT_res, DinvT)
      isEqual(gij_res, gij)
      isAlmostEqual(gijinv_res, gijinv, digit=13)
      isAlmostEqual(J_res, J)

   def test02_BndryExchange1D_FDM(self):
      ndims = 1
      psi, psi_res = self.VertifyBndryExchange_FDM(ndims)
      self.diff1D.BndryExchange(psi)
      isAlmostEqual(psi_res, psi)

   def test02_BndryExchange2D_FDM(self):
      ndims = 2
      psi, psi_res = self.VertifyBndryExchange_FDM(ndims)
      self.diff2D.BndryExchange(psi)
      isAlmostEqual(psi_res, psi)

   def test03_Gradient1D_FDM(self):
      ndims = 1
      psi, gpsi_res = self.VertifyGradient_FDM(ndims)
      self.diff1D.BndryExchange(psi)
      gpsi = self.diff1D.Gradient(psi)
      isAlmostEqual(gpsi_res, gpsi, digit=13)

   def test03_Gradient2D_FDM(self):
      ndims = 2
      psi, gpsi_res = self.VertifyGradient_FDM(ndims)
      self.diff2D.BndryExchange(psi)
      gpsi = self.diff2D.Gradient(psi)
      isAlmostEqual(gpsi_res, gpsi, digit=13)

   def test03_Gradient3D_FDM_check(self):
      ndims = 3
      psi, gpsi_res = self.VertifyGradient_FDM(ndims)
      self.diff3D.BndryExchange(psi)
      gpsi = self.diff3D.Gradient(psi)
      isAlmostEqual(gpsi_res, gpsi, digit=13)

   def test04_Divergence1D_FDM(self):
      ndims = 1
      vec, dvec_res = self.VertifyDivergence_FDM(ndims)
      self.diff1D.BndryExchange(vec)
      dvec = self.diff1D.Divergence(vec)
      isAlmostEqual(dvec_res, dvec, digit=13)

   def test04_Divergence2D_FDM(self):
      ndims = 2
      vec, dvec_res = self.VertifyDivergence_FDM(ndims)
      self.diff2D.BndryExchange(vec)
      dvec = self.diff2D.Divergence(vec)
      isAlmostEqual(dvec_res, dvec, digit=13)

   def VertifyDefault(self, ndims):
      nelems  = 1
      method  = 'FDM'
      order   = 2
      halo    = 1
      if ndims == 1:
         dimsize    = 3
         xmin       = 0.0
         xmax       = 2.0
      elif ndims == 2:
         dimsize    = [3,3]
         xmin       = [0.0, 0.0]
         xmax       = [2.0, 2.0]
      elif ndims == 3:
         dimsize    = [3,3,3]
         xmin       = [0.0, 0.0, 0.0]
         xmax       = [2.0, 2.0, 2.0]
      #diffM = DifferentialOperator(ndims,nelems,dimsize,xmin,xmax,method,order,halo,True)
      diff  = DifferentialOperator(ndims,nelems,dimsize,xmin,xmax,method,order,halo)
      return diff

   def VertifyMetric(self, ndims):

      if ndims == 1:
         D_res_out      = 2.0
         DT_res_out     = 2.0
         Dinv_res_out   = 0.5
         DinvT_res_out  = 0.5
         gij_res_out    = 4.0
         gijinv_res_out = 0.25
         J_res_out      = 2.0
      elif ndims == 2:
         D_res_out      = [[2.0,1.0],[3.0,2.0]]
         DT_res_out     = [[2.0,3.0],[1.0,2.0]]
         Dinv_res_out   = [[2.0,-1.0],[-3.0,2.0]]
         DinvT_res_out  = [[2.0,-3.0],[-1.0,2.0]]
         gij_res_out    = [[13.0,8.0],[8.0,5.0]]
         gijinv_res_out = [[5.0,-8.0],[-8.0,13.0]]
         J_res_out      = 1.0
      elif ndims == 3:
         D_res_out      = [[2.0,1.0,0.0],[0.0,3.0,2.0],[0.0,1.0,1.0]]
         DT_res_out     = [[2.0,0.0,0.0],[1.0,3.0,1.0],[0.0,2.0,1.0]]
         Dinv_res_out   = [[0.5,-0.5,1.0],[0.0,1.0,-2.0],[0.0,-1.0,3.0]]
         DinvT_res_out  = [[0.5,0.0,0.0],[-0.5,1.0,-1.0],[1.0,-2.0,3.0]]
         gij_res_out    = [[4.0,2.0,0.0],[2.0,11.0,7.0],[0.0,7.0,5.0]]
         gijinv_res_out = [[1.5,-2.5,3.5],[-2.5,5.0,-7.0],[3.5,-7.0,10.0]]
         J_res_out      = 2.0


      return D_res_out, DT_res_out, Dinv_res_out, DinvT_res_out, gij_res_out, gijinv_res_out, J_res_out


   def VertifyBndryExchange_FDM(self, ndims):

      if ndims == 1:
         before_out = [0.0,1.0,2.0,3.0,0.0]
         after_out  = [3.0,1.0,2.0,3.0,1.0]
      elif ndims == 2:
         before_out = [[0.0,0.0,0.0,0.0,0.0],[0.0,1.0,2.0,3.0,0.0],[0.0,4.0,5.0,6.0,0.0],[0.0,7.0,8.0,9.0,0.0],[0.0,0.0,0.0,0.0,0.0]]
         after_out  = [[0.0,7.0,8.0,9.0,0.0],[3.0,1.0,2.0,3.0,1.0],[6.0,4.0,5.0,6.0,4.0],[9.0,7.0,8.0,9.0,7.0],[0.0,1.0,2.0,3.0,0.0]]
      elif ndims == 3:
         before_out = [[0.0,0.0,0.0,0.0,0.0],[0.0,1.0,2.0,3.0,0.0],[0.0,4.0,5.0,6.0,0.0],[0.0,7.0,8.0,9.0,0.0],[0.0,0.0,0.0,0.0,0.0]]
         after_out  = [[0.0,7.0,8.0,9.0,0.0],[3.0,1.0,2.0,3.0,1.0],[6.0,4.0,5.0,6.0,4.0],[9.0,7.0,8.0,9.0,7.0],[0.0,1.0,2.0,3.0,0.0]]
      return before_out, after_out


   def VertifyGradient_FDM(self, ndims):


      if ndims == 1:
         dimsize = self.diff1D.dimsize
         halo    = self.diff1D.halo
         psi_out        = [0.0, 0.0, 1.0, 0.0, 0.0]
         gpsi_out       = [0.0, 0.5, 0.0,-0.5, 0.0]
      elif ndims == 2:
         dimsize = self.diff2D.dimsize
         halo    = self.diff2D.halo
         psi_out        = [[0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0]]
         gpsi_out       = [[[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0]],[[0.0,0.0],[0.0,0.5],[0.0,0.0],[0.0,-0.5],[0.0,0.0]],\
                           [[0.0,0.0],[0.0,0.5],[0.0,0.0],[0.0,-0.5],[0.0,0.0]],[[0.0,0.0],[0.0,0.5],[0.0,0.0],[0.0,-0.5],[0.0,0.0]],\
                           [[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0]]]
      elif ndims == 3:
         dimsize = self.diff3D.dimsize
         halo    = self.diff3D.halo
         psi_out   = [[[0.0 for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] \
                     for i in range(dimsize[0]+2*halo)]
         gpsi_out  = [[[[0.0,0.0,0.0] for k in range(dimsize[2]+2*halo)] for j in range(dimsize[1]+2*halo)] \
                     for i in range(dimsize[0]+2*halo)]
         for i in range(dimsize[0]+2*halo):
            for j in range(dimsize[1]+2*halo):
               for k in range(dimsize[2]+2*halo):
                  psi_out[i][j][k]  = float(i+j+k)
         for i in range(halo,dimsize[0]+halo):
            for j in range(halo,dimsize[1]+halo):
               for k in range(halo,dimsize[2]+halo):
                  gpsi_out[i][j][k] = [1.0,1.0,1.0]
      return psi_out, gpsi_out

   def VertifyDivergence_FDM(self, ndims):

      if ndims == 1:
         dimsize  = self.diff1D.dimsize
         halo     = self.diff1D.halo
         vec_out  = [0.0, 0.0, 1.0, 2.0, 0.0]
         dvec_out = [0.0, -0.5, 1.0,-0.5, 0.0]
      elif ndims == 2:
         dimsize  = self.diff2D.dimsize
         halo     = self.diff2D.halo
         vec_out  = [[[0.0,0.0] for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         dvec_out = [[0.0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         for i in range(halo,dimsize[0]+halo):
            for j in range(halo,dimsize[1]+halo):
               vec_out[i][j][0] = float(i+1)
               vec_out[i][j][1] = float(j+1)
               dvec_out[i][j]   = 2.0
               if i == halo:              dvec_out[i][j] = -1.0
               if i == dimsize[0]+halo-1: dvec_out[i][j] = -1.0
               if j == halo:              dvec_out[i][j] = -1.0
               if j == dimsize[1]+halo-1: dvec_out[i][j] = -1.0
               if i == halo and j == 2:              dvec_out[i][j] = 0.5
               if i == dimsize[0]+halo-1 and j == 2: dvec_out[i][j] = 0.5
               if j == halo and i == 2:              dvec_out[i][j] = 0.5
               if j == dimsize[1]+halo-1 and i == 2: dvec_out[i][j] = 0.5
      elif ndims == 3:
         xmax_out       = [2.0, 2.0, 2.0]
      return vec_out, dvec_out



if __name__ == '__main__':
   unittest.main()
