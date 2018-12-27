#!/usr/bin/env python
import os
import sys
import random
import unittest
sys.path.append('/home/jhkim/work/share/python')
from LibTimeIntegration import *




class TestTimeIntegration(unittest.TestCase):

   def setUp(self):
      ''' @brief this method is called before every test methods '''
      self.ndims  = 2
      self.dt     = 0.5
      self.nstep  = 10

      self.nparts = 2
      self.LagFEtl = TimeIntegration('Lagrangian', self.dt, self.nstep, 'ForwardEuler', self.ndims, self.nparts)
      self.LagLFtl = TimeIntegration('Lagrangian', self.dt, self.nstep, 'Leapfrog', self.ndims, self.nparts)

      self.dimsize = [2, 2]
      self.halo    = 1
      self.EulFEtl = TimeIntegration('Eulerian', self.dt, self.nstep, 'ForwardEuler', self.ndims, self.dimsize, self.halo)
      self.EulLFtl = TimeIntegration('Eulerian', self.dt, self.nstep, 'Leapfrog', self.ndims, self.dimsize, self.halo)

   def tearDown(self):
      ''' @brief this method is called after every test methods '''
      pass

   def test01_UpdateLevel(self):

      # lup = True: default o.tlon
      self.LagFEtl.UpdateLevel()
      self.assertTrue(self.LagFEtl.nm1 == 1 and self.LagFEtl.n0 == 2 and self.LagFEtl.np1 == 0)
      self.assertEqual(self.LagFEtl.istep, 1)

      self.LagFEtl.UpdateLevel()
      self.assertTrue(self.LagFEtl.nm1 == 2 and self.LagFEtl.n0 == 0 and self.LagFEtl.np1 == 1)
      self.assertEqual(self.LagFEtl.istep, 2)

      self.LagFEtl.UpdateLevel()
      self.assertTrue(self.LagFEtl.nm1 == 0 and self.LagFEtl.n0 == 1 and self.LagFEtl.np1 == 2)
      self.assertEqual(self.LagFEtl.istep, 3)

      # lup = False
      self.LagFEtl.UpdateLevel(lup=False)
      self.assertTrue(self.LagFEtl.nm1 == 1 and self.LagFEtl.n0 == 2 and self.LagFEtl.np1 == 0)
      self.assertEqual(self.LagFEtl.istep, 3)

   def test02_Update_LagrangianForward(self):
      x, v, xRHS, vRHS, res_x, res_v = self.VertifyValues(self.LagFEtl)
      self.LagFEtl.Update(x, xRHS)
      self.LagFEtl.Update(v, vRHS)
      self.assertEqual(x[self.LagFEtl.n0][:][:], res_x[self.LagFEtl.n0][:][:])
      self.assertEqual(v[self.LagFEtl.n0][:][:], res_v[self.LagFEtl.n0][:][:])

   def test02_Update_LagrangianLeapfrog(self):
      x, v, xRHS, vRHS, res_x, res_v = self.VertifyValues(self.LagLFtl)
      self.LagLFtl.Update(x, xRHS)
      self.LagLFtl.Update(v, vRHS)
      self.assertEqual(x[self.LagLFtl.n0][:][:], res_x[self.LagLFtl.n0][:][:])
      self.assertEqual(x[self.LagLFtl.nm1][:][:], res_x[self.LagLFtl.nm1][:][:])
      self.assertEqual(v[self.LagLFtl.n0][:][:], res_v[self.LagLFtl.n0][:][:])
      self.assertEqual(v[self.LagLFtl.nm1][:][:], res_v[self.LagLFtl.nm1][:][:])

   def test03_Update_EulerianForward(self):
      x, xRHS, res_x = self.VertifyValues(self.EulFEtl)
      self.EulFEtl.Update(x, xRHS)
      self.assertEqual(x[self.EulFEtl.n0][:][:], res_x[self.EulFEtl.n0][:][:])

   def test03_Update_EulerianLeapfrog(self):
      x, xRHS, res_x = self.VertifyValues(self.EulLFtl)
      self.EulLFtl.Update(x, xRHS)
      self.assertEqual(x[self.EulLFtl.n0][:][:], res_x[self.EulLFtl.n0][:][:])
      self.assertEqual(x[self.EulLFtl.nm1][:][:], res_x[self.EulLFtl.nm1][:][:])

   def VertifyValues(self,tl_i):
      dt     = self.dt
      ndims  = self.ndims
      ntls   = tl_i.ntls
      
      if tl_i.method == 'Lagrangian':
         nparts = self.nparts

         x_out     = [[[1.0 for idim in range(ndims)] for ipart in range(nparts)] for itl in range(ntls)]
         v_out     = [[[2.0 for idim in range(ndims)] for ipart in range(nparts)] for itl in range(ntls)]
         xRHS_out  = [[1.0 for idim in range(ndims)] for iparts in range(nparts)]
         vRHS_out  = [[2.0 for idim in range(ndims)] for iparts in range(nparts)]
         res_x_out = [[[1.0 for idim in range(ndims)] for ipart in range(nparts)] for itl in range(ntls)]
         res_v_out = [[[2.0 for idim in range(ndims)] for ipart in range(nparts)] for itl in range(ntls)]

         if tl_i.uptype == 'ForwardEuler':
            for ipart in range(nparts):
               for idim in range(ndims):
                  res_x_out[tl_i.np1][ipart][idim] = x_out[tl_i.n0][ipart][idim]+dt*xRHS_out[ipart][idim]
                  res_v_out[tl_i.np1][ipart][idim] = v_out[tl_i.n0][ipart][idim]+dt*vRHS_out[ipart][idim]
         elif tl_i.uptype == 'Leapfrog':
            for ipart in range(nparts):
               for idim in range(ndims):
                  res_x_out[tl_i.np1][ipart][idim] = x_out[tl_i.nm1][ipart][idim]+2.0*dt*xRHS_out[ipart][idim]
                  res_v_out[tl_i.np1][ipart][idim] = v_out[tl_i.nm1][ipart][idim]+2.0*dt*vRHS_out[ipart][idim]

         return x_out, v_out, xRHS_out, vRHS_out, res_x_out, res_v_out


      elif tl_i.method == 'Eulerian':
         dimsize = self.dimsize
         halo    = self.halo

         psi_out     = [[[1.0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)] for itl in range(ntls)]
         psiRHS      = [[1.0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)]
         res_psi_out = [[[1.0 for j in range(dimsize[1]+2*halo)] for i in range(dimsize[0]+2*halo)] for itl in range(ntls)]

         if tl_i.uptype == 'ForwardEuler':
            for i in range(halo,dimsize[0]+halo):
               for j in range(halo,dimsize[1]+halo):
                  res_psi_out[tl_i.np1][i][j] = psi_out[tl_i.n0][i][j]+dt*psiRHS[i][j]
         elif tl_i.uptype == 'Leapfrog':
            for i in range(halo,dimsize[0]+halo):
               for j in range(halo,dimsize[1]+halo):
                  res_psi_out[tl_i.np1][i][j] = psi_out[tl_i.nm1][i][j]+2.0*dt*psiRHS[i][j]

         return psi_out, psiRHS, res_psi_out

   
if __name__ == '__main__':
   unittest.main()
