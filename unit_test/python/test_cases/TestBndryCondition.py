#!/usr/bin/env python
import os
import sys
import random
import unittest
sys.path.append('/home/jhkim/work/shared/python')
from LibBndryCondition  import *
from LibTimeIntegration import *




class TestBndryCondition(unittest.TestCase):

   def setUp(self):
      ''' @brief this method is called before every test methods '''
      self.nparts = 3
      self.ndims  = 2
      self.idim  = 0
      self.min    = -10.0
      self.max    =  10.0

      self.dt     = 0.1
      self.nstep  = 2

   def tearDown(self):
      ''' @brief this method is called after every test methods '''
      pass

   def test01_LagrangianPeriodicForward(self):

      method = 'Lagrangian'
      uptype = 'ForwardEuler'
      tl = TimeIntegration(method, self.dt, self.nstep, uptype, self.ndims, self.nparts)

      # Test1: Lagrangian + Periodic
      type        = 'Periodic'
      Bndry = BndryCondition(method,type,self.nparts,self.ndims,self.idim,self.min,self.max)
      pos, vel, vRHS, res_pos, res_vel = self.VertifyPosVel(tl, method, type, uptype)
      Bndry.CheckBoundary(tl, pos, vel, vRHS)
      self.assertEqual(pos[tl.n0][:][:], res_pos[tl.n0][:][:])
      self.assertEqual(vel[tl.n0][:][:], res_vel[tl.n0][:][:])

   def test01_LagrangianReflectionForward(self):

      method = 'Lagrangian'
      uptype = 'ForwardEuler'
      tl = TimeIntegration(method, self.dt, self.nstep, uptype, self.ndims, self.nparts)

      # Test2: Lagrangian + Reflection
      type        = 'Reflection'
      Bndry = BndryCondition(method,type,self.nparts,self.ndims,self.idim,self.min,self.max)
      pos, vel, vRHS, res_pos, res_vel = self.VertifyPosVel(tl, method, type, uptype)
      Bndry.CheckBoundary(tl, pos, vel, vRHS)
      self.assertEqual(pos[tl.n0][:][:], res_pos[tl.n0][:][:])
      self.assertEqual(vel[tl.n0][:][:], res_vel[tl.n0][:][:])


   def VertifyPosVel(self,tl_in,method_in,type_in,uptype_in):
      dt     = self.dt
      nparts = self.nparts
      ndims  = self.ndims
      idim   = self.idim
      ntls   = tl_in.ntls
      
      x_out     = [[[-15.0,0.0], [0.0,0.0], [15.0,0.0]] for itl in range(ntls)]
      v_out     = [[[  2.0,0.0], [2.0,0.0], [ 2.0,0.0]] for itl in range(ntls)]
      if method_in == 'Lagrangian':
         if type_in == 'Periodic':
            vRHS_out  = [[1.0 for id in range(ndims)] for iparts in range(nparts)]
            res_x_out = [[[5.0,0.0], [0.0,0.0], [-5.0,0.0]] for itl in range(ntls)]
            res_v_out = [[[2.0,0.0], [2.0,0.0], [ 2.0,0.0]] for itl in range(ntls)]
         elif type_in == 'Reflection':
            vRHS_out  = [[1.0 for id in range(ndims)] for iparts in range(nparts)]
            res_x_out = [[[-5.0,0.0], [0.0,0.0], [5.0,0.0]] for itl in range(ntls)]
            res_v_out = [[[-2.0,0.0], [2.0,0.0], [-2.0,0.0]] for itl in range(ntls)]

      return x_out, v_out, vRHS_out, res_x_out, res_v_out

   
if __name__ == '__main__':
   unittest.main()
