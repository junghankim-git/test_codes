#!/usr/bin/env python
import os
import sys
import random
import unittest
sys.path.append('/home/jhkim/work/share/python')
from control_mesh import *




class test_control_mesh(unittest.TestCase):


   def setUp(self):
      ''' @brief this method is called before every test methods '''
      self.nparts = 3
      self.ndims  = 2
      self.idim   = 0
      self.min    = -10.0
      self.max    =  10.0

      self.dt     = 0.1
      self.nstep  = 2

   def tearDown(self):
      ''' @brief this method is called after every test methods '''
      pass

   def test01_direction_axis(self):
      w, e, s, n, ws, wn, es, en = get_direction_index()
      x, y, xy, mxy, org         = get_axis_index()
      self.assertEqual(w,   0)
      self.assertEqual(e,   1)
      self.assertEqual(s,   2)
      self.assertEqual(n,   3)
      self.assertEqual(ws,  4)
      self.assertEqual(es,  5)
      self.assertEqual(wn,  6)
      self.assertEqual(en,  7)
      self.assertEqual(x,   0)
      self.assertEqual(y,   1)
      self.assertEqual(xy,  2)
      self.assertEqual(mxy, 3)
      self.assertEqual(org, 4)


   def test02_symmetry_mesh(self):
      for iaxis in range(5):
         n, mat, res_mat = self.vertify_symmetry_mesh(iaxis)
         self.assertEqual(res_mat, symmetry_mesh(n,iaxis,mat))

   def test03_vector_to_matrix_2d(self):
      nx, ny, vec, res_mat = self.vertify_vector_to_matrix_2d()
      self.assertEqual(res_mat, vector_to_matrix_2d(nx,ny,vec))

   def vertify_symmetry_mesh(self,axis_in):

      mat     = [[1.0, 2.0],\
                 [3.0, 4.0]]
      res_mat = [[0.0, 0.0],\
                 [0.0, 0.0]]
      if axis_in == 0:
         res_mat = [[2.0, 1.0],\
                    [4.0, 3.0]]
      elif axis_in == 1:
         res_mat = [[3.0, 4.0],\
                    [1.0, 2.0]]
      elif axis_in == 2:
         res_mat = [[1.0, 3.0],\
                    [2.0, 4.0]]
      elif axis_in == 3:
         res_mat = [[4.0, 2.0],\
                    [3.0, 1.0]]
      elif axis_in == 4:
         res_mat = [[4.0, 3.0],\
                    [2.0, 1.0]]

      return 2, mat, res_mat 



   def vertify_vector_to_matrix_2d(self):

      nx_out = 2
      ny_out = 2

      vec_out  = [[[1.0,2.0], [3.0,4,0]],\
                  [[5.0,6.0], [7.0,8.0]]]
      res_mat1 = [[1.0, 3.0],\
                  [5.0, 7.0]]
      res_mat2 = [[2.0, 4.0],\
                  [6.0, 8.0]]
      res_mat_out = [res_mat1, res_mat2]
      return nx_out, ny_out, vec_out, res_mat_out

   
if __name__ == '__main__':
   unittest.main()
