#!/usr/bin/env python
import os
import sys
import random
import unittest
sys.path.append('/home/jhkim/Study/Library/Shared/Python')
from LibTimeIntegration import *

class TestSequenceFunctions(unittest.TestCase):

   def setUp(self):
      self.nparts = 2
      self.ndims  = 2
      self.dt     = 0.5
      self.nstep  = 10
      self.uptype = 'ForwardEuler'

   def test_shuffle(self):
      ti = TimeIntegration(self.nparts, self.ndims, self.dt, self.nstep, self.uptype)
      ti.UpdateLevel()
      self.assertFalse(ti.nm1 != 1 or ti.n0 != 2 or ti.np1 != 0)
      ti.UpdateLevel()
      self.assertFalse(ti.nm1 != 2 or ti.n0 != 0 or ti.np1 != 1)
      ti.UpdateLevel()
      self.assertFalse(ti.nm1 != 0 or ti.n0 != 1 or ti.np1 != 2)
      ti.UpdateLevel()
      self.assertFalse(ti.nm1 != 2 or ti.n0 != 2 or ti.np1 != 0)
      # make sure the shuffled sequence does not lose any elements
      self.seq = range(10)
      random.shuffle(self.seq)
      self.seq.sort()
      self.assertEqual(self.seq, range(10))

      # should raise an exception for an immutable sequence
      self.assertRaises(TypeError, random.shuffle, (1,2,3))

   def test_choice(self):
      self.seq = range(10)
      random.shuffle(self.seq)
      element = random.choice(self.seq)
      self.assertTrue(element in self.seq)

   def test_sample(self):
      self.seq = range(10)
      random.shuffle(self.seq)
      with self.assertRaises(ValueError):
         random.sample(self.seq, 20)
      for element in random.sample(self.seq, 5):
         self.assertTrue(element in self.seq)

if __name__ == '__main__':
   unittest.main()
