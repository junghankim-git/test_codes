#!/usr/bin/env python
import os
import sys
import random
import unittest
sys.path.append('/home/jhkim/work/share/python')
from class_xxxx import *




class Test????(unittest.TestCase):

    def setUp(self):
        ''' @brief this method is called before every test methods '''
        pass

    def tearDown(self):
        ''' @brief this method is called after every test methods '''
        pass

    def test_01_xxx(self):
        self.assertEqual(0,0)

   
if __name__ == "__main__":
    unittest.main()
