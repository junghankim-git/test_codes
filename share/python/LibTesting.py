import os
import sys
#from numpy.testing import assert_true
from numpy.testing import assert_equal
from numpy.testing import assert_almost_equal
from inspect import currentframe, getframeinfo

DIGIT = 15

#def isTrue(A_in):
#   assert_true(A_in)

def isEqual(A_in, B_in):
   assert_equal(A_in, B_in)

def isAlmostEqual(A_in, B_in, digit=DIGIT):
   # default DIGIT: single precision(5-7)
   assert_almost_equal(A_in, B_in, digit)
