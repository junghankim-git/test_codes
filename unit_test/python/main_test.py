#!/usr/bin/env python
import os
import sys
import random
import unittest
from optparse import OptionParser
sys.path.append('/home/jhkim/work/shared/python')
sys.path.append('/home/jhkim/work/unit_test/python/test_cases')
# Classes
from TestTimeIntegration      import *
from TestBndryCondition       import *
from test_control_mesh        import *
from TestDifferentialOperator import *

parser = OptionParser()
# action: 'store', 'store_const', 'append', 'count', 'callback'
parser.add_option('-q', '--quiet',   dest='quiet',   default=False, action='store_true', help='Addtional message.')
parser.add_option('-v', '--verbose', dest='verbose', default=2,     action='store',      help='Unittest\'s verbose.', metavar='LEVEL')
parser.add_option('-l', '--list',    dest='list',    default=False, action='store_true', help='Listing the testcases.')
parser.add_option('-c', '--case',    dest='cases',   type='int',    action='append',     help='Name of testcases for unittest.', metavar='CASE')

(options, args) = parser.parse_args()

ncases   = 4
casename = ['control_mesh','TimeIntegration', 'BndryExchange','DifferentialOperator']


if options.list:
   for icase in range(ncases):
      print 'case '+str(icase)+': '+casename[icase]
   quit()


cases    = [unittest.TestSuite() for icase in range(ncases)]
runcases = []


if options.cases == None:
   options.cases = [icase for icase in range(ncases)]
print options.cases


for icase in range(ncases):
   if icase in options.cases:
   
      if icase == 0:
         # Class:: control_mesh
         cases[icase].addTest(test_control_mesh('test01_direction_axis'))
         cases[icase].addTest(test_control_mesh('test02_symmetry_mesh'))
         cases[icase].addTest(test_control_mesh('test03_vector_to_matrix_2d'))
         runcases.append(cases[icase])

      if icase == 1:
         # Class:: TimeIntegration
         cases[icase].addTest(TestTimeIntegration('test01_UpdateLevel'))
         cases[icase].addTest(TestTimeIntegration('test02_Update_LagrangianForward'))
         cases[icase].addTest(TestTimeIntegration('test02_Update_LagrangianLeapfrog'))
         cases[icase].addTest(TestTimeIntegration('test03_Update_EulerianForward'))
         cases[icase].addTest(TestTimeIntegration('test03_Update_EulerianLeapfrog'))
         runcases.append(cases[icase])
   
      if icase == 2:
         # Class:: BndryCondition
         cases[icase].addTest(TestBndryCondition('test01_LagrangianPeriodicForward'))
         cases[icase].addTest(TestBndryCondition('test01_LagrangianReflectionForward'))
         runcases.append(cases[icase])
   
      if icase == 3:
         # Class:: DifferentialOperator
         cases[icase].addTest(TestDifferentialOperator('test01_OnlyMetric1D'))
         cases[icase].addTest(TestDifferentialOperator('test01_OnlyMetric2D'))
         cases[icase].addTest(TestDifferentialOperator('test01_OnlyMetric3D'))
         cases[icase].addTest(TestDifferentialOperator('test02_BndryExchange1D_FDM'))
         cases[icase].addTest(TestDifferentialOperator('test02_BndryExchange2D_FDM'))
         cases[icase].addTest(TestDifferentialOperator('test03_Gradient1D_FDM'))
         cases[icase].addTest(TestDifferentialOperator('test03_Gradient2D_FDM'))
         cases[icase].addTest(TestDifferentialOperator('test03_Gradient3D_FDM_check'))
         cases[icase].addTest(TestDifferentialOperator('test04_Divergence1D_FDM'))
         cases[icase].addTest(TestDifferentialOperator('test04_Divergence2D_FDM'))
         runcases.append(cases[icase])


allcase = unittest.TestSuite(runcases)
#print 'verbose = ', options.verbose
runner = unittest.TextTestRunner(verbosity=options.verbose)
runner.run(allcase)
