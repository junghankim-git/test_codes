#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
#
# Brief
#   - 
#-------------------------------------------------------------------------------

#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import cm
import os
import sys
import random
from optparse import OptionParser
SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
sys.path.append(SHARE_DIR)
from LibErrorHandler import *




Logger('test1')
Logger('test2', -1)
Logger('test3', -1, 'method')
Logger('test4', -2, 'method2')
Logger('test4', -2, 'method2')
