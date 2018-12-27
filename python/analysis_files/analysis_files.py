#!/usr/bin/env python
#-------------------------------------------------------------------------------
# MODULE KiapsGrid
#
#> @brief
#>  - show status of batch jobs
#>
#> @date 13JAN2015
#>  - JH KIM(jh.kim@kiaps.org) : First written

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import cm
import os
import sys
import random
from optparse import OptionParser
#SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
SHARE_DIR='/home/jhkim/work/share/python'
sys.path.append(SHARE_DIR)
from class_analysis_files import *

opt = OptionParser()

# action: 'store', 'store_const', 'append', 'count', 'callback'
opt.add_option('-d', '--dir',     dest='dir',       default='./',  action='store',       help='searching drectory')
opt.add_option('-e', '--ext',     dest='ext',       default=[],    action='append',      help='extension of files')
opt.add_option('-x', '--exc',     dest='exc',       default=['.svn','.mod','.o'], action='append',      help='excluded files')
opt.add_option('-q', '--quiet',   dest='on_quiet',   default=True,    action='store_false', help='print files')
opt.add_option('-p', '--print',   dest='on_print_ptn',  default=True,  action='store_false', help='print pattern')
opt.add_option('-r', '--replace', dest='on_replace_ptn', default=False, action='store_true',  help='replace word')
opt.add_option('-m', '--move',    dest='on_move_file',    default=False, action='store_true',  help='move file')

(options, args) = opt.parse_args()


#print options
#print args

#if not options.

#if len(args) != 1:
#  print 'check argument....'
#  quit()

print 'Options = ', options
print 'Args    = ', args

analysis = class_analysis_files(options.dir,options.ext,options.exc,on_quiet=options.on_quiet,on_print_ptn=options.on_print_ptn,on_replace_ptn=options.on_replace_ptn,on_move_file=options.on_move_file)
analysis.run_process(args)

print ' '
print 'for Args :', args
