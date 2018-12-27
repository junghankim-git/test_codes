#!/usr/bin/env python
#-------------------------------------------------------------------------------
# MODULE KiapsGrid
#
#> @brief
#>  - show status of batch jobs
#>
#> @date 13JAN2015
#>  - JH KIM(jh.kim@kiaps.org) : First written

import os
import sys
from optparse import OptionParser
#SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
#sys.path.append(SHARE_DIR)
sys.path.append('/home/jhkim/work/share/python')
from class_diff_files import *

opt = OptionParser()

# action: 'store', 'store_const', 'append', 'count', 'callback'
opt.add_option('-s', '--sdir',    dest='sdir',      default='./',  action='store',       help='Searching source drectory')
opt.add_option('-t', '--tdir',    dest='tdir',      default='./',  action='store',       help='Searching target drectory')
opt.add_option('-e', '--ext',     dest='ext',       default=[],    action='append',      help='Extension of files')
opt.add_option('-x', '--xext',    dest='exc',    default=['.svn'],    action='append',      help='Exception of files')
opt.add_option('-q', '--quiet',   dest='isQuiet',   default=True,  action='store_false', help='Print files')
opt.add_option('-p', '--pattern', dest='isPrtPat',  default=True,  action='store_false', help='Print pattern')
opt.add_option('-m', '--move',    dest='isMove',    default=False, action='store_true',  help='Move file')

(options, args) = opt.parse_args()


#print options
#print args

#if not options.

#if len(args) != 1:
#  print 'check argument....'
#  quit()

print 'Options = ', options

analysis = diff_files(options.sdir, options.tdir, options.ext, options.exc, onQuiet=options.isQuiet, onPrtPattern=options.isPrtPat, isMoveFile=options.isMove)
analysis.run_process(args)

print ' '
