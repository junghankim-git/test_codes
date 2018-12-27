#!/usr/bin/env python

import os
import sys
SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
sys.path.append(SHARE_DIR)
from LibKIMOUT import *
from optparse import OptionParser

opt = OptionParser()
### action: 'store', 'store_const', 'append', 'count', 'callback'
opt.add_option('-v', '--version',  dest='version',  type='string',  default='',       action='store',      help='Version')
opt.add_option('-n', '--ne',       dest='ne',       type='int',     default=30,       action='store',      help='number of elements (default: ne30)')
opt.add_option('-c', '--compiler', dest='compiler', type='string',  default='gnu',    action='store',      help='Compiler (default: gnu)')
#opt.add_option('-m', '--make',     dest='make',                     default=False,    action='store_true', help='Make directory (def. False)')

opt.add_option('-s', '--source',   dest='source',   type='string',  default=None,     action='store',      help='Source version')
opt.add_option('-l', '--link',     dest='link',                     default=False,    action='store_true', help='on/off link')
opt.add_option('-o', '--copy',     dest='copy',                     default=False,    action='store_true', help='on/off copy')
opt.add_option('-d', '--dynamics', dest='dynamics', type='string',  default='all',    action='store',      help='Dynamics (SH or SW or all)')
opt.add_option('-p', '--physics',  dest='physics',  type='string',  default='all',    action='store',      help='Physics  (G or K or all)')
###
(options, args) = opt.parse_args()
##


nargs = len(args)

print options.version

if nargs == 0 and options.version == '':
  print 'check version'
  exit()

if options.version == '':
  version = args[0]
else:
  version = options.version
ne        = options.ne
compiler  = options.compiler

if os.path.isdir(version):
  print 'WARNING: directory was existed...'
#  exit()

manager = KIMOUT(version, ne, compiler)

manager.CreateDirs()


if options.link:
   if options.source == None:
      print 'FATAL: need source version (-s)'
      exit()
   if options.copy: options.copy = False
   manager.ActionDirs(options.source, options.dynamics, options.physics, 2)

if options.copy:
   manager.ActionDirs(options.source, options.dynamics, options.physics, 3)


