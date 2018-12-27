#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
sys.path.append(SHARE_DIR)
from class_kim_out import *

opt = OptionParser()
### action: 'store', 'store_const', 'append', 'count', 'callback'
opt.add_option('-v', '--version',  dest='version',  type='string',  default='',       action='store',      help='Version')
opt.add_option('-n', '--ne',       dest='ne',       type='int',     default=30,       action='store',      help='number of elements (default: ne30)')
opt.add_option('-c', '--compiler', dest='compiler', type='string',  default='gnu',    action='store',      help='Compiler (default: gnu)')
#
opt.add_option('-b', '--basedir',  dest='basedir',  type='string',  default='none',   action='store',      help='Base directory (/scratch/jhkim...) for version')
opt.add_option('-w', '--cversion', dest='cversion', type='string',  default='',       action='store',      help='Check Version')
opt.add_option('-s', '--specific', dest='tpath',    type='string',  default='',       action='store',      help='Specific target directory for cversion')
#
opt.add_option('-i', '--integration', dest='integration', type='string',  default='10h',    action='store', help='Integration (10h or 2d or all)')
opt.add_option('-d', '--dynamics', dest='dynamics', type='string',  default='SW',    action='store',      help='Dynamics (SH or SW or all)')
opt.add_option('-p', '--physics',  dest='physics',  type='string',  default='G',    action='store',      help='Physics  (G or K or all)')
opt.add_option('-o', '--option',   dest='option',   type='int',     default=0,        action='store',      help='Options (0:diff, 1:ncdump -h, )')
opt.add_option('-f', '--field',    dest='field',    type='string',  default='u',      action='store',      help='Field ')
opt.add_option('-l', '--list',     dest='ifile',    type='int',     default=-1,       action='store',      help='File index')

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

ne = options.ne
compiler = options.compiler

manager = kim_out(version, ne, compiler)



if options.cversion == '':
  cversion = version
else:
  cversion = options.cversion

integration = options.integration
basedir     = options.basedir
dynamics    = options.dynamics
physics     = options.physics
opt         = options.option
field       = options.field
ifile       = options.ifile
tpath       = options.tpath


manager.check_difference(path=basedir, version=cversion, integration=integration, dynamics=dynamics, physics=physics, opt=opt, var=field, ifile=ifile, tpath=tpath)
