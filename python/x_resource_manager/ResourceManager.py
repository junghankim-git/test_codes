#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
#
# Brief
#   - 
#
# Date 24MAR2015
#   - JH KIM(jh.kim@kiaps.org) : 초기 작성
#
# Date: 01FEB2015
#   - JH KIM: Functions for viewing difference between resources
#                       and setting value of variables for specific resource
#-------------------------------------------------------------------------------

import sys
from optparse import OptionParser
SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
sys.path.append(SHARE_DIR)
from LibResourceManager import *


opt = OptionParser()
## action: 'store', 'store_const', 'append', 'count', 'callback'
opt.add_option('-n', '--names',   dest='names',      type='string', action='append',             help='Names for profiling')
opt.add_option('-i', '--input',   dest='inputpath',  type='string', action='store',              help='Input path')
opt.add_option('-o', '--output',  dest='outputpath', type='string', action='store',              help='Output path')
opt.add_option('-e', '--ne',      dest='ne',         type='int', default=30, action='store',     help='Number of elements')
opt.add_option('-p', '--print',   dest='isPrint',    default=False, action='store_true',         help='Print Resources')
opt.add_option('-w', '--wiki',    dest='isWiki',     default=False, action='store_true',         help='Print Wiki version')
opt.add_option('-d', '--diff',    dest='isDiff',     default=False, action='store_true',         help='Difference')
opt.add_option('-g', '--gen',     dest='isGen',      default=False, action='store_true',         help='Generate resources')
opt.add_option('-r', '--replace', dest='replace',    type='string', action='append',             help='')

##
(options, args) = opt.parse_args()
#

def Error(string='', values=None):
   if values == None:
      print string
   else:
      print string, values
   quit()


#
# Before (rcstring): 
#    ['Shared:A:1', 'Shared:B:2,3', 'KiapsHOMME:C:4']
# After  (names, variables, nvalues, values): 
#    ['Shared', 'Shared', 'KiapsHOMME'], ['A','B','C'], [1,2,1], [['1'],['2','3'],['4']]
def ParsingRCs(rcstring):
   names     = []
   variables = []
   nvalues   = []
   values    = []
   nrcs = len(rcstring)
   for i in range(nrcs):
      string = rcstring[i].split(':')
      ns     = len(string)
      if ns != 3:
         Error('check the resource options...')
      names.append(string[0])
      variables.append(string[1])
      val  = string[2].split(',')
      nval = len(val)
      nvalues.append(nval)
      values.append(val)
   return names, variables, values


# Before (names, variables, nvalues, values): 
#    ['Shared', 'Shared', 'KiapsHOMME'], ['A','B','C'], [1,2,1], [['1'],['2','3'],['4']]
# After  (names, variables, nvalues, values): 
#    ['Shared', 'KiapsHOMME'], [['A','B'],'C'], [[1,2],1], [[['1'],['2','3']],['4']]
def ReducingRCs(names, variables, values):
   rnames     = []
   rvariables = []
   rvalues    = []

   nrcs = len(names)
   for irc in range(nrcs):
      n = rnames.count(names[irc])
      if n == 0:
         rnames.append(names[irc])
         rvariables.append([variables[irc]])
         rvalues.append([values[irc]])
      elif n == 1:
         idx = rnames.index(names[irc])
         rvariables[idx].append(variables[irc])
         rvalues[idx].append(values[irc])
      else:
         Error('Error in ReducingRCs...', n)
   return rnames, rvariables, rvalues

   


# Make a profiler of KIM
if options.inputpath == None:
   rcmng = ResourceManager(names=options.names, ne=options.ne)
   rcmng.IniDefault()
else:
   rcmng = ResourceManager(names=options.names)
   rcmng.IniFiles(options.inputpath)

# Specfic variables and values
if options.replace!= None:
   rcmng.RcStringToSetValue(options.replace)
#   names, variables, values = ParsingRCs(options.replace)
#   rnames, rvariables, rvalues = ReducingRCs(names, variables, values)
#   rcmng.SetRCandValues(rnames, rvariables, rvalues)



# Print Resource
if options.isPrint:
   rcmng.PrintResources()


# Print Wiki
if options.isWiki:
   rcmng.PrintWiki()


# Check the difference between resources
nargs = len(args)  # arguments was list
if options.isDiff:
   if nargs == 0:
      # need to coded with error handle
      #path1 = '/home/jhkim/Study/Library/Python/ResourceManager/test1'
      #path2 = '/home/jhkim/Study/Library/Python/ResourceManager/test2'
      #rcmng.CheckResources(path1, path2)
      Error('need one more arguments for...')
   elif nargs == 1:
      rcmng.CheckResources(args[0])
   elif nargs == 2:
      rcmng.CheckResources(args[0], args[1])
   else:
      Error('too many arguments ...', args[2:])
      rcmng.CheckResources(args[0], args[1])
else:
   if nargs > 0:
      print 'arguments was ignored :', args

if options.isGen:
   if options.outputpath == None:
      rcmng.GenResources('./')
   else:
      rcmng.GenResources(options.outputpath)


print ' '
print ' '
print '------------------------------------------------------------------------------'
print '# Contact: Junghan Kim (jh.kim@kiaps.org)'
print '------------------------------------------------------------------------------'
