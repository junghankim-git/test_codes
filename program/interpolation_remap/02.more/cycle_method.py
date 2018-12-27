#!/usr/bin/env python
import os
import sys
import time



mono  = 1
bndry = 2
cases = [1,2,3,4,5,6]

#os.system('rm -rf ./figs/*')
for icase in cases:
   command = './methods.py -b {} -m {} -c {} -f'.format(bndry,mono,icase)
   print command
   os.system(command)
