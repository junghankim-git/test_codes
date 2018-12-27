#!/usr/bin/env python
import os
import sys
import time



testcases = [3, 5, 10]
orders    = [0, 1, 2]


#os.system('rm -rf ./figs/*')
nlfilename = 'inputs.nl'
for icase in testcases:
   for iorder in orders:
      os.system('rm -rf ./'+nlfilename)
      # STEP 1: start - namelist file
      infile = open('namelist/'+nlfilename, 'r')
      oufile = open('./'+nlfilename,      'w')
      for line in infile:
         if len(line.split()) > 0:
            if line.split()[0] == 'order':
               wline = 'order = '+str(iorder)
            elif line.split()[0] == 'testcase':
               wline = 'testcase = '+str(icase)
            else:
               wline = line
         else:
            wline = line
         oufile.write(wline)
      infile.close()
      oufile.close()
      os.system('./a.out; ./draw.py')
