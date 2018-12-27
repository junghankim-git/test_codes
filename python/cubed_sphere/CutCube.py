#!/usr/bin/env python
from numpy import *
from scipy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
import matplotlib.patches as patches
import time
from matplotlib import cm
import os
import sys
#font
import matplotlib.font_manager as fontm

# load My Class
sys.path.append('/home/jhkim/Study/Library/Shared/Python')

#os.system('rm -rf ./FigCube')
#os.system('mkdir ./FigCube')
pxcut = 250
pycut = 150
mxcut = 200
mycut = 150

inDIR  = './FigCube'
outDIR = './FigCube_Cut'
os.system('rm -rf '+outDIR)
os.system('mkdir -p '+outDIR)

inFileNames = os.listdir(inDIR)

nFiles = len(inFileNames)
print 'nFiles = ', nFiles

command1 = 'convert -crop +'+str(pxcut)+'+'+str(pycut)
command2 = 'convert -crop -'+str(mxcut)+'-'+str(mycut)
for ifile in range(nFiles):
   if inFileNames[ifile][-4:] == '.png':
      inFile  = inDIR +'/'+inFileNames[ifile]
      outFile = outDIR+'/'+inFileNames[ifile]
      tmpFile = outDIR+'/'+'tmp.png'
      command = command1 + ' ' + inFile + ' ' + tmpFile
      print command
      os.system(command)
      command = command2 + ' ' + tmpFile + ' ' + outFile
      print command
      os.system(command)
      print ' '

os.system('rm '+tmpFile)
