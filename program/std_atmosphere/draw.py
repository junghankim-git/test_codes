#!/usr/bin/env python
import os
import sys
import time
from numpy import *
from scipy import *
import matplotlib.pyplot as plt
sys.path.append('/home/jhkim/work/share/python')
from utilities import create_plot, axis_init



nlayers = 100
nlevs   = 137+1

H  = [0.0 for i in range(nlayers)]
T  = [0.0 for i in range(nlayers)]

P  = [0.0 for i in range(nlevs)]
H1 = [0.0 for i in range(nlevs)]


filename = 'H-T.ascii'
infile = open(filename, 'r')

k = 0
for iline in infile:
   line = iline.strip()
   splitline = str(line).split()
   H[k] = float(splitline[0])/1000.0
   T[k] = float(splitline[1])
   k = k + 1
infile.close()

print k, nlayers
if k != nlayers:
   print 'check nlayer'
   quit()


filename = 'P-H.ascii'
infile = open(filename, 'r')

k = 0
for iline in infile:
   line = iline.strip()
   splitline = str(line).split()
   P[k]  = float(splitline[0])/100.0
   H1[k] = float(splitline[1])/1000.0
   k = k + 1
infile.close()

print k, nlevs
if k != nlevs:
   print 'check nlayer'
   quit()

######################################
## Plotting
######################################

nrows = 1
ncols = 2


nAtms   = 8
Atms    = ['Troposphere', 'Tropopause', 'Stratosphere', 'Stratosphere', 'Stratopause', 'Mesosphere', 'Mesosphere', 'Mesopause']
AtmLay  = [0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852]
AtmCol  = ['b', 'g', 'g', 'g', 'r', 'r', 'r', 'c']
AtmPosX = [165, 2.0e-2]
AtmPosH = [0.0 for i in range(nAtms)]
for i in range(nAtms):
   if i != nAtms-1:
      AtmPosH[i] = 0.5*(AtmLay[i]+AtmLay[i+1])-1.0
   else:
      AtmPosH[i] = AtmLay[i]+1.2


xmins   = [[160.0],[1.0e-2]]
xmaxs   = [[300.0],[1.0e3]]
ymins   = [[0.0],[0.0]]
ymaxs   = [[90.0],[90.0]]
titles  = [['Temperature-Height'],['Pressure-Height']]
xlabels = [['Temperature [K]'],['Pressure [hPa]']]
ylabels = [['Height [km]'],['Height [km]']]

pltfig, axis = create_plot(ncols,nrows,title='Standard Atmosphere',yscale=2.0)

for j in range(nrows):
    for i in range(ncols):
      axis_init(axis[i][j],titles[i][j],xlabels[i][j],ylabels[i][j],xmins[i][j],xmaxs[i][j],ymins[i][j],ymaxs[i][j])
      if i == 0:
         axis[i][j].plot(T,H,'k-',lw=2.0)
      if i == 1:
         axis[i][j].plot(P,H1,'k-',lw=2.0)
         axis[i][j].set_xscale('log')
      for ilayer in range(nAtms):
         if ilayer != 0: axis[i][j].axhline(y=AtmLay[ilayer],xmin=0.0,xmax=1.0,c=AtmCol[ilayer],ls='--')
         axis[i][j].text(AtmPosX[j], AtmPosH[ilayer], Atms[ilayer], size=15, color=AtmCol[ilayer], style='italic')


outfilename = './std_atmosphere.png'
pltfig.savefig(outfilename)
plt.show()





