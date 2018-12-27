#!/usr/bin/env python
import os
import sys
import time
from numpy import *
from scipy import *
import matplotlib.pyplot as plt



filename = 'Result.dat'
infile = open(filename, 'r')

n        = int(infile.readline())
nout     = int(infile.readline())
nwave    = int(infile.readline())
testcase = int(infile.readline())
message  = infile.readline()
print message

print 'n    = ' + str(n)

####  Initialization variables  ####
eta_i       = array([0.0 for i in arange(n)])
psi_i       = array([0.0 for i in arange(n)])
reta_i      = array([0.0 for i in arange(nout)])
rpsi_i      = array([0.0 for i in arange(nout)])
apsi_i      = array([0.0 for i in arange(nout)])
diff_i      = array([0.0 for i in arange(nout)])

for i in arange(n):
  line      = infile.readline()
  splitline = str(line).split()
  eta_i[i]  = float(splitline[0])
  psi_i[i]  = float(splitline[1])
message= infile.readline()
print message

for i in arange(nout):
  line      = infile.readline()
  splitline = str(line).split()
  reta_i[i] = float(splitline[0])
  rpsi_i[i] = float(splitline[1])
message= infile.readline()
print message

# analytic solution
for i in arange(nout):
  line      = infile.readline()
  splitline = str(line).split()
  reta_i[i] = float(splitline[0])
  apsi_i[i] = float(splitline[1])
message= infile.readline()
print message

infile.close()

for i in arange(nout):
  diff_i[i] = rpsi_i[i] - apsi_i[i] 


######################################
## Plotting
######################################

nrows = 1
ncols = 3

mmin    = min(apsi_i)
mmax    = max(apsi_i)
ext     = 0.1*(mmax-mmin)

xmins   = [[-0.1, -0.1, -0.1] for i in range(nrows)] 
xmaxs   = [[1.1, 1.1, 1.1] for i in range(nrows)] 
ymins   = [[mmin-ext, mmin-ext, -0.0001] for i in range(nrows)] 
ymaxs   = [[mmax+ext, mmax+ext, 0.0001] for i in range(nrows)] 
titles  = [['Input Points', 'Interpolation', 'Difference'] for i in range(nrows)]
etat    = r'$\eta$'
etat    = r'$\eta$'
xlabels = [[etat, etat, etat] for i in range(nrows)]
ylabels = [[r'$\psi$', r'$\psi^\prime$', r'$\psi^{\prime}-\psi$'] for i in range(nrows)]


def axis_init(axis_in,xmin,xmax,ymin,ymax,title,xlabel,ylabel):
    axis_in.clear()
    axis_in.set_title(title)
    axis_in.set_xlabel(xlabel)
    axis_in.set_ylabel(ylabel)
    axis_in.set_xlim(xmin,xmax)
    axis_in.set_ylim(ymin,ymax)
    axis_in.grid(True)

pltfig = plt.figure(figsize=(18,5))
pltfig.suptitle('Cubic Spline Interpolation',fontsize=18)
axis = [[pltfig.add_subplot(nrows,ncols,j+1) for j in range(ncols)] for i in range(nrows)]


for i in range(nrows):
   for j in range(ncols):
      axis_init(axis[i][j],xmins[i][j],xmaxs[i][j],ymins[i][j],ymaxs[i][j],titles[i][j],xlabels[i][j],ylabels[i][j])
      if j == 0:
         axis[i][j].plot(reta_i,apsi_i,'--')
         axis[i][j].plot(eta_i,psi_i,'ro')
      if j == 1:
         axis[i][j].plot(reta_i,apsi_i,'--')
         axis[i][j].plot(reta_i,rpsi_i,'r--')
      if j == 2:
         axis[i][j].plot(reta_i,diff_i,'r--')

#      if j == 0:
#         canvas1.Plot(i, j, eta_i, psi_i, 'ro')
#         canvas1.Plot(i, j, reta_i, rpsi_i, 'r--')
outfilename = 'Figs/CubicSpline_'+str(testcase)+'.png'
pltfig.savefig(outfilename)
plt.show()






