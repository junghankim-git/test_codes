#!/usr/bin/env python
import os
import sys
import time


print '#### into Postprocessing ... ####'

filename = 'Result.dat'
infile = open(filename, 'r')

np   = int(infile.readline())
ne   = int(infile.readline())
nUPs = int(infile.readline())
x  = [0.0 for i in range(nUPs)]
dx = [0.0 for i in range(nUPs)]
y  = [0.0 for i in range(nUPs)]
cx = [0.0 for i in range(nUPs+1)]
cy = [0.0 for i in range(nUPs+1)]

for i in range(nUPs):
  line  = infile.readline()
  x[i]  = float(line)
for i in range(nUPs+1):
  line = infile.readline()
  cx[i] = float(line)
for i in range(nUPs):
  dx[i] = cx[i+1]-cx[i]

infile.close()


print nUPs
print x
print cx
mincx = min(cx)
maxcx = max(cx)
lencx = maxcx-mincx

ymin = -0.1
ymax =  0.1


title   = 'Grid (quadrature points, nUPs={0:d}) & its control volume'.format(nUPs)
pltfig = plt.figure(figsize=(12,4))
pltfig.suptitle('',fontsize=18)
axis = pltfig.add_subplot(1,1,1)
axis.clear()
axis.set_title(title)
axis.set_xlabel('x')
axis.set_ylabel('')
axis.set_xlim(min(cx)+0.1*min(cx),max(cx)+0.1*max(cx))
axis.set_ylim(ymin,ymax)
axis.grid(True)

axis.plot(x,y,'r*',label='grid points')
axis.axhline(y=0.0,xmin=0,xmax=1,c='k',ls='-')
for i in range(nUPs):
  # x
  str = 'x{0:d}'.format(i+1)
  # dx
  str = 'dx{0:d}={1:4.2f}'.format(i+1,dx[i])
for i in range(nUPs+1):
  str = 'x{0:d}/2'.format(2*i+1)
  axis.axvline(x=cx[i],ymin=0.45,ymax=0.7,c='k',ls='--')

axis.set_yticks([])
axis.set_yticklabels([])

oufigname = './fig_{0:02d}np_{1:02}ne.png'.format(np,ne)
pltfig.savefig(oufigname)
plt.show()
