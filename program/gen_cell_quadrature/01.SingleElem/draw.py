#!/usr/bin/env python
import os
import sys
import time
import matplotlib.pyplot as plt


print '#### into Postprocessing ... ####'

filename = 'Result.dat'
infile = open(filename, 'r')

nn = int(infile.readline())
x  = [0.0 for i in range(nn)]
dx = [0.0 for i in range(nn)]
y  = [0.0 for i in range(nn)]
cx = [0.0 for i in range(nn+1)]
cy = [0.0 for i in range(nn+1)]

for i in range(nn):
  line  = infile.readline()
  x[i]  = float(line)
for i in range(nn+1):
  line = infile.readline()
  cx[i] = float(line)
for i in range(nn):
  dx[i] = cx[i+1]-cx[i]

infile.close()


print nn
print x
print cx

ymin = -0.1
ymax =  0.1


title   = 'Grid (quadrature points, np={0:d}) & its control volume'.format(nn)
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
for i in range(nn):
    # x
    str = 'x{0:d}'.format(i+1)
    axis.text(x[i],ymax*0.1,str,fontsize='9',color='r',horizontalalignment='center',verticalalignment='center')
    # dx
    str = 'dx{0:d}={1:4.2f}'.format(i+1,dx[i])
    axis.text(x[i],ymax*0.2,str,fontsize='9',color='r',horizontalalignment='center',verticalalignment='center')
for i in range(nn+1):
    str = 'x{0:d}/2'.format(2*i+1)
    axis.text(cx[i],ymin*0.2,str,fontsize='9',color='k',horizontalalignment='center',verticalalignment='center')
    axis.axvline(x=cx[i],ymin=0.45,ymax=0.7,c='k',ls='--') 

axis.set_yticks([])
axis.set_yticklabels([])

oufigname = './fig_{0:02d}.png'.format(nn)
pltfig.savefig(oufigname)
plt.show()
