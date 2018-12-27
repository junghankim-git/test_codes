#!/usr/bin/env python

from numpy import *
from scipy import *
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.numerix as nu


def speedup(nproc, f):
        speed = 1.0/(f+(1.0+f)/nproc)
        return speed

maxnproc = 4000.0
ndata    = 100
d_nproc  = maxnproc/real(ndata)

frac1 = 0.01
frac2 = 0.001
frac3 = 0.0001

nproc = [int(d_nproc*(i+1)) for i in range(ndata)]
exp1  = [speedup(nproc[i],frac1) for i in range(ndata)]
exp2  = [speedup(nproc[i],frac2) for i in range(ndata)]
exp3  = [speedup(nproc[i],frac3) for i in range(ndata)]


###   Plotting   ###
pltfig = plt.figure(figsize=(8,6))
axis   = pltfig.add_subplot(1,1,0)
axis.set_xlim(0,maxnproc*1.01)
axis.set_ylim(0,maxnproc*1.01)
axis.set_title('Scalability (Amdahl\'s Law)',size=14)
axis.set_xlabel('Number of Processes [#]')
axis.set_ylabel('Speed-up')

id  = axis.plot(nproc,nproc,'k--',lw=2,label='ideal')
ex1 = axis.plot(nproc,exp1,'g-',lw=2,label='frac:{}'.format(frac1))
ex2 = axis.plot(nproc,exp2,'b-',lw=2,label='frac:{}'.format(frac2))
ex3 = axis.plot(nproc,exp3,'r-',lw=2,label='frac:{}'.format(frac3))
axis.legend(shadow=True, loc=2)
axis.grid(True)

pltfig.savefig('ideal.png')
plt.show()
