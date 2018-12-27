#!/usr/bin/env python

from numpy import *
from scipy import *
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as ptl
#import matplotlib.numerix as nu

ptlfig = ptl.figure(figsize=(8, 6))
fig    = ptl.gcf()

def speedup(n, p):
        w=1.0/(p+(1.0+p)/n)
        return w

maxxval       = 260
maxyval       = 260
nres          = 'ne30np4'
ndata         = 13
nproc         = array([1, 4, 8, 12, 16, 20, 24, 32, 64, 80, 128, 192, 256])

avgAll    = array([ 0.0 for i in arange(ndata) ])
stdAll    = array([ 0.0 for i in arange(ndata) ])

infilename = 'data_'+nres+'/result.txt'
infile = open(infilename, 'r')
for idata in arange(ndata):
   line = infile.readline()
   splitline = str(line).split()
   avgAll[idata] = float(splitline[0])
   stdAll[idata] = float(splitline[1])
infile.close()


avgAll0        = avgAll[0]

speedAll       = avgAll0/avgAll
speedAll_std   = speedAll*stdAll


iproc         = array([i for i in arange(1., maxxval+1, 10)])
ideal         = iproc


###   Fitting   ###
popt, pcov = curve_fit(speedup, nproc, speedAll)
speed = speedup(iproc, popt)
eff   = speed/iproc*100.0


###   Plotting   ###

# Plot (1, 1)
ptl.subplot(1,1,1)
ptl.title('Scalability of KIAPS-GM Testbed', size=14)
ptl.xlabel('Number of Processes [#]')
ptl.ylabel('Speed-up')

leg1, = ptl.plot(iproc, ideal, 'k-')
leg2 = ptl.errorbar(nproc, speedAll, yerr=speedAll_std, fmt='rs')
ptl.legend([leg1, leg2], ['Ideal Case', 'Experiment'], shadow=True, loc=2)
ptl.axis([0, maxxval, 0, maxyval])
ptl.grid(True)


outfilename='KIAPSGM_Scalability.png'
fig.savefig(outfilename)
#ptl.show()
