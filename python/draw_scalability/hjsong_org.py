#!/usr/bin/env python

#from numpy import *
#from scipy import *
import numpy as np
import scipy as sp
import math  as mt
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as ptl

ptlfig = ptl.figure(figsize=(8, 6))
fig    = ptl.gcf()

def speedup(n, p):
  w=1.0/(p+(1.0+p)/n)
  return w

def Rfunc(r, R0, L):
  return R0*(1.0+r/L)*np.exp(-r/L)


ndata         = 12
x             = np.array([100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0])
y             = np.array([0.35,   0.33,  0.31,  0.29,  0.27,  0.25,  0.20,  0.16,  0.12,   0.08,   0.04,  0.001])

# test
R00 = 0.38
L0  = 200.0
y   = Rfunc(x, R00, L0)
# test

fun_y = np.array([0.0 for i in range(ndata)])

###   Fitting   ###
popt, pcov = curve_fit(Rfunc, x, y)
print popt, pcov
fun_y = Rfunc(x, popt[0], popt[1])

#print fun_y

###   Plotting   ###
xmax       = 1210.0
ymax       = 1.0
#ymax       = max(y)*1.1

ptl.subplot(1,1,1)
ptl.title('title', size=14)
ptl.xlabel('Distance [km]')
ptl.ylabel('R(r)')

leg1, = ptl.plot(x, y, 'k.')
leg2, = ptl.plot(x, fun_y, 'r-')
ptl.legend([leg1, leg2], ['points', 'function'], shadow=True, loc=2)
ptl.axis([0, xmax, 0, ymax])
ptl.grid(True)


outfilename='test.png'
fig.savefig(outfilename)
ptl.show()
