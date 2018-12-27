#!/usr/bin/env python 
from numpy import *
from scipy import *
#from scipy.io import FortranFile
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import cm
import os
import sys
import struct




x = [0.00000,0.00990,0.01980,0.40590,0.98020,1.00000]
y = [0.00000,0.00060,0.00100,0.10000,0.99763,1.00000]

#####################
# plot
#####################
margin = 0.08
wspace = margin*1.3
hspace = 0.25
pltfig = plt.figure(figsize=(12,8))
pltfig.subplots_adjust(left=margin,right=1.0-margin,bottom=1.2*margin,top=1.0-1.2*margin,wspace=wspace,hspace=hspace)

fmts = ['-',':','--','-.']

# ecmwf
axis1 = pltfig.add_subplot(1,1,1)
axis1.set_title('test')
axis1.set_xlabel('nomalized level [a.u.]')
axis1.set_ylabel(r'$\eta$ [a.u.]')
axis1.set_xlim(0,1)
axis1.set_ylim(0,1)
axis1.plot(x,y,'k-')



#pltfig.savefig('test.png')
plt.show()



def polynomial(x, arg):
    w=arg[0]+arg[1]*x+arg[2]*power(x,2)+arg[3]*power(x,3)+arg[4]*power(x,4)
    return w

def l_polynomial(x, arg):
    w=arg[0]*x+arg[1]*power(x,2)+arg[2]*power(x,3)+arg[3]*power(x,4)
    return w

def high_polynomial(x, arg):
    order = len(arg)-1
    w = 0.0
    for i in range(order+1):
        w = w + arg[order-i]*power(x, i)
    return w

def gaussian(x, a, mu, sigma):
    return a*exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma))

def gaussian2(x, a, sigma):
    return a*exp(-(x-1.0)*(x-1.0)/(2.0*sigma*sigma))
