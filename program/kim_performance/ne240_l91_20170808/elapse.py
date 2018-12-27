#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
sys.path.append('/home/jhkim/work/share/python')
from utilities import *


carts00 = ''
comps00 = ['elapse time']
xxx00   = [00.0,45.0,60.0,75.0,90.0,135.0,150.0,165.0,180.0,195.0,210.0,225.0]
yyy00   = [17.0,25.0,30.0,28.0,27.0,038.0,032.0,033.0,025.0,020.0,019.0,017.0]

pltfig,axis = create_plot(ncols=1,nrows=1)
draw_plots(axis,'elapse time',carts00,comps00,xxx00,yyy00, \
           xlabel='time [min]',ylabel='elapse [min]',ymin=0.0,leg_loc=2)
pltfig.savefig('fig.png')
plt.show()

