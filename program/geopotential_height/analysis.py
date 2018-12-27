#!/usr/bin/env python
import sys
sys.path.append('/home/jhkim/work/share/python')
from utilities import *


nlev = 50

pre = [0.0 for i in range(nlev)]
h_p = [0.0 for i in range(nlev)]
h_e = [0.0 for i in range(nlev)]

file = open('result.txt','r')
k=0
for line in file:
    sline = line.split()
    pre[k] = float(sline[0])
    h_p[k] = float(sline[1])
    h_e[k] = float(sline[2])
    print pre[k], h_p[k], h_e[k]  
    k = k+1
file.close()


pltfig,axis = create_plot(base=7)
#axis_init(axis,'test','pressure','height',-100.,100000.,0.,80000.)

draw_plots(axis,'TT','',['p','exner'],pre,[h_p,h_e], \
           xlabel='pressure',ylabel='height',xmin=-100.,xmax=100000.,ymin=0.,ymax=80000.,leg_loc=3)
#axis.set_yscale('log')

plt.show()
