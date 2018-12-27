#!/usr/bin/env python
import os
import sys
from collections import OrderedDict
sys.path.append('/home/jhkim/work/share/python')
from utilities import *


class class_draw_scalability:



    def __init__(self, axis, title, carts, names, nprocs, walltimes):
        ''' construct '''
        carts = tolist_nd(carts,1)
        names = tolist_nd(names,2)
        ncarts = len(carts)
        if len(carts)!=len(names):
            print('check dimensions of carts and names')
            quit()

        nprocs    = tolist_nd(nprocs,3)
        walltimes = tolist_nd(walltimes,3)
        xmin = 1000000000.0
        xmax = -1000000000.0
        ymin = 1000000000.0
        ymax = -1000000000.0
        if len(nprocs)!=ncarts or len(walltimes)!=ncarts:
            print('check dimensions of nprocs({}) and walltimes({})'.          \
                                            format(len(nprocs),len(walltimes)))
            quit()

        for i in range(ncarts):
            nopts = len(names[i])
            if len(nprocs[i])!=nopts or len(walltimes[i])!=nopts:
                print('check {}th dimensions of nprocs({}) and walltimes({})'. \
                                  format(i+1,len(nprocs[i]),len(walltimes[i])))
                quit()
            for j in range(nopts):
                if len(nprocs[i][j])!=len(walltimes[i][j]):
                    print('check ({},{})th dimension of nprocs({}) and walltimes({})'.\
                                   format(i,j,len(nprocs[i][j]),len(walltimes[i][j])))
                xmin = min(xmin,min(nprocs[i][j]))
                xmax = max(xmax,max(nprocs[i][j]))
                ymin = min(ymin,min(walltimes[i][j]))
                ymax = max(ymax,max(walltimes[i][j]))
 
        halo = xmax*0.05
        xmin = xmin-halo
        xmax = xmax+halo
        halo = ymax*0.05
        ymin = ymin-halo
        ymax = ymax+halo

        fmt_cart = ['k','r','b','g','c']
        fmt_name = ['o-','o--','o:']

        xlabel = 'number of cores [#]'
        ylabel = 'wall-clock times [s]'

        axis_init(axis,title,xlabel,ylabel,xmin,xmax,ymin,ymax)
        for i in range(ncarts):
           for j in range(len(names[i])):
               label = carts[i]+': '+names[i][j]
               axis.plot(nprocs[i][j],walltimes[i][j],fmt_cart[i]+fmt_name[j],label=label)
        y_labels = axis.get_yticklabels()
        for label in y_labels:
            label.set_rotation(90)
        axis.grid(True)
        axis.legend(loc=1,fontsize=11,shadow=True)



