#!/usr/bin/env python
import os
import sys
import time
import numpy
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
sys.path.append('/home/jhkim/work/share/python')
from utilities import create_plot,axis_init
from class_statistics import get_lerror


def coord2box(si, ei, sj, ej, color, hatch=False):
 
    verts = [[si,sj], [si,ej], [ei,ej], [ei,sj], [si,sj]]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    path = Path(verts, codes)
    if color!='none':
        lw   =0.0
        alpha=0.7
    else:
        lw   =1.5
        alpha=1.0
    if hatch:
        return patches.PathPatch(path,facecolor=color,lw=lw,hatch='/',alpha=alpha)
    else:
        return patches.PathPatch(path,facecolor=color,lw=lw,alpha=alpha)

print ' '
print ' '
print '#### Draw ####'



n   =  5

sx  = [float(i+1) for i in range(n)]
scx = [float(i)+0.5 for i in range(n+1)]
sy  = [1.0, 2.0, 1.0, 3.0, 2.0]
tx  = [0.0]
tcx = [0.0]
ty  = [0.0]
gi  = 3.2
dgi = 0.6


print sx
print scx

xmin = 0.0
xmax = 6.0
ymin =-0.5
ymax = 5.0
dx   = 0.1
dy   = (ymax-ymin)*0.07

pnx = 1
pny = 1
pltfig,axis = create_plot(pnx,pny,base=6.0,yscale=0.6,title='')
axis_init(axis,'','x','y',xmin,xmax,ymin,ymax)
axis.axhline(y=0.0,xmin=xmin,xmax=xmax,c='k',ls='-')
axis.axvline(x=0.0,ymin=ymin,ymax=ymax,c='k',ls='-')
axis.plot(sx,sy,'ko',label='samples')
for i in range(n):
    axis.text(sx[i]+dx,-0.5*dy,'$x_{:1d}$'.format(i+1), \
              fontsize=12,ha='center',va='center')
    axis.text(sx[i],sy[i]+dy,r'$\bar{{f}}_{0:}={1:}$'.format(i+1,sy[i]), \
              fontsize=12,ha='center',va='center')
    # gi
    axis.text(gi,-0.9*dy,r'$x_{i}^\prime$',color='b',fontsize=12,ha='center',va='center')
    axis.axvline(x=gi,ymin=0,ymax=1,c='b',ls='--')
axis.legend(loc=2,fontsize=12,shadow=True)
pltfig.savefig('figs/present/01.png')

for i in range(n+1):
    axis.text(scx[i]+dx,-0.5*dy,r'$x_\frac{{{}}}{{{}}}$'.format(2*i+1,2),color='r', \
             fontsize=12,ha='center',va='center')
    axis.axvline(x=scx[i],ymin=0,ymax=1,c='r',ls='--')
pltfig.savefig('figs/present/02.png')

for i in range(n):
    patch = coord2box(scx[i],scx[i+1],0.0,sy[i],'g')
    axis.add_patch(patch)
pltfig.savefig('figs/present/03.png')

for i in range(n):
    axis.text(sx[i],sy[i]+2*dy,r'$f_{0:}(x)$'.format(i+1),color='g', \
              fontsize=12,ha='center',va='center')
pltfig.savefig('figs/present/04.png')


axis.text(gi-dgi+dx,-0.9*dy,r'$x_{i-1/2}^\prime$',color='b',fontsize=12,ha='center',va='center')
axis.text(gi+dgi+dx,-0.9*dy,r'$x_{i+1/2}^\prime$',color='b',fontsize=12,ha='center',va='center')
if True:
    axis.axvline(x=gi-dgi,ymin=0,ymax=1,c='k',ls='--')
    axis.axvline(x=gi+dgi,ymin=0,ymax=1,c='k',ls='--')
    patch = coord2box(gi-dgi,gi+dgi,0.0,5.0,'none',hatch=True)
    axis.add_patch(patch)
else:
    patch = coord2box(gi-dgi,scx[3],0.0,sy[2],'b',hatch=True)
    axis.add_patch(patch)
    patch = coord2box(scx[3],gi+dgi,0.0,sy[3],'b',hatch=True)
    axis.add_patch(patch)
pltfig.savefig('figs/present/05.png')


#==================================================================

def gen_fun(nn, xi, xm, ym, yi, mono, ll):
    ou_x = [[0.0 for i in range(ll)] for j in range(nn)]
    ou_y = [[0.0 for i in range(ll)] for j in range(nn)]
    for i in range(nn):
        dx = xi[i+1]-xi[i]
        ou_dx = dx/float(ll)
        for j in range(ll):
            ou_x[i][j] = xi[i]+ou_dx*float(j)
            xxi = ou_x[i][j]-xi[i]/dx
            if not mono:
                aa = yi[i][0]
                bb = -4.0*yi[i][0]-2.0*yi[i][1]+6.0*ym[i]
                cc =  3.0*yi[i][0]+3.0*yi[i][1]-6.0*ym[i]
            else:
                if   i==0:
                  aa = 0.5
                  bb = 0.40669856459330145
                  cc = 0.88995215311004827
                elif i==1:
                  aa = 2.0
                  bb = 0.0
                  cc = 0.0
                elif i==2:
                  aa = 1.0
                  bb = 0.0
                  cc = 0.0
                elif i==3:
                  aa = 3.0
                  bb = 0.0
                  cc = 0.0
                elif i==4:
                  aa = 2.8875598086124405
                  bb = -2.5502392344497622
                  cc = 1.1626794258373216
            ou_y[i][j] = aa+bb*xxi+cc*xxi*xxi
    return ou_x, ou_y


l        = 100
scy      = [[0.0 for i in range(n+1)] for j in range(2)]
mcy      = [[[0.0, 0.0] for i in range(n+1)] for j in range(2)]
ismethod = True
ismono   = True
if ismethod:
    scy[0] = [1.0/3.0,5.0/3.0,1.41667,2.0,2.75,1.25]
    scy[1] = [0.5,1.79665,1.313397,1.94967,2.887559,1.5]
    if not ismono:
        # PPM & PSM
        #scy[0] = [1.0/3.0,5.0/3.0,1.41667,2.0,2.75,1.25]
        mcy[0][0] = [1.0/3.0,5.0/3.0]
        mcy[0][1] = [5.0/3.0,1.41667]
        mcy[0][2] = [1.41667,2.0]
        mcy[0][3] = [2.0,2.75]
        mcy[0][4] = [2.75,1.25]
        #scy[1] = [0.5,1.79665,1.313397,1.94967,2.887559,1.5]
        mcy[1][0] = [0.5,1.79665]
        mcy[1][1] = [1.79665,1.313397]
        mcy[1][2] = [1.313397,1.94967]
        mcy[1][3] = [1.94967,2.887559]
        mcy[1][4] = [2.887559,1.5]
    else:
        # PPM & PSM (monotone)
        #scy[0] = [1.0/3.0,5.0/3.0,2.0,1.0,3.0,1.33333333]
        mcy[0][0] = [1.0/3.0,5.0/3.0]
        mcy[0][1] = [2.0,2.0]
        mcy[0][2] = [1.0,1.0]
        mcy[0][3] = [3.0,3.0]
        mcy[0][4] = [2.75,1.25]
        #mcy[0][4] = [2.66667,1.33333333]
        #scy[1] = [0.5,1.79665,1.313397,1.94967,2.887559,1.5]
        mcy[1][0] = [0.5,1.79665]
        mcy[1][1] = [1.79665,1.313397]
        mcy[1][2] = [1.313397,1.94967]
        mcy[1][3] = [1.94967,2.887559]
        mcy[1][4] = [2.887559,1.5]
else:
    # no constraint
    scy[0] = [0.5,2.0,2.5,3.5,1.0,2.0]
    mcy[0][0] = [0.5,2.0]
    mcy[0][1] = [2.0,2.5]
    mcy[0][2] = [2.5,3.5]
    mcy[0][3] = [3.5,1.0]
    mcy[0][4] = [1.0,2.0]
    scy[1] = [2.0,0.5,1.5,1.5,2.5,3.0]
    mcy[1][0] = [2.0,0.5]
    mcy[1][1] = [0.5,1.5]
    mcy[1][2] = [1.5,1.5]
    mcy[1][3] = [1.5,2.5]
    mcy[1][4] = [2.5,3.0]
ttx    = [[] for j in range(2)]
tty    = [[] for j in range(2)]
for i in range(2):
    if i==0:
        ttx[i],tty[i] = gen_fun(n,scx,sx,sy,mcy[i],False,l)
    else:
        ttx[i],tty[i] = gen_fun(n,scx,sx,sy,mcy[i],ismono,l)
pnx  = 2
pny  = 1
pltfig2,axis2 = create_plot(pnx,pny,base=6.0,yscale=0.6,title='')
if ismethod:
    if ismono:
        titles = ['PPM (monotone)','PSM (monotone)']
    else:
        titles = ['PPM','PSM']
else:
    titles = ['','']
for j in range(2):
    axis_init(axis2[j][0],titles[j],'x','y',xmin,xmax,ymin,ymax)
    axis2[j][0].axhline(y=0.0,xmin=xmin,xmax=xmax,c='k',ls='-')
    axis2[j][0].axvline(x=0.0,ymin=ymin,ymax=ymax,c='k',ls='-')
    axis2[j][0].plot(sx,sy,'ko',label='samples')
    for i in range(n):
        axis2[j][0].text(sx[i]+dx,-0.5*dy,'$x_{:1d}$'.format(i+1), \
                  fontsize=12,ha='center',va='center')
        axis2[j][0].text(sx[i],sy[i]+dy,r'$\bar{{f}}_{0:}={1:}$'.format(i+1,sy[i]), \
                  fontsize=12,ha='center',va='center')
    
    for i in range(n+1):
        axis2[j][0].text(scx[i]+dx,-0.5*dy,r'$x_\frac{{{}}}{{{}}}$'.format(2*i+1,2),color='r', \
                 fontsize=12,ha='center',va='center')
        axis2[j][0].axvline(x=scx[i],ymin=0,ymax=1,c='r',ls='--')

    # my cy
    for i in range(n+1):
        axis2[j][0].plot(scx,scy[j],'ro')
        axis2[j][0].text(scx[i],scy[j][i]+dy,r'$f_\frac{{{}}}{{{}}}={:5.2f}$'.format(2*i+1,2,scy[j][i]), \
                    color='r',fontsize=12,ha='center',va='center')
    for i in range(n):
        if i==0:
            axis2[j][0].plot(ttx[j][i],tty[j][i],'g-',label='remap function')
        else:
            axis2[j][0].plot(ttx[j][i],tty[j][i],'g-')
    axis2[j][0].legend(loc=2,fontsize=12,shadow=True)
pltfig2.savefig('figs/present/06.png')


#==================================================================









plt.show()

