#!/usr/bin/env python
import os
import sys
import time
import numpy
import netCDF4 as nc
import matplotlib.pyplot as plt
sys.path.append('../../..//share/python')
from utilities import create_plot,axis_init
from class_statistics import get_lerror


print ' '
print ' '
print '#### into Postprocessing ... ####'


filename   = 'result/result.nc'
infile     = nc.Dataset(filename, mode='r')
n          = len(infile.dimensions['n'])
m          = len(infile.dimensions['m'])
l          = len(infile.dimensions['l'])
use_spline = infile.use_spline
ismono     = infile.ismono
order      = infile.order
bndry      = infile.bndry
testcase   = infile.testcase
print 'use_spline    = '+str(use_spline)
print 'order         = '+str(order)
print 'bndry         = '+str(bndry)
print 'monotonic     = '+str(ismono)
print 'testcase      = '+str(testcase)
print 'n             = '+str(n)
print 'l             = '+str(l)
print 'm             = '+str(m)
sx     = infile.variables['sx'][:].tolist()
scx    = infile.variables['scx'][:].tolist()
sy     = infile.variables['sy'][:].tolist()
siy    = infile.variables['siy'][:].tolist()
tx     = infile.variables['tx'][:].tolist()
tcx    = infile.variables['tcx'][:].tolist()
ty     = infile.variables['ty'][:].tolist()
tiy    = infile.variables['tiy'][:].tolist()
ty_ana = infile.variables['ty_ana'][:].tolist()
x_ana  = infile.variables['x_ana'][:].tolist()
y_ana  = infile.variables['y_ana'][:].tolist()
ty_fun = infile.variables['ty_fun'][:].tolist()
ty_int = infile.variables['ty_int'][:].tolist()
infile.close()

l2_t= get_lerror(ty_ana,ty,'L2')
l2_a= get_lerror(y_ana,ty_fun,'L2')
print 'error = ', l2_t
print 'error = ', l2_a



#### Plot #####
fx_string    = r'$f(x)$'
if testcase == 1:
   fx_string = fx_string+r'$=2(x-\mu)$'
elif testcase == 2:
   fx_string = fx_string+r'$=4(x-\mu)^2-1$'
elif testcase == 3:
   fx_string = fx_string+r'$=4(x-\mu)^3-1$'
elif testcase == 4:
   fx_string = fx_string+r'$=sin(x)$'
elif testcase == 5:
   fx_string = fx_string+r'$=exp(-\alpha(x-\mu)^2)$'
elif testcase == 6:
   fx_string = fx_string+r'$=0,1$'


xmin  = min(x_ana)
xmax  = max(x_ana)
xmins = xmin-(xmax-xmin)*0.1
xmaxs = xmax+(xmax-xmin)*0.1

ymin  = min(y_ana)
ymax  = max(y_ana)
ymins = ymin-(ymax-ymin)*0.7
ymaxs = ymax+(ymax-ymin)*0.9

pnx = 3
pny = 2

if use_spline==1:
   if order == 2:
      MainString = 'PSM (Parabolic Spline Method)'
   elif order == 4:
      MainString = 'QSM (Quartic Spline Method)'
else:
   if order == 0:
      MainString = 'PCoM (Piecewise Constant Method)'
   elif order == 1:
      MainString = 'PLM (Piecewise Linear Method)'
   elif order == 2:
      MainString = 'PPM (Piecewise Parabolic Method)'
   '''
   elif order == 3:
      MainString = 'PCM (Piecewise Cubic Method)'
   elif order == 4:
      MainString = 'QSM (Quartic Spline Method)'
   '''

pltfig,axis = create_plot(pnx,pny,title=MainString)

# Analytic
axis_init(axis[0][0],'Analytic','x','y',xmins,xmaxs,ymins,ymaxs)
axis[0][0].plot(x_ana,y_ana,'k-',label=fx_string)
axis[0][0].legend(loc=2,fontsize=12,shadow=True)

# Samples
axis_init(axis[1][0],'Samples (source)','x','y',xmins,xmaxs,ymins,ymaxs)
axis[1][0].plot(x_ana,y_ana,'k--',label='analytic')
axis[1][0].plot(sx,sy, 'k*',label='samples(source)')
for i in range(n+1):
   axis[1][0].axvline(x=scx[i],ymin=0,ymax=1,c='k',ls='--')
axis[1][0].legend(loc=2,fontsize=12,shadow=True)

# Integral values at cells
axis_init(axis[2][0],'Integral values at each cells','x','y',xmins,xmaxs,ymins,ymaxs)
axis[2][0].plot(sx,sy, 'k*',label='samples(source)')
for i in range(n+1):
   axis[2][0].axvline(x=scx[i],ymin=0,ymax=2,c='k',ls='--')
axis[2][0].plot(scx,siy,'b.',label='integral ('+str(siy[-1])+')')
axis[2][0].legend(loc=2,fontsize=12,shadow=True)

# Spline
label = r'${0:6s},error = {1:4.4f} \%$'.format('y(x)',l2_a*100.0)
axis_init(axis[0][1],'Functions','x','y',xmins,xmaxs,ymins,ymaxs)
axis[0][1].plot(sx,sy,'k*',label='smaples(source)')
axis[0][1].plot(x_ana,y_ana,'k-',label=fx_string)
for i in range(n):
   if i == 0:
      axis[0][1].plot(x_ana[i:l/n],ty_fun[i:l/n],'r-',label=label)
   else:
      axis[0][1].plot(x_ana[i*l/n:(i+1)*l/n],ty_fun[i*l/n:(i+1)*l/n],'r-')
for i in range(n+1):
   axis[0][1].axvline(x=scx[i],ymin=0,ymax=1,c='k',ls='--')
axis[0][1].legend(loc=2,fontsize=12,shadow=True)

# Integral
axis_init(axis[1][1],'Integration of function','x','y',xmins,xmaxs,ymins,ymaxs)
for i in range(n):
   axis[1][1].axvline(x=sx[i],ymin=0,ymax=2,c='r',ls='--')
axis[1][1].plot(sx,sy,'k*',label='samples(source)')
axis[1][1].plot(scx,siy,'b.',label='integral ('+str(siy[-1])+')')
axis[1][1].plot(x_ana,ty_int,'g-',label=r'$\int{y_i(x)}dx$ '+'('+str(ty_int[-1])+')')
axis[1][1].legend(loc=2,fontsize=12,shadow=True)

# Value
label = r'target(${0:6s},error = {1:4.4f} \%$)'.format('y_i(x)',l2_t*100.0)
axis_init(axis[2][1],'Target','x','y',xmins,xmaxs,ymins,ymaxs)
axis[2][1].plot(x_ana,y_ana,'k--',label=fx_string)
axis[2][1].plot(sx,sy,'k*',label='samples')
axis[2][1].plot(tx,ty,'r.',label=label)
axis[2][1].legend(loc=2,fontsize=12,shadow=True)


if use_spline==1:
   oufilename = 'figs/spline'+'_o'+str(order)+'_t'+str(testcase)+'_m'+str(ismono)+'.png'
else:
   oufilename = 'figs/piecewise'+'_o'+str(order)+'_t'+str(testcase)+'_m'+str(ismono)+'.png'
pltfig.savefig(oufilename)
plt.show()



