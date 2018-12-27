#!/usr/bin/env python
import os
import sys
import time
import netCDF4 as nc
import matplotlib.pyplot as plt
sys.path.append('/home/jhkim/work/share/python')
from class_statistics import *
from utilities import create_plot, axis_init


print '#### into Postprocessing ... ####'


filename   = 'result/result.nc'
infile     = nc.Dataset(filename, mode='r')
nfun       = len(infile.dimensions['nfun'])
nder       = len(infile.dimensions['norder'])
n          = len(infile.dimensions['n'])
m          = len(infile.dimensions['m'])
testcase   = infile.testcase
isgll      = infile.isgll
sx     = infile.variables['sx'][:].tolist()
sy     = infile.variables['sy'][:].tolist()
tx     = infile.variables['tx'][:].tolist()
ty_a = [[] for j in range(nder)]
ty_a[0] = infile.variables['ty_a'][:].tolist()
ty_a[1] = infile.variables['typ_a'][:].tolist()
ty_a[2] = infile.variables['typp_a'][:].tolist()
ty_a[3] = infile.variables['tint_a'][:].tolist()
ty   = [[[] for j in range(nder)] for k in range(nfun)]
ty_0 = infile.variables['ty'][:].tolist()
ty_1 = infile.variables['typ'][:].tolist()
ty_2 = infile.variables['typp'][:].tolist()
ty_3 = infile.variables['tint'][:].tolist()
for i in range(nfun):
   ty[i][0][:] = ty_0[i][:]
   ty[i][1][:] = ty_1[i][:]
   ty[i][2][:] = ty_2[i][:]
   ty[i][3][:] = ty_3[i][:]
infile.close()



ty_e = [[0.0 for i in range(nder)] for k in range(nfun)]
for i in range(nfun):
   for j in range(nder):
      ty_e[i][j] = get_lerror(ty_a[j],ty[i][j],'L2')
ty_e[0][2] = 9.9999




labels_poly = [['' for i in range(nder)] for k in range(nfun)]
for i in range(nfun):
  labels_poly[i][0] = '{0:6s}, error = {1:4.4f} %'.format('y(x)', ty_e[i][0]*100.0)
  labels_poly[i][1] = '{0:6s}, error = {1:4.4f} %'.format('y\'(x)', ty_e[i][1]*100.0)
  labels_poly[i][2] = '{0:6s}, error = {1:4.4f} %'.format('y\'\'(x)', ty_e[i][2]*100.0)
  labels_poly[i][3] = '{0:6s}, error = {1:4.4f} %'.format('\int{y(x)}', ty_e[i][3]*100.0)
  if i == 0:
    labels_poly[i][2] = '{0:6s}, error = None'.format('y\'\'(x)')
    labels_poly[i][3] = '{0:6s}, error = None'.format('y\'\'(x)')


#### Plot #####
fx_string    = r'$f(x)$'
fpx_string   = r'$f^{\prime}(x)$'
fppx_string  = r'$f^{\prime\prime}(x)$'
intfx_string = r'$\int{f(x)}$'
if testcase == 1:
   fx_string    = fx_string+r'$=2(x-5)$'
   fpx_string   = fpx_string+r'$=2$'
   fppx_string  = fppx_string+r'$=0$'
   intfx_string = intfx_string+r'$=(x-5)^2-25$'
elif testcase == 2:
   fx_string    = fx_string+r'$=(x-5)^2-1$'
   fpx_string   = fpx_string+r'$=2(x-5)$'
   fppx_string  = fppx_string+r'$=2$'
   intfx_string = intfx_string+r'$=\frac{1}{3}(x-5)^3-\frac{125}{3}$'
elif testcase == 3:
   fx_string    = fx_string+r'$=(x-5)^3-1$'
   fpx_string   = fpx_string+r'$=3(x-5)^2$'
   fppx_string  = fppx_string+r'$=6(x-5)$'
   intfx_string = intfx_string+r'$=\frac{1}{4}(x-5)^4-\frac{625}{4}$'
elif testcase == 4:
   fx_string    = fx_string+r'$=sin(x)$'
   fpx_string   = fpx_string+r'$=cos(x)$'
   fppx_string  = fppx_string+r'$=-sin(x)$'
   intfx_string = intfx_string+r'$=-cos(x)$'
elif testcase == 5:
   fx_string    = fx_string+r'$=exp(-(x-5)^2)$'
   fpx_string   = fpx_string+r'$=-2(x-5)exp(-(x-5)^2)$'
   fppx_string  = fppx_string+r'$=2((x-5)^2-1)exp(-(x-5)^2)$'

xmin  = min(tx)
xmax  = max(tx)
xmins = xmin-(xmax-xmin)*0.1
xmaxs = xmax+(xmax-xmin)*0.1

ymin  = min(ty_a[0])
ymax  = max(ty_a[0])
ymins = ymin-(ymax-ymin)*0.7
ymaxs = ymax+(ymax-ymin)*0.9
titles = ['Lagrangian','Linear Spline','Quadratic Spline','Cubic Spline']


nx = 3
ny = 2
nplots = nx*ny


pltfig, axis = create_plot(nx,ny,title='Interpolation')

# Analytic
axis_init(axis[0][0],'Analytic','x','y',xmins,xmaxs,ymins,ymaxs)
axis[0][0].plot(tx,ty_a[0],'k-', label=fx_string)
axis[0][0].plot(tx,ty_a[1],'k--',label=fpx_string)
axis[0][0].plot(tx,ty_a[2],'k-.',label=fppx_string)
axis[0][0].plot(tx,ty_a[3],'k:',label=intfx_string)
axis[0][0].legend(loc=2,fontsize=12,shadow=True)

# Samples
axis_init(axis[1][0],'Samples','x','y',xmins,xmaxs,ymins,ymaxs)
axis[1][0].plot(tx,ty_a[0],'k-', label='analytic')
axis[1][0].plot(sx,sy,     'k*', label='samples')
axis[1][0].legend(loc=2,fontsize=12,shadow=True)

# Lagrangian
axis_init(axis[2][0],titles[0],'x','y',xmins,xmaxs,ymins,ymaxs)
axis[2][0].plot(sx,sy,      'k*',label='samples')
axis[2][0].plot(tx,ty[0][0],'r-',label=labels_poly[0][0])
axis[2][0].plot(tx,ty[0][1],'b-',label=labels_poly[0][1])
axis[2][0].plot(tx,ty[0][2],'g-',label=labels_poly[0][2])
axis[2][0].plot(tx,ty[0][3],'c-',label=labels_poly[0][3])
axis[2][0].legend(loc=2,fontsize=12,shadow=True)

# Splines
for ix in range(nx):
    axis_init(axis[ix][1],titles[ix+1],'x','y',xmins,xmaxs,ymins,ymaxs)
    axis[ix][1].plot(sx,sy,         'k*',label='samples')
    axis[ix][1].plot(tx,ty[ix+1][0],'r-',label=labels_poly[ix+1][0])
    axis[ix][1].plot(tx,ty[ix+1][1],'b-',label=labels_poly[ix+1][1])
    axis[ix][1].plot(tx,ty[ix+1][2],'g-',label=labels_poly[ix+1][2])
    axis[ix][1].plot(tx,ty[ix+1][3],'c-',label=labels_poly[ix+1][3])
    axis[ix][1].legend(loc=2,fontsize=12,shadow=True)

oufilename = 'testcase'+str(testcase)+'.png'
pltfig.savefig(oufilename)
plt.show()





