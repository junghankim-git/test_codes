#!/usr/bin/env python
import os
import sys
import time
import numpy
import netCDF4 as nc
import matplotlib.pyplot as plt
sys.path.append('/home/jhkim/work/share/python')
from utilities import create_plot,axis_init


print ' '
print ' '
print '#### into Postprocessing ... ####'


filename   = 'result.nc'
infile     = nc.Dataset(filename, mode='r')
nexp       = len(infile.dimensions['nexp'])
nmet       = len(infile.dimensions['nmethod'])
tri        = infile.tridiagonal
print 'nexp          = '+str(nexp)
print 'nmethod       = '+str(nmet)
if tri==1:
    maintitle = 'Tridiagonal matrix'
else:
    maintitle = 'Pentadiagonal matrix'
nsize  = infile.variables['nsize'][:].tolist()
l2err  = infile.variables['l2err'][:].tolist()
elapse = infile.variables['elapse'][:].tolist()
infile.close()

if nmet==2:
  print('{:>6s}: {:>10s} {:>10s}'.format('n','Fast','LAPACK'))
  for i in range(nexp):
      print('{:6d}: {:10.2e} {:10.2e}'.format(nsize[i], elapse[0][i], elapse[1][i]))
else:
  print('{:>6s}: {:>10s} {:>10s} {:>10s}'.format('n','Fast','LAPACK (ge)','LAPACK (gb)'))
  for i in range(nexp):
      print('{:6d}: {:10.2e} {:10.2e} {:10.2e}'.format(nsize[i], elapse[0][i], elapse[1][i], elapse[2][i]))

xmin  = min(nsize)
xmax  = max(nsize)
xmin  = xmin-(xmax-xmin)*0.1
xmax  = xmax+(xmax-xmin)*0.1
ymint = 0.0
#ymaxt = max(elapse[-1])*1.2
ymaxt = 0.012
ymine = 1.0e-16
ymaxe = 1.0e4


pnx = 2
pny = 1


pltfig,axis = create_plot(pnx,pny,title=maintitle)
if tri==1:
    labels = ['Fast algorithm','LAPACK (general)','LAPACK (tridiagonal)']
else:
    labels = ['Fast algorithm','LAPACK (general)','LAPACK (band)']
forms  = ['k*-','ro-','go-']

# Elapse time
axis_init(axis[0][0],'Elapse time','n [#]','elapse time [s]',xmin,xmax,ymint,ymaxt)
for im in range(nmet):
    axis[0][0].plot(nsize,elapse[im],forms[im],label=labels[im])
axis[0][0].legend(loc=2,fontsize=12,shadow=True)

# L2 error
axis_init(axis[1][0],'L2 error','n [#]','L2 error',xmin,xmax,ymine,ymaxe)
for im in range(nmet):
    axis[1][0].plot(nsize,l2err[im],forms[im],label=labels[im])
axis[1][0].set_yscale('log',nonposx='clip')
axis[1][0].legend(loc=1,fontsize=12,shadow=True)


oufilename = 'post.png'
pltfig.savefig(oufilename)
plt.show()



