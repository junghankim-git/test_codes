#!/usr/bin/env python
import os
import sys
import time
import numpy
import netCDF4 as nc
import matplotlib.pyplot as plt
from optparse import OptionParser
sys.path.append('/home/jhkim/work/share/python')
from utilities import create_plot,axis_init
from class_statistics import get_lerror

opt = OptionParser()
   ### action: 'store', 'store_const', 'append', 'count', 'callback'
opt.add_option('-b', '--bndry',action='store',dest='bndry',default=2,help='bndry(0,1,2*)',type='int')
opt.add_option('-m', '--monotone',action='store',dest='mono',default=0,help='monotone(0*,1)',type='int')
opt.add_option('-c', '--case',action='store',dest='case',default=6,help='case(1~6*)',type='int')
opt.add_option('-f', '--output',action='store_true',dest='output',default=False,help='figure')

(options, args) = opt.parse_args()

bndry  = options.bndry
mono   = options.mono
case   = options.case
output = options.output

def readmeta(filename):
   infile     = nc.Dataset(filename, mode='r')
   n          = len(infile.dimensions['n'])
   m          = len(infile.dimensions['m'])
   l          = len(infile.dimensions['l'])
   use_spline = infile.use_spline
   ismono     = infile.ismono
   order      = infile.order
   bndry      = infile.bndry
   testcase   = infile.testcase
   sx         = infile.variables['sx'][:].tolist()
   scx        = infile.variables['scx'][:].tolist()
   sy         = infile.variables['sy'][:].tolist()
   tx         = infile.variables['tx'][:].tolist()
   ty_ana     = infile.variables['ty_ana'][:].tolist()
   x_ana      = infile.variables['x_ana'][:].tolist()
   y_ana      = infile.variables['y_ana'][:].tolist()
   infile.close()
   print 'spline = ',use_spline
   print 'order  = ',order
   print 'bndry  = ',bndry
   print 'mono   = ',ismono
   print 'case   = ',testcase
   print 'n      = ',n
   print 'm      = ',m
   print 'l      = ',l

   return n, m, l, sx, scx, sy, tx, ty_ana, x_ana, y_ana



def readfile(filename):
   infile = nc.Dataset(filename, mode='r')
   ty     = infile.variables['ty'][:].tolist()
   ty_fun = infile.variables['ty_fun'][:].tolist()
   infile.close()
   return ty, ty_fun


if mono==0 or mono==1:
   nmon = 1
   monos = [mono]
else:
   nmon = 2
   monos = [0,1]

nexp      = 3
exps      = [[0,2], [1,2], [1,4]]
expnames  = ['PPM', 'PSM', 'QSM']
filenames = [['' for j in range(nexp)] for i in range(nmon)]
for i in range(nmon):
   for j in range(nexp):
      filenames[i][j] = 'result/s{}_o{}_b{}_m{}_t{}.nc'.format(exps[j][0],exps[j][1],bndry,monos[i],case)
      print filenames[i][j]

tys     = [[[] for j in range(nexp)] for i in range(nmon)]
ty_funs = [[[] for j in range(nexp)] for i in range(nmon)]
l2_a    = [[0.0 for j in range(nexp)] for i in range(nmon)]
l2_t    = [[0.0 for j in range(nexp)] for i in range(nmon)]


n, m, l, sx, scx, sy, tx, ty_ana, x_ana, y_ana = readmeta(filenames[0][0])
for i in range(nmon):
   for j in range(nexp):
      tys[i][j], ty_funs[i][j] = readfile(filenames[i][j])
      l2_a[i][j] = get_lerror(y_ana,ty_funs[i][j],'L2')
      l2_t[i][j] = get_lerror(ty_ana,tys[i][j],'L2')


#### Plot #####
fx_string    = 'analytic: '+r'$f(x)$'
if case == 1:
   fx_string = fx_string+r'$=2(x-\mu)$'
elif case == 2:
   fx_string = fx_string+r'$=4(x-\mu)^2-1$'
elif case == 3:
   fx_string = fx_string+r'$=4(x-\mu)^3-1$'
elif case == 4:
   fx_string = fx_string+r'$=sin(x)$'
elif case == 5:
   fx_string = fx_string+r'$=exp(-\alpha(x-\mu)^2)$'
elif case == 6:
   fx_string = fx_string+r'$=0,1$'


xmin  = min(x_ana)
xmax  = max(x_ana)
xmins = xmin-(xmax-xmin)*0.1
xmaxs = xmax+(xmax-xmin)*0.1

ymin  = min(y_ana)
ymax  = max(y_ana)
ymins = ymin-(ymax-ymin)*0.7
ymaxs = ymax+(ymax-ymin)*1.1


# PLOTs
pnx = 3
pny = 2

main_title = ''
pltfig,axis = create_plot(pnx,pny,title=main_title)


'''
# Analytic
axis_init(axis[0][0],'Analytic','x','y',xmins,xmaxs,ymins,ymaxs)
axis[0][0].plot(x_ana,y_ana,'k-',label=fx_string)
axis[0][0].legend(loc=2,fontsize=12,shadow=True)
'''


# Functions
for ie in range(nexp):
   label = 'remap func: '+r'${0:6s},error = {1:4.4f} \%$'.format('y(x)',l2_a[0][ie]*100.0)
   axis_init(axis[ie][0],expnames[ie],'x','y',xmins,xmaxs,ymins,ymaxs)
   axis[ie][0].plot(x_ana,y_ana,'k-',label=fx_string)
   axis[ie][0].plot(sx,sy,'k*',label='source: samples')
   for i in range(n):
      if i == 0:
         axis[ie][0].plot(x_ana[i:l/n],ty_funs[0][ie][i:l/n],'r-',label=label)
      else:
         axis[ie][0].plot(x_ana[i*l/n:(i+1)*l/n],ty_funs[0][ie][i*l/n:(i+1)*l/n],'r-')
   for i in range(n+1):
      axis[ie][0].axvline(x=scx[i],ymin=0,ymax=1,c='k',ls='--')
   axis[ie][0].legend(loc=2,fontsize=12,shadow=True)


'''
# Samples
axis_init(axis[0][1],'Samples (source)','x','y',xmins,xmaxs,ymins,ymaxs)
axis[0][1].plot(x_ana,y_ana,'k--',label='analytic')
axis[0][1].plot(sx,sy,'k*',label='samples(source)')
for i in range(n+1):
   axis[0][1].axvline(x=scx[i],ymin=0,ymax=1,c='k',ls='--')
axis[0][1].legend(loc=2,fontsize=12,shadow=True)
'''

# Values
if nmon==1:
   for ie in range(nexp):
      label = r'target: ${0:6s},error = {1:4.4f} \%$'.format('y_i(x)',l2_t[0][ie]*100.0)
      axis_init(axis[ie][1],expnames[ie],'x','y',xmins,xmaxs,ymins,ymaxs)
      axis[ie][1].plot(x_ana,y_ana,'k--',label=fx_string)
      axis[ie][1].plot(tx,tys[0][ie],'r.',label=label)
      axis[ie][1].legend(loc=2,fontsize=12,shadow=True)
   oufilename = 'figs/t{}_m{}_b{}.png'.format(case,mono,bndry)
else:
   for ie in range(nexp):
      label = 'remap func: '+r'${0:6s},error = {1:4.4f} \%$'.format('y(x)',l2_a[1][ie]*100.0)
      axis_init(axis[ie][1],expnames[ie],'x','y',xmins,xmaxs,ymins,ymaxs)
      axis[ie][1].plot(x_ana,y_ana,'k-',label=fx_string)
      axis[ie][1].plot(sx,sy,'k*',label='source: samples')
      for i in range(n):
         if i == 0:
            axis[ie][1].plot(x_ana[i:l/n],ty_funs[1][ie][i:l/n],'r-',label=label)
         else:
            axis[ie][1].plot(x_ana[i*l/n:(i+1)*l/n],ty_funs[1][ie][i*l/n:(i+1)*l/n],'r-')
      for i in range(n+1):
         axis[ie][1].axvline(x=scx[i],ymin=0,ymax=1,c='k',ls='--')
      axis[ie][1].legend(loc=2,fontsize=12,shadow=True)
   oufilename = 'figs/t{}_b{}.png'.format(case,bndry)

if output:
   print 'save file: ',oufilename
   pltfig.savefig(oufilename)
else:
   plt.show()


