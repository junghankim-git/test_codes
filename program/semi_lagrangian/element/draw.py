#!/usr/bin/env python
import os
import sys
import time
import math
import numpy
import matplotlib.pyplot as plt
sys.path.append('/home/jhkim/work/share/python')
from class_statistics import *


print '#### into Postprocessing ... ####'

filename   = 'Result.dat'
sefilename = '/home/jhkim/TestBed/KSM/1.02/01.dev/Exp_1D_Advec/Result.txt'
if not os.path.isfile(filename):
  print 'check file ...', filename
  quit()

infile = open(filename, 'r')

#=======================================
# Read meta data
#=======================================
np       = int(infile.readline())
ne       = int(infile.readline())
nUPs     = int(infile.readline())
nstep    = int(infile.readline())
dt       = float(infile.readline())
velocity = float(infile.readline())
iexp     = int(infile.readline())
amp      = float(infile.readline())
mu       = float(infile.readline())
sigma    = float(infile.readline())


#=======================================
# Spectral-Element
#=======================================
if os.path.isfile(sefilename):
  print 'checking spectral-element...'
  sefile = open(sefilename, 'r')
  ndim_se  = int(sefile.readline())
  np_se    = int(sefile.readline())
  ne_se    = int(sefile.readline())
  nstep_se = int(sefile.readline())
  dt_se    = float(sefile.readline())
  xmin_se  = float(sefile.readline())
  xmax_se  = float(sefile.readline())
  v_se     = float(sefile.readline())
  amp_se   = float(sefile.readline())
  mu_se    = float(sefile.readline())
  sigma_se = float(sefile.readline())
  line = sefile.readline()  # uptype
  ofreq_se = float(sefile.readline())
  nostep_se= int(sefile.readline())
  line = sefile.readline()  # x

  if ndim_se != 1 or np_se != np or ne_se != ne or nstep_se != nstep or dt_se != dt or v_se != velocity or amp_se != amp or mu_se != mu or sigma_se != sigma:
    print 'check spectral-element setting ... '
    sefile.close()
    nvars = 1
  else:
    nvars = 2
else:
  nvars = 1


#=======================================
# Define variables
#=======================================
x      = [0.0 for i in range(nUPs)]
y      = [0.0 for i in range(nUPs)]
dx     = [0.0 for i in range(nUPs)]
cx     = [0.0 for i in range(nUPs+1)]
# out (Semi-Lagrangian)
T      = [[[0.0 for i in range(nUPs)] for j in range(nstep+1)] for ivar in range(nvars)]
# theory
T_t    = [0.0 for i in range(nUPs)]
mass   = [[0.0 for j in range(nstep+1)] for ivar in range(nvars)]
mass_m = [[0.0 for j in range(nstep+1)] for ivar in range(nvars)]


#=======================================
# Read Lines
#=======================================

line   = infile.readline()
slines = line.split()
for i in range(nUPs):
  x[i]  = float(slines[i])

line   = infile.readline()
slines = line.split()
for i in range(nUPs+1):
  cx[i] = float(slines[i])

for i in range(nstep+1):
  line  = infile.readline()
  slines = line.split()
  for j in range(nUPs):
    T[0][i][j]  = float(slines[j])
  line  = infile.readline()
  mass[0][i] = float(line)

infile.close()


# spectral-element
if nvars > 1:
  for i in range(nstep+1):
    line  = sefile.readline()
    slines = line.split()
    for j in range(nUPs):
      T[1][i][j]  = float(slines[j])
  
  sefile.close()


for i in range(nUPs):
  dx[i] = cx[i+1]-cx[i]

for ivar in range(nvars):
  for i in range(nstep+1):
    mass[ivar][i]   = 0.0
    mass_m[ivar][i] = 0.0
    for j in range(nUPs-1):
      mass[ivar][i] = mass[ivar][i] + T[ivar][i][j]*dx[j]
      if T[ivar][i][j]*dx[j] < 0.0:
        mass_m[ivar][i] = mass_m[ivar][i] + T[ivar][i][j]*dx[j]


print ' * np       = ', np
print ' * ne       = ', ne
print ' * nUPs     = ', nUPs
print ' * nstep    = ', nstep
print ' * dt       = ', dt
print ' * velocity = ', velocity
print ' * mass(SL) = ', mass[0][0]
if nvars > 1:
  print ' * mass(SE) = ', mass[1][0]

#for i in range(nUPs):
#  if T[0][0][i] != T[1][0][i]:
#    print T[0][0][i], T[1][0][i], '<='
#  else:
#    print T[0][0][i], T[1][0][i]


def GetExp(iexp_in, u_in, dt_in, nn, x_in, istep_in):

  psi = [0.0 for i in range(nn)]
  l_xmin = x[0]
  l_xmax = x[nn-1]

  if iexp == 1:
    for i in range(nn):
      m_x   = x_in[i]-u_in*dt_in*istep_in
      while m_x < l_xmin:
        m_x = l_xmax - (l_xmin - m_x)
      if m_x > mu-2.0*sigma and m_x < mu+2.0*sigma:
        psi[i] = amp
      else:
        psi[i] = 0.0
  elif iexp == 2:
    for i in range(nn):
      m_x   = x_in[i]-u_in*dt_in*istep_in
      while m_x < l_xmin:
        m_x = l_xmax - (l_xmin - m_x)
      if m_x > mu-5.0*sigma and m_x < mu+5.0*sigma:
        psi[i] = amp*math.exp(-1.0*(m_x-mu)*(m_x-mu)/(2.0*sigma*sigma))
      else:
        psi[i] = 0.0
  return psi




#   CANVAS 1
mincx = min(cx)
maxcx = max(cx)
lencx = maxcx-mincx
xmin  = min(cx)-0.1*lencx
xmax  = max(cx)+0.1*lencx

ymin = -0.1
ymax =  0.1


'''
title   = 'Grid (quadrature points, nUPs={0:d}) & its control volume'.format(nUPs)
pltfig = plt.figure(figsize=(12,4))
pltfig.suptitle('',fontsize=18)
axis = pltfig.add_subplot(1,1,1)
axis.clear()
axis.set_title(title)
axis.set_xlabel('x')
axis.set_ylabel('')
axis.set_xlim(xmin,xmax)
axis.set_ylim(ymin,ymax)
axis.grid(True)
axis.plot(x,y,'r*',label='grid points')
#axis.plot(cx,cy,'r*',label='control volumes')
axis.axhline(y=0,xmin=0,xmax=1,c='k',ls='-')



for i in range(nUPs):
  # x
  str = 'x{0:d}'.format(i+1)

  # dx
  str = 'dx{0:d}={1:4.2f}'.format(i+1,dx[i])

for i in range(nUPs+1):
  str = 'x{0:d}/2'.format(2*i+1)
  axis.axvline(x=cx[i],ymin=0.45,ymax=0.7,c='k',ls='--')

axis.set_yticks([])
axis.set_yticklabels([])

oufigname = 'Figs/grid_{0:02d}np_{1:02}ne.png'.format(np,ne)
pltfig.savefig(oufigname)
'''






#   CANVAS 2
minx = min(x)
maxx = max(x)
lenx = maxx-minx
xmin = minx-0.1*lenx
xmax = maxx+0.1*lenx

miny = min(T[0][0])
maxy = max(T[0][0])
leny = maxy-miny
ymin = miny-0.1*leny
ymax = maxy+0.4*leny



T_t = GetExp(iexp, velocity, dt, nUPs, x, 0)



# Plot (t = 0)
title = 'Semi-Lagrangian (1D)'
pltfig = plt.figure(figsize=(12,8))
pltfig.suptitle('',fontsize=18)
axis = pltfig.add_subplot(1,1,1)
axis.clear()
axis.set_title(title)
axis.set_xlabel('x')
axis.set_ylabel('')
axis.set_xlim(xmin,xmax)
axis.set_ylim(ymin,ymax)
axis.grid(True)



l2_err     = [0.0 for i in range(nvars)]
mass_err   = [0.0 for i in range(nvars)]
mass_m_err = [0.0 for i in range(nvars)]
labels     = ['' for i in range(nvars)]
labels_h   = ['Semi-Lagrangian  ', 'Spectral-Element ']
labels_msg = '[L2 = {0:6.2E}, mass = {1:11.8f}, {2:4.2f}%, mass(-) = {3:11.8f})]'
for ivar in range(nvars):
  l2_err[ivar]  = get_lerror(T_t, T[ivar][0], 'L2')
  labels[ivar] = labels_h[ivar]+labels_msg.format(l2_err[ivar], mass[ivar][0], 0.0, mass_m[ivar][0])

plot_t, = axis.plot(x,T_t,'k--',label='ideal')
plot,   = axis.plot(x,T[0][0],'r*',label=labels[0])
if nvars > 1:
    plot_se, = axis.plot(x,T[1][0],'b*',label=labels[1])



plt.ion()
oufigname = 'Figs/T_{0:02d}np_{1:02d}ne_{2:003d}.png'.format(np,ne,0)
plt.show()


# Plot (t != 0)
for istep in range(nstep+1):
    # data
    plot.set_ydata(T[0][istep])
    if nvars > 1:
        plot_se.set_ydata(T[1][istep])
  
    # theory
    T_t = GetExp(iexp, velocity, dt, nUPs, x, istep)
    plot_t.set_ydata(T_t)
  
  
  
    # plots
    lplots  = []
    if nvars == 1:
        lplots  = [plot_t, plot]
    if nvars == 2:
        lplots = [plot_t, plot, plot_se]
  
    # labels
    llabels = ['ideal']
    for ivar in range(nvars):
        mass_err[ivar] = abs(mass[ivar][0]-mass[ivar][istep])/mass[ivar][0]*100.0
        l2_err[ivar]   = get_lerror(T_t, T[ivar][istep], 'L2')
        labels[ivar]   = labels_h[ivar]+labels_msg.format(l2_err[ivar], mass[ivar][istep], mass_err[ivar], mass_m[ivar][istep])
        llabels.append(labels[ivar])
  
    plt.legend(lplots, llabels)
    oufigname = 'Figs/T_{0:02d}np_{1:02d}ne_{2:003d}.png'.format(np,ne,istep)
    pltfig.savefig(oufigname)
    if istep == nstep:
        plt.show(True)
    else:
        plt.draw()









