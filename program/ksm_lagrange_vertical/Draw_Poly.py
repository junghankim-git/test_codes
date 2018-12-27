#!/usr/bin/env python
from numpy import *
from scipy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time
from matplotlib import cm
import os
import sys


filename = 'Result.dat'
infile = open(filename, 'r')

n       = int(infile.readline())
frac    = int(infile.readline())
testcae = int(infile.readline())
nx_i = frac*(n+1)
nx_m = frac*n
message= infile.readline()
print message

print 'n    = ' + str(n)
print 'nx_i = ' + str(nx_i)
print 'nx_m = ' + str(nx_m)

####  Initialization variables  ####
eta_i       = array([0.0 for i in arange(n+1)])
psi_i       = array([0.0 for i in arange(n+1)])
x_i         = array([0.0 for i in arange(nx_i)])
lag_i       = array([[0.0 for i in arange(nx_i)] for j in arange(n+1)])
dlag_i      = array([[0.0 for i in arange(nx_i)] for j in arange(n+1)])
ilag_i      = array([[0.0 for i in arange(n+1)]  for j in arange(n+1)])
wlag_i      = array([[0.0 for i in arange(nx_i)] for j in arange(n+1)])
wdlag_i     = array([[0.0 for i in arange(nx_i)] for j in arange(n+1)])
psi_inter_i = array([0.0 for i in arange(nx_i)])
diff_psi_i  = array([0.0 for i in arange(nx_i)])
int_psi_i   = array([0.0 for i in arange(nx_i)])
#analytic
ana_diff_i  = array([0.0 for i in arange(nx_i)])
ana_int_i   = array([0.0 for i in arange(nx_i)])

eta_m       = array([0.0 for i in arange(n)])
psi_m       = array([0.0 for i in arange(n)])
x_m         = array([0.0 for i in arange(nx_m)])
lag_m       = array([[0.0 for i in arange(nx_m)] for j in arange(n)])
dlag_m      = array([[0.0 for i in arange(nx_m)] for j in arange(n)])
ilag_m      = array([[0.0 for i in arange(n)]    for j in arange(n)])
wlag_m      = array([[0.0 for i in arange(nx_m)] for j in arange(n)])
wdlag_m     = array([[0.0 for i in arange(nx_m)] for j in arange(n)])
psi_inter_m = array([0.0 for i in arange(nx_m)])
diff_psi_m  = array([0.0 for i in arange(nx_m)])
int_psi_m   = array([0.0 for i in arange(nx_m)])
#analytic
ana_diff_m  = array([0.0 for i in arange(nx_m)])
ana_int_m   = array([0.0 for i in arange(nx_m)])

for i in arange(n+1):
  line      = infile.readline()
  splitline = str(line).split()
  eta_i[i]  = float(splitline[0])
  psi_i[i]  = float(splitline[1])

for i in arange(n):
  line      = infile.readline()
  splitline = str(line).split()
  eta_m[i]  = float(splitline[0])
  psi_m[i]  = float(splitline[1])

message= infile.readline()
print message

# half level
for j in arange(nx_i):
  x_i[j] = float(infile.readline())
message= infile.readline()
print message

for i in arange(n+1):
  for j in arange(nx_i):
    lag_i[i][j] = float(infile.readline())
message= infile.readline()
print message

for j in arange(nx_i):
  psi_inter_i[j] = float(infile.readline())
message= infile.readline()
print message

for i in arange(n+1):
  for j in arange(nx_i):
    dlag_i[i][j] = float(infile.readline())
message= infile.readline()
print message

for j in arange(nx_i):
  diff_psi_i[j] = float(infile.readline())
message= infile.readline()
print message

#for i in arange(n+1):
#  for j in arange(n+1):
#    ilag_i[i][j] = float(infile.readline())
#message= infile.readline()
#print message
#
for j in arange(nx_i):
  int_psi_i[j] = float(infile.readline())
message= infile.readline()
print message

# full level
for j in arange(nx_m):
  x_m[j] = float(infile.readline())
message= infile.readline()
print message

for i in arange(n):
  for j in arange(nx_m):
    lag_m[i][j] = float(infile.readline())
message= infile.readline()
print message

for j in arange(nx_m):
  psi_inter_m[j] = float(infile.readline())
message= infile.readline()
print message

for i in arange(n):
  for j in arange(nx_m):
    dlag_m[i][j] = float(infile.readline())
message= infile.readline()
print message

for j in arange(nx_m):
  diff_psi_m[j] = float(infile.readline())
message= infile.readline()
print message

#for i in arange(n):
#  for j in arange(n):
#    ilag_m[i][j] = float(infile.readline())
#message= infile.readline()
#print message
#
for j in arange(nx_m):
  int_psi_m[j] = float(infile.readline())
message= infile.readline()
print message


# analytic solution
for j in arange(nx_i):
  ana_diff_i[j] = float(infile.readline())
message= infile.readline()
print message

for j in arange(nx_i):
  ana_int_i[j] = float(infile.readline())
message= infile.readline()
print message

for j in arange(nx_m):
  ana_diff_m[j] = float(infile.readline())
message= infile.readline()
print message

for j in arange(nx_m):
  ana_int_m[j] = float(infile.readline())
message= infile.readline()
print message

infile.close()


for i in arange(n+1):
  for j in arange(nx_i):
    wlag_i[i][j]  = psi_i[i]*lag_i[i][j]
    wdlag_i[i][j] = psi_i[i]*dlag_i[i][j]

for i in arange(n):
  for j in arange(nx_m):
    wlag_m[i][j]  = psi_m[i]*lag_m[i][j]
    wdlag_m[i][j] = psi_m[i]*dlag_m[i][j]


#### Plot #####

styles = array([' '*35]*11)
styles[0]  = 'r-'
styles[1]  = 'g-'
styles[2]  = 'b-'
styles[3]  = 'y-'
styles[4]  = 'k-'
styles[4]  = 'c-'
styles[5]  = 'r--'
styles[6]  = 'g--'
styles[7]  = 'b--'
styles[8]  = 'y--'
styles[9]  = 'k--'
styles[10] = 'k--'


## Half Level Plot

pltfig_i   = plt.figure(figsize=(14,18))
pltfig_i.suptitle('Lagrange Polynomial', fontsize=18)
pltgrid = gridspec.GridSpec(4, 2)


#- plot (0,0)
axis = plt.subplot(pltgrid[0,0])
#axis.set_title(r'$\psi_i\ -\ \eta_i$')
axis.set_xlabel(r'$\eta_j$')
axis.set_ylabel(r'$\psi_j$')
axis.set_xlim(0,1)
axis.plot(eta_i, psi_i, 'ro')
#axis.set_ylim(ymin[i][j],ymax[i][j])


#- plot (0,1)
axis = plt.subplot(pltgrid[0,1])
#axis.set_title(r'Lagrange Polynomial($L_j (\eta)$)')
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$L_j (\eta)$')
axis.set_xlim(0,1)
axis.axhline(y=0.0, color='k', ls='--')
axis.axhline(y=1.0, color='k', ls='--')
for i in range(n+1):
   axis.axvline(x=eta_i[i], color='k', ls='--')

for i in range(n+1):
   axis.plot(x_i, lag_i[i], styles[i])

#- plot (1,0)
axis = plt.subplot(pltgrid[1,0])
#axis.set_title(r'$\psi_j L_j (\eta)$')
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$\psi_j L_j (\eta)$')
axis.set_xlim(0,1)
axis.axhline(y=0.0, color='k', ls='--')
for i in range(n+1):
   axis.axhline(y=psi_i[i], color='k', ls='--')
   axis.axvline(x=eta_i[i], color='k', ls='--')

for i in range(n+1):
   axis.plot(x_i, wlag_i[i], styles[i])
 

#- plot (1,1)
axis = plt.subplot(pltgrid[1,1])
#axis.set_title(r'$\psi(\eta)$')
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$\psi (\eta)$')
axis.set_xlim(0,1)
axis.plot(eta_i, psi_i, 'ro')
axis.plot(x_i, psi_inter_i, 'g-')


#- plot (2,0)
axis = plt.subplot(pltgrid[2,0])
#axis.set_title(r'Lagrange Polynomial($L_j (\eta)$)')
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$\psi_j \frac{\partial{L_j (\eta)}} {\partial \eta}$')
axis.set_xlim(0,1)
axis.axhline(y=0.0, color='k', ls='--')
axis.axhline(y=1.0, color='k', ls='--')
for i in range(n+1):
   axis.axvline(x=eta_i[i], color='k', ls='--')

for i in range(n+1):
   axis.plot(x_i, wdlag_i[i], styles[i])


#- plot (2,1)
axis = plt.subplot(pltgrid[2,1])
#axis.set_title(r'$\frac{d \psi(\eta)}{d \eta}$')
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$\frac{\partial{\psi (\eta)}} {\partial \eta}$')
axis.set_xlim(0,1)
axis.plot(x_i, diff_psi_i, 'r-', label='Numeriacal')
axis.plot(x_i, ana_diff_i, 'g--', label='Analytic')
axis.legend(shadow=True, loc=2)




#- plot (3,0)
axis = plt.subplot(pltgrid[3,0])
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$\int L_j (\eta) \mathbf{d} x$')
axis.set_xlim(0,1)
#axis.axhline(y=0.0, color='k', ls='--')
#axis.axhline(y=1.0, color='k', ls='--')
for i in range(n+1):
   axis.axvline(x=eta_i[i], color='k', ls='--')

for i in range(n+1):
   axis.plot(eta_i, ilag_i[i], styles[i])


#- plot (3,1)
axis = plt.subplot(pltgrid[3,1])
#axis.set_title(r'$\psi(\eta)$')
#axis.set_ylim(0, 10)
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$\int \psi (\eta) \mathbf{d} x$')
axis.set_xlim(0,1)
axis.plot(x_i, int_psi_i, 'r-', label='Numeriacal')
axis.plot(x_i, ana_int_i, 'g--', label='Analytic')
axis.legend(shadow=True, loc=2)




## Full Level Plot

pltfig_m   = plt.figure(figsize=(14,18))
pltfig_m.suptitle('Lagrange Polynomial', fontsize=18)
pltgrid = gridspec.GridSpec(4, 2)


#- plot (0,0)
axis = plt.subplot(pltgrid[0,0])
#axis.set_title(r'$\psi_m\ -\ \eta_m$')
axis.set_xlabel(r'$\eta_j$')
axis.set_ylabel(r'$\psi_j$')
axis.set_xlim(0,1)
axis.plot(eta_m, psi_m, 'ro')
#axis.set_ylim(ymin[i][j],ymax[i][j])


#- plot (0,1)
axis = plt.subplot(pltgrid[0,1])
#axis.set_title(r'Lagrange Polynomial($L_j (\eta)$)')
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$L_j (\eta)$')
axis.set_xlim(0,1)
axis.axhline(y=0.0, color='k', ls='--')
axis.axhline(y=1.0, color='k', ls='--')
for i in range(n):
   axis.axvline(x=eta_m[i], color='k', ls='--')

for i in range(n):
   axis.plot(x_m, lag_m[i], styles[i])

#- plot (1,0)
axis = plt.subplot(pltgrid[1,0])
#axis.set_title(r'$\psi_j L_j (\eta)$')
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$\psi_j L_j (\eta)$')
axis.set_xlim(0,1)
axis.axhline(y=0.0, color='k', ls='--')
for i in range(n):
   axis.axhline(y=psi_m[i], color='k', ls='--')
   axis.axvline(x=eta_m[i], color='k', ls='--')

for i in range(n):
   axis.plot(x_m, wlag_m[i], styles[i])
 

#- plot (1,1)
axis = plt.subplot(pltgrid[1,1])
#axis.set_title(r'$\psi(\eta)$')
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$\psi (\eta)$')
axis.set_xlim(0,1)
axis.plot(eta_m, psi_m, 'ro')
axis.plot(x_m, psi_inter_m, 'g-')


#- plot (2,0)
axis = plt.subplot(pltgrid[2,0])
#axis.set_title(r'Lagrange Polynomial($L_j (\eta)$)')
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$\psi_j \frac{\partial{L_j (\eta)}} {\partial \eta}$')
axis.set_xlim(0,1)
axis.axhline(y=0.0, color='k', ls='--')
axis.axhline(y=1.0, color='k', ls='--')
for i in range(n):
   axis.axvline(x=eta_m[i], color='k', ls='--')

for i in range(n):
   axis.plot(x_m, wdlag_m[i], styles[i])


#- plot (2,1)
axis = plt.subplot(pltgrid[2,1])
#axis.set_title(r'$\frac{d \psi(\eta)}{d \eta}$')
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$\frac{\partial{\psi (\eta)}} {\partial \eta}$')
axis.set_xlim(0,1)
axis.plot(x_m, diff_psi_m, 'r-', label='Numerical')
axis.plot(x_m, ana_diff_m, 'g--', label='Analytic')
axis.legend(shadow=True, loc=2)




#- plot (3,0)
axis = plt.subplot(pltgrid[3,0])
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$\int L_j (\eta) \mathbf{d} x$')
axis.set_xlim(0,1)
#axis.axhline(y=0.0, color='k', ls='--')
#axis.axhline(y=1.0, color='k', ls='--')
for i in range(n):
   axis.axvline(x=eta_m[i], color='k', ls='--')

for i in range(n):
   axis.plot(eta_m, ilag_m[i], styles[i])


#- plot (3,1)
axis = plt.subplot(pltgrid[3,1])
#axis.set_title(r'$\psi(\eta)$')
#axis.set_ylim(0, 10)
axis.set_xlabel(r'$\eta$')
axis.set_ylabel(r'$\int \psi (\eta) \mathbf{d} x$')
axis.set_xlim(0,1)
axis.plot(x_m, int_psi_m, 'r-', label='Numerical')
axis.plot(x_m, ana_int_m, 'g--', label='Analytic')
axis.legend(shadow=True, loc=2)






# Save
#outfilename = ''+str(testcase)+'.png'
outfilename = 'PolynomialInterpolation_i'+'.png'
pltfig_i.savefig(outfilename)
outfilename = 'PolynomialInterpolation_m'+'.png'
pltfig_m.savefig(outfilename)

plt.show(True)
