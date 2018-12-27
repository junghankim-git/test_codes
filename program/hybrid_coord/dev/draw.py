#!/usr/bin/env python 
from numpy import *
from scipy import *
#from scipy.io import FortranFile
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import cm
import os
import sys
import struct



# ecmwf
#nlevs_e = [62, 91, 137]
nlevs_e = [91, 137]
#nlevs_e = [137]
# grims
nlevs_g = [64]


nlevs_t = [100]
files_t = ['output/half_0100nlevs_0041s_0099e_00222K_00030pa_00057km.dat']
#files_t = ['output/half_0100nlevs_0041s_0099e_00222K_00001pa_00077km.dat']
#files_t = ['output/half_0100nlevs_0041s_0099e_00222K_00029pa_00057km.dat']
#files_t = ['test.dat']
#nlevs_t = [100,100]
#files_t = ['./output/half_0100nlevs_0041s_0099e_00222K_00030pa_00057km.dat','/home/jhkim/TestBed/KIM/3.0.micros/3.0.04/03.dev/Exp_10h/SW_G/vcoord/Half_0100nlevs_0030plev_0101slev_00222K_00030pa_00057km.dat']

basedir = './coefficient'


# ecmwf
necmwf = len(nlevs_e)
print('necmwf = %d'%necmwf)
lev_e  = [[0.0 for j in range(nlevs_e[i]+1)] for i in range(necmwf)]
coa_e  = [[0.0 for j in range(nlevs_e[i]+1)] for i in range(necmwf)]
cob_e  = [[0.0 for j in range(nlevs_e[i]+1)] for i in range(necmwf)]
eta_e  = [[0.0 for j in range(nlevs_e[i]+1)] for i in range(necmwf)]

ecmwf_file = ['' for i in range(necmwf)]
for i in range(necmwf):
    levstr        = '{0:04d}'.format(nlevs_e[i])+'nlevs'
    ecmwf_file[i] = basedir+'/ecmwf/ecmwf_half_'+levstr+'.ascii'
    file = open(ecmwf_file[i], 'r')
    ilev  = 0
    for iline in file:
        line = iline.strip()
        splitline = str(line).split()
        coa_e[i][ilev] = float(splitline[0])/100000.0
        cob_e[i][ilev] = float(splitline[1])
        eta_e[i][ilev] = coa_e[i][ilev]+cob_e[i][ilev]
        lev_e[i][ilev] = float(ilev)/float(nlevs_e[i])
        ilev = ilev + 1
    file.close()
# grims
ngrims = len(nlevs_g)
print('ngrims = %d'%ngrims)
lev_g  = [[0.0 for j in range(nlevs_g[i]+1)] for i in range(ngrims)]
coa_g  = [[0.0 for j in range(nlevs_g[i]+1)] for i in range(ngrims)]
cob_g  = [[0.0 for j in range(nlevs_g[i]+1)] for i in range(ngrims)]
eta_g  = [[0.0 for j in range(nlevs_g[i]+1)] for i in range(ngrims)]

grims_file = ['' for i in range(ngrims)]
for i in range(ngrims):
    levstr        = '{0:04d}'.format(nlevs_g[i])+'nlevs'
    grims_file[i] = basedir+'/grims/grims_half_'+levstr+'.ascii'
    file = open(grims_file[i], 'r')
    ilev  = 0
    for iline in file:
        line = iline.strip()
        splitline = str(line).split()
        coa_g[i][ilev] = float(splitline[0])/100000.0
        cob_g[i][ilev] = float(splitline[1])
        eta_g[i][ilev] = coa_g[i][ilev]+cob_g[i][ilev]
        lev_g[i][ilev] = float(ilev)/float(nlevs_g[i])
        ilev = ilev + 1
    file.close()


# my test
ntest = len(nlevs_t)
print('ntest = %d'%ntest)
if ntest!=len(files_t):
    print 'some error'
    quit()
lev_t  = [[0.0 for j in range(nlevs_t[i]+1)] for i in range(ntest)]
coa_t  = [[0.0 for j in range(nlevs_t[i]+1)] for i in range(ntest)]
cob_t  = [[0.0 for j in range(nlevs_t[i]+1)] for i in range(ntest)]
eta_t  = [[0.0 for j in range(nlevs_t[i]+1)] for i in range(ntest)]
lev_m  = [[0.0 for j in range(nlevs_t[i]  )] for i in range(ntest)]
deta_m = [[0.0 for j in range(nlevs_t[i]  )] for i in range(ntest)]


for i in range(ntest):
    file = open(files_t[i],'rb')
    print('reading file: %s'%(files_t[i]))
    for k in range(nlevs_t[i]+1):
        byte1 = file.read(4)
        byte2 = file.read(8)
        byte3 = file.read(4)
        byte4 = file.read(4)
        byte5 = file.read(8)
        byte6 = file.read(4)
        coa_t[i][k] = struct.unpack('d',byte2)[0]
        cob_t[i][k] = struct.unpack('d',byte5)[0]
        eta_t[i][k] = coa_t[i][k]+cob_t[i][k]
        lev_t[i][k] = float(k)/float(nlevs_t[i])
    file.close()
    for k in range(nlevs_t[i]):
        lev_m[i][k]  = 0.5*(lev_t[i][k]+lev_t[i][k+1])
        deta_m[i][k] = eta_t[i][k+1]-eta_t[i][k]


if ntest==2 and nlevs_t[0]==nlevs_t[1]:
    for k in range(nlevs_t[0]+1):
        print('%4d,%10.7f,%10.7f,%10.7f,%10.7f,%10.7f,%10.7f'%(k,coa_t[0][k],coa_t[1][k],cob_t[0][k],cob_t[1][k],eta_t[0][k],eta_t[1][k]))


#####################
# plot
#####################
margin = 0.08
wspace = margin*1.3
hspace = 0.25
pltfig = plt.figure(figsize=(12,8))
pltfig.subplots_adjust(left=margin,right=1.0-margin,bottom=1.2*margin,top=1.0-1.2*margin,wspace=wspace,hspace=hspace)

fmts = ['-',':','--','-.']
xmin =  0.0
ymin = -0.1

# ecmwf
axis1 = pltfig.add_subplot(2,2,1)
axis1.set_title('test')
axis1.set_xlabel('nomalized level [a.u.]')
axis1.set_ylabel(r'$\eta$ [a.u.]')
axis1.set_xlim(xmin,1)
axis1.set_ylim(ymin,1)

for i in range(necmwf):
    label = 'ecmwf-%04d ('%(nlevs_e[i])
    axis1.plot(lev_e[i],coa_e[i],'r'+fmts[i],label=label+r'$A$'+')')
    axis1.plot(lev_e[i],cob_e[i],'b'+fmts[i],label=label+r'$B$'+')')
    axis1.plot(lev_e[i],eta_e[i],'k'+fmts[i],label=label+r'$\eta$'+')')
for i in range(ngrims):
    label = 'grims-%04d ('%(nlevs_g[i])
    axis1.plot(lev_g[i],coa_g[i],'r'+'-.',label=label+r'$A$'+')')
    axis1.plot(lev_g[i],cob_g[i],'b'+'-.',label=label+r'$B$'+')')
    axis1.plot(lev_g[i],eta_g[i],'k'+'-.',label=label+r'$\eta$'+')')
axis1.legend(loc=2,shadow=True,prop={'size':9})

# test
axis2 = pltfig.add_subplot(2,2,2)
axis2.set_title('test')
axis2.set_xlabel('nomalized level [a.u.]')
axis2.set_ylabel(r'$\eta$ [a.u.]')
axis2.set_xlim(xmin,1)
axis2.set_ylim(ymin,1)

for i in range(ntest):
    label = 'test(%d)-%04d ('%(i+1,nlevs_t[i])
    axis2.plot(lev_t[i],coa_t[i],'r'+fmts[i],label=label+r'$A$'+')')
    axis2.plot(lev_t[i],cob_t[i],'b'+fmts[i],label=label+r'$B$'+')')
    axis2.plot(lev_t[i],eta_t[i],'k*'+fmts[i],label=label+r'$\eta$'+')')
axis2.legend(loc=2,shadow=True,prop={'size':9})

# deta
axis3 = pltfig.add_subplot(2,2,3)
axis3.set_title('deta')
axis3.set_xlabel('nomalized level [a.u.]')
axis3.set_ylabel(r'$d\eta$ [a.u.]')
axis3.set_xlim(xmin,1)
axis3.set_ylim(ymin,1.2*max(max(deta_m)))

for i in range(ntest):
    label = 'test(%d)-%04d ('%(i+1,nlevs_t[i])
    axis3.plot(lev_m[i],deta_m[i],'k'+fmts[i],label=label+r'd$\eta$'+')')
axis3.legend(loc=2,shadow=True,prop={'size':9})


pltfig.savefig('test.png')
plt.show()



def polynomial(x, arg):
    w=arg[0]+arg[1]*x+arg[2]*power(x,2)+arg[3]*power(x,3)+arg[4]*power(x,4)
    return w

def l_polynomial(x, arg):
    w=arg[0]*x+arg[1]*power(x,2)+arg[2]*power(x,3)+arg[3]*power(x,4)
    return w

def high_polynomial(x, arg):
    order = len(arg)-1
    w = 0.0
    for i in range(order+1):
        w = w + arg[order-i]*power(x, i)
    return w

def gaussian(x, a, mu, sigma):
    return a*exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma))

def gaussian2(x, a, sigma):
    return a*exp(-(x-1.0)*(x-1.0)/(2.0*sigma*sigma))
