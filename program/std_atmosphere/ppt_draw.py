#!/usr/bin/env python
import os
import sys
import time
from numpy import *
from scipy import *
import matplotlib.pyplot as plt
sys.path.append('/home/jhkim/work/share/python')
from utilities import create_plot,axis_init



nlayers = 100
nlevs   = 137+1
p0      = 1000.0
ps      =  920.0
ps2     =  950.0

H     = [0.0 for i in range(nlayers)]
T     = [0.0 for i in range(nlayers)]

P     = [0.0 for i in range(nlevs)]
H1    = [0.0 for i in range(nlevs)]
H1_e  = [0.0 for i in range(nlevs)]

Pts   = [0.5 for i in range(nlevs)]


filename = 'H-T.ascii'
infile = open(filename,'r')

k = 0
for iline in infile:
   line = iline.strip()
   splitline = str(line).split()
   H[k] = float(splitline[0])/1000.0
   T[k] = float(splitline[1])
   k = k + 1
infile.close()

print k,nlayers
if k != nlayers:
   print 'check nlayer'
   quit()


filename = 'P-H.ascii'
infile = open(filename,'r')

k = 0
for iline in infile:
   line = iline.strip()
   splitline = str(line).split()
   P[k]  = float(splitline[0])/100.0
   H1[k] = float(splitline[1])/1000.0
   k = k + 1
infile.close()


print k,nlevs
if k != nlevs:
   print 'check nlayer'
   quit()


for k in range(nlevs):
   H1_e[k]  = H1[0]/float(nlevs-1)*float(nlevs-1-k)



######################################
## Plotting
######################################

# Canvas 1

ncols = 3
nrows = 1

nAtms   = 3
Atms    = ['Troposphere','Stratosphere','Mesosphere']
AtmCol  = ['b','g','r']
AtmLay  = [20.0,47.0,84.852]
AtmPre  = [54.7, 1.1, 0.37]


xmins   = [[0.0],[0.0],[0.0]]
xmaxs   = [[1.0],[1.0],[1.0]]
ymins   = [[0.0],[0.0],[0.0]]
ymaxs   = [[90.0],[90.0],[1000.0]]
titles  = [['Height'],['Height'],['Pressure']]
xlabels = [[''],[''],['']]
ylabels = [['Height [km]'],['Height [km]'],['Pressure [hPa]']]

pltfig,axis = create_plot(ncols,nrows,'vertical grid',yscale=2.0)

for i in range(ncols):
    for j in range(nrows):
        axis_init(axis[i][j],titles[i][j],xlabels[i][j],ylabels[i][j],xmins[i][j],xmaxs[i][j],ymins[i][j],ymaxs[i][j])
        axis[i][j].set_xticks([])
        if i == 0:
           axis[i][j].plot(Pts[::3],H1_e[::3],'r.')
        if i == 1:
           axis[i][j].plot(Pts[::3],H1[::3],'r.')
        if i == 2:
           axis[i][j].invert_yaxis()
           axis[i][j].plot(Pts[::3],P[::3],'r.')
        if i == 2:
           for iAtm in range(nAtms):
              axis[i][j].axhline(y=AtmPre[iAtm],xmin=0.0,xmax=1.0,c=AtmCol[iAtm],ls='--')
        else:
           for iAtm in range(nAtms):
              axis[i][j].axhline(y=AtmLay[iAtm],xmin=0.0,xmax=1.0,c=AtmCol[iAtm],ls='--')

outfilename = './vertical_grid.png'
pltfig.savefig(outfilename)



# Canvas 2

nHyb   = 30
iHyb   = 7
iSig   = 20
P_x    = [0.5 for i in range(nHyb)]
P_pre  = [0.0 for i in range(nHyb)]
P_sig  = [0.0 for i in range(nHyb)]
P_sig2 = [0.0 for i in range(nHyb)]
P_hyb  = [0.0 for i in range(nHyb)]

con_width = 0.06
con_left  = 0.5-con_width/2.0
con_right = 0.5+con_width/2.0
connect   = [[[con_left,con_right],[0.0,0.0]] for i in range(nHyb)]
connect2  = [[[con_left,con_right],[0.0,0.0]] for i in range(nHyb)]

for k in range(nHyb):
   P_pre[k] = p0*k/(nHyb-1)
   P_sig[k] = ps*k/(nHyb-1)
   P_sig2[k]= ps2*k/(nHyb-1)
   connect[k][1][0] = P_pre[k]
   connect[k][1][1] = P_sig[k]

dp   = P_sig[iSig] - P_pre[iHyb-1]
dlev = iSig - iHyb + 1

for k in range(nHyb):
   if k < iHyb:
      P_hyb[k] = P_pre[k]
   elif k >= iHyb and k < iSig:
      P_hyb[k] = P_pre[iHyb-1]+dp/float(dlev)*(k-iHyb+1)
   elif k >= iSig:
      P_hyb[k] = P_sig[k]
   print P_hyb[k]

for k in range(nHyb):
   connect2[k][1][0] = P_pre[k]
   connect2[k][1][1] = P_hyb[k]

nrows = 1
ncols = 4


xmins   = [[0.0],[0.0],[0.0],[0.0]]
xmaxs   = [[1.0],[1.0],[1.0],[1.0]]
ymins   = [[-3.0],[-3.0],[-3.0],[-3.0]]
ymaxs   = [[1003.0],[1003.0],[1003.0],[1003.0]]
titles  = [['Pure Pressure'],['Sigma'],['Hybrid(p0 and ps)'],['Hybrid']]
xlabels = [[''],[''],[''],['']]
ylabels = [['Pressure [hPa]'],['Pressure [hPa]'],['Pressure [hPa]'],['Pressure [hPa]']]



pltfig, axis = create_plot(ncols,nrows,'Vertical Grid',yscale=2.0)

for i in range(ncols):
   for j in range(nrows):
      axis_init(axis[i][j],titles[i][j],xlabels[i][j],ylabels[i][j],xmins[i][j],xmaxs[i][j],ymins[i][j],ymaxs[i][j])
      axis[i][j].set_xticks([])
      axis[i][j].invert_yaxis()
      if i == 0:
         for iline in range(nHyb):
            axis[i][j].axhline(y=P_pre[iline],xmin=0.0,xmax=1.0,c='g',ls='--')
      if i == 1:
         for iline in range(nHyb):
            axis[i][j].axhline(y=P_pre[iline],xmin=0.0,xmax=con_left,c='g',ls='--')
            axis[i][j].axhline(y=P_sig[iline],xmin=con_right,xmax=1.0,c='b',ls='--')
            axis[i][j].plot(connect[iline][0],connect[iline][1],'c--')
         rect = plt.Rectangle((con_right+0.02,P_sig[nHyb-1]+2.0),con_left-0.04,1000.0-P_sig[nHyb-1]-2.0,facecolor="#aaaaaa")
         axis[i][j].add_patch(rect)
      if i == 2:
         for iline in range(nHyb):
            if iline < iHyb:
               axis[i][j].axhline(y=P_pre[iline],xmin=0.0,xmax=1.0,c='g',ls='--')
               P_hyb[iline] = P_pre[iline]
            elif iline >= iHyb and iline < iSig:
               axis[i][j].axhline(y=P_pre[iline],xmin=0.0,xmax=con_left,c='g',ls='--')
               axis[i][j].axhline(y=P_sig[iline],xmin=con_right,xmax=1.0,c='b',ls='--')
               axis[i][j].plot(connect[iline][0],connect[iline][1],'c--')
            elif iline >= iSig:
               axis[i][j].axhline(y=P_sig[iline],xmin=0.0,xmax=1.0,c='b',ls='--')
               P_hyb[iline] = P_sig[iline]
         rect = plt.Rectangle((0.02,P_sig[nHyb-1]+2.0),1.0-0.04,1000.0-P_sig[nHyb-1]-2.0,facecolor="#aaaaaa")
         axis[i][j].add_patch(rect)
      if i == 3:
         for iline in range(nHyb):
            axis[i][j].axhline(y=P_pre[iline],xmin=0.0,xmax=con_left,c='g',ls='--')
            axis[i][j].axhline(y=P_hyb[iline],xmin=con_right,xmax=1.0,c='k',ls='--')
            axis[i][j].plot(connect2[iline][0],connect2[iline][1],'c--')
         rect = plt.Rectangle((con_right+0.02,P_sig[nHyb-1]+2.0),con_left-0.04,1000.0-P_sig[nHyb-1]-2.0,facecolor="#aaaaaa")
         axis[i][j].add_patch(rect)

outfilename = './vertical_grid_2.png'
pltfig.savefig(outfilename)




# Canvas 4

pmin    = 0.0
pmax    = 10.0
nlevs   = 11
P_x     = [0.5 for i in range(nlevs)]
P_HOMME = [0.0 for i in range(nlevs)]
P_ECMWF = [0.0 for i in range(nlevs)]

dpc = (pmax-pmin)/float(nlevs-1)
dp  = [0.0,1.5,1.3,1.2,1.0,0.9,0.8,0.6,0.8,0.9,1.0]
#dp  = [0.0,0.5,0.3,0.2,0.0,-0.1,-0.2,-0.4,-0.2,-0.1,0.0]

P_ECMWF[0] = pmin
for k in range(nlevs):
   P_HOMME[k] = pmin+float(pmax-pmin)/float(nlevs-1)*float(k)
   if k != 0: P_ECMWF[k] = P_ECMWF[k-1] + dp[k]

for k in range(nlevs):
   print P_HOMME[k],P_ECMWF[k]



nrows,ncols = 1,1

xmins   = [[0.0]]
xmaxs   = [[1.0]]
ymins   = [[-1.0]]
ymaxs   = [[11.0]]
titles  = [['KIAPSGM vs. ECMWF']]
xlabels = [['']]
ylabels = [['Pressure [hPa]']]

pltfig, axis = create_plot(1,1,'KIASPGM vs. ECMWF',yscale=2.0)

for i in range(nrows):
   for j in range(ncols):
      axis_init(axis,titles[i][j],xlabels[i][j],ylabels[i][j],xmins[i][j],xmaxs[i][j],ymins[i][j],ymaxs[i][j])
      axis.set_xticks([])
      axis.invert_yaxis()
      if j == 0:
         for iline in range(nlevs):
            axis.axhline(y=P_HOMME[iline],xmin=0.0,xmax=con_left,c='r',ls='--')
            axis.axhline(y=P_ECMWF[iline],xmin=con_right,xmax=1.0,c='b',ls='--')

outfilename = './vertical_grid_3.png'
pltfig.savefig(outfilename)
plt.show()


