#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
sys.path.append('/home/jhkim/work/share/python')
from utilities import *



#====================================================================
# KIM v2.5
#====================================================================
print ' '
print '# 1. KIM v2.5'
print ' - 1.1. variables'
datafile = './walltime.dat'
nps      = len_file(datafile)
#comps    = ['Total','Computation (dyn+phy)','Others (include I/O)']
comps    = ['Total','Computation (dyn+phy)','Write','Initialize']
ncomps   = len(comps)
nprocs   = []
walltime = [[] for i in range(ncomps)]
speedup  = [[] for i in range(ncomps)]

print ' - 1.2. read walltime and expectation'
'''
infile   = open(datafile, 'r')
for line in infile:
    sline = line.split()
    if sline[0] == '#': continue
    nprocs.append(int(sline[0]))
    walltime[2].append(float(sline[1]))                  # i/o
    walltime[0].append(float(sline[2]))                  # total
    walltime[1].append(float(sline[2])-float(sline[1]))  # computation
infile.close()
nps = len(nprocs)
for i in range(nps):
    print('    {:5d}: {:7.1f}, {:7.1f}, {:7.1f}'.format(nprocs[i],walltime[0][i],walltime[1][i],walltime[2][i]))
'''
infile   = open(datafile, 'r')
for line in infile:
    sline = line.split()
    if sline[0]=='#': continue
    nprocs.append(int(sline[0]))
    walltime[0].append(float(sline[1])) # total
    walltime[1].append(float(sline[2])) # comp
    walltime[2].append(float(sline[3])) # i/o
    print 30./float(sline[3])
    walltime[3].append(float(sline[4])) # init
    walltime[0][-1] = walltime[1][-1]+walltime[2][-1]+walltime[3][-1]  # remove bug
infile.close()
nps = len(nprocs)
for i in range(nps):
    print('    {:5d}: {:7.1f}, {:7.1f}, {:7.1f}'.format(nprocs[i],walltime[0][i],walltime[1][i],walltime[2][i]))


print ' - 1.3. speedup'
basewall = [0.0 for i in range(ncomps)]
for ic in range(ncomps):
    basewall[ic] = walltime[ic][0]*nprocs[0]
    for ip in range(nps):
        speedup[ic].append(basewall[ic]/walltime[ic][ip])

print ' '
print '# 2. KIM v3.0.09'
print ' - 2.1. variables'
#comps_   = ['Total','Computation (dyn+phy)','Output','Set+Ini','Others']
comps_   = ['Total','Computation (dyn+phy)','Write','Others (set+ini+...)']
ncomps   = len(comps_)
#nprocs_   = [240,10008,19968,40032,60000,69696]
nprocs_   = [10008,19968,40032,60000,69696]
nios      = [139,104,139,125,144]
nps       = len(nprocs_)
walltime_ = [[0.0 for j in range(nps)] for i in range(ncomps)]
speedup_  = [[0.0 for j in range(nps)] for i in range(ncomps)]
tags = ['total','set_atm_comp_driver','ini_atm_comp_driver','run_dynamics','run_physics','write_io_system(run)']

walltime_s = [0.0 for j in range(nps)]
speedup_s  = [0.0 for j in range(nps)]

for ip in range(nps):
    filename = 'prof/{:05d}/profile_summary.prf'.format(nprocs_[ip])
    values = get_kim_profile(filename,tags)
    walltime_[0][ip] = values[0]
    walltime_[1][ip] = values[3]+values[4]
    walltime_[2][ip] = values[5]
    walltime_[3][ip] = walltime_[0][ip]-walltime_[1][ip]-walltime_[2][ip]
#    walltime_[3][ip] = values[1]+values[2]
#    walltime_[4][ip] = walltime_[0][ip]-walltime_[1][ip]-walltime_[2][ip]-walltime_[3][ip]
    walltime_s[ip]   = values[0]-values[1]-values[2]
for i in range(nps):
    print('    {:5d}: {:7.1f}, {:7.1f}, {:7.1f}'.format(nprocs_[i],walltime_[0][i],walltime_[1][i],walltime_[2][i]))

basewall = [0.0 for i in range(ncomps)]
for ic in range(ncomps):
    basewall[ic] = walltime_[ic][0]*nprocs_[0]
    for ip in range(nps):
        speedup_[ic][ip] = basewall[ic]/walltime_[ic][ip]
basewall[0] = walltime_s[0]*nprocs_[0]
for ip in range(nps):
    speedup_s[ip] = basewall[0]/walltime_s[ip]


# fitting
popt1, pcov1 = curve_fit(speedup_func, nprocs_, speedup_[0])
popt2, pcov2 = curve_fit(speedup_func, nprocs_, speedup_s)
print('fitting {}, {}: {}, {}'.format(popt1,popt2,1.0/popt1,1.0/popt2))
nfs = 100
nprocs_f = [(i+1)*7000 for i in range(nfs)]
sps1_f   = [0.0 for i in range(nfs)]
sps2_f   = [0.0 for i in range(nfs)]
efs1_f   = [0.0 for i in range(nfs)]
efs2_f   = [0.0 for i in range(nfs)]
for i in range(nfs):
    np = nprocs_f[i]
    sps1_f[i] = speedup_func(np,popt1)[0]
    sps2_f[i] = speedup_func(np,popt2)[0]
    efs1_f[i] = sps1_f[i]/np*100.
    efs2_f[i] = sps2_f[i]/np*100.
    #print('{:06d}: {:7.1f}({:2.1f}%), {:7.1f}({:2.1f}%)'.format(np,sps1_f[i],efs1_f[i],sps2_f[i],efs2_f[i]))




print ' '
print '# 3. plots'
xmax1 = 35000.0
xmax2 = 75000.0
nprocs_id = [0.0,xmax2]
#====================================================================
# Plots
#====================================================================
print ' - plot 1'
xlabel  = 'number of cores [#]'
ylabel1 = 'wall-clock times [s]'
ylabel2 = 'speed-up [a.u.]'
carts = ''
pltfig,axis = create_plot(ncols=2,nrows=1,title='KIM scalability',yscale=0.9)
title1 = 'KIM v2.5 (NE240L100, w/o I/O decomposition)'
title2 = 'KIM v2.5 (NE240L100, w/o I/O decomposition)'
#draw_plots(axis[0][0],title1,carts,comps,nprocs,walltime,xlabel=xlabel,ylabel=ylabel1,xmax=xmax1)
#--
draw_plots(axis[0][0],title2,carts,comps[0],nprocs,speedup[0],xlabel=xlabel,ylabel=ylabel2,xmin=0.0,xmax=xmax1,ymin=0.0,ymax=xmax1)
add_plots(axis[0][0],'Ideal',nprocs_id,nprocs_id,['r-'],leg_loc=2)

print ' - plot 2'
title1 = 'KIM v3.0.09 (NE240L91, w/ I/O decomposition)'
title2 = 'KIM v3.0.09 (NE240L91, w/ I/O decomposition)'
#draw_plots(axis[0][1],title1,carts,comps_,nprocs_,walltime_,xlabel=xlabel,ylabel=ylabel1,xmin=0.0,xmax=xmax2,ymin=0.0)
#add_plots(axis[0][1],'Total (w/o set, ini)',nprocs_,walltime_s,['bo-'])
#for i in range(nps):
#    axis[0][1].text(nprocs_[i],80,str(nios[i]),color='k',size=10,va='center',ha='center')
#--
draw_plots(axis[1][0],title2,carts,comps_[0],nprocs_,speedup_[0],xlabel=xlabel,ylabel=ylabel2,xmin=0.0,xmax=xmax2,ymin=0.0,ymax=xmax2)
add_plots(axis[1][0],'Total (w/o set, ini)',nprocs_,speedup_s,['bo-'])
add_plots(axis[1][0],'Ideal',nprocs_id,nprocs_id,['r-'],leg_loc=2)

pltfig.savefig('fig_scal.png')

plt.show()



