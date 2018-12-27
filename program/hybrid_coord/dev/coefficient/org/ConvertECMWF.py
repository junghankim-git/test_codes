#!/usr/bin/env python
from numpy import *
from scipy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time
from matplotlib import cm
import os
import sys
# load My Class
sys.path.append('/home/jhkim/Study/Library/Shared/Python')


#####################
###  Levels
#####################
# 19, 31, 40, 50, 60, 62, 91, 137
nlevs  = 137
lecmwf = True
coefA_h = [0.0 for i in range(nlevs+1)]
coefB_h = [0.0 for i in range(nlevs+1)]
coefA_f = [0.0 for i in range(nlevs)]
coefB_f = [0.0 for i in range(nlevs)]

p0 = 100000.0

if lecmwf:
   basedir = './ECMWF'
else:
   basedir = './KIAPSGM'

levstr     = '{0:04d}'.format(nlevs)+'nlevs.ascii'
inputfile  = basedir+'/web/ecmwf_'+levstr

file       = open(inputfile, 'r')

ilev = 0
for iline in file:
   line          = iline.strip()
   splitline     = str(line).split()
   coefA_h[ilev] = float(splitline[1])
   coefB_h[ilev] = float(splitline[2])
   ilev = ilev + 1
file.close()
if ilev != nlevs+1:
   print 'Check: nlevs'
   quit()

for ilev in range(nlevs):
   coefA_f[ilev] = 0.5*(coefA_h[ilev+1]+coefA_h[ilev])
   coefB_f[ilev] = 0.5*(coefB_h[ilev+1]+coefB_h[ilev])

hvifile = basedir+'/ecmwf_half_'+levstr
hvmfile = basedir+'/ecmwf_full_'+levstr

outifile = open(hvifile, 'w')
for ilev in range(nlevs+1):
   outifile.write(str(coefA_h[ilev])+' '+str(coefB_h[ilev])+'\n')
outifile.close()

outmfile = open(hvmfile, 'w')
for ilev in range(nlevs):
   outmfile.write(str(coefA_f[ilev])+' '+str(coefB_f[ilev])+'\n')
outmfile.close()





# Statistic

Rgas    = 287.04
#Tmean   = 222.7
Tmean   = 222.56
gravity = 9.80616
Hs      = Rgas*Tmean/gravity

TherMin =      0.0
TherMax =      1.0
MesoMin =      1.0
MesoMax =    100.0
StraMin =    100.0
StraMax =  10000.0
TropMin =  10000.0
TropMax = 100000.0
nTher = 0
nMeso = 0
nStra = 0
nTrop = 0

p0 = 100000.0
ps = 100000.0
p_h = [0.0 for i in range(nlevs+1)]
p_f = [0.0 for i in range(nlevs)]
h_h = [0.0 for i in range(nlevs+1)]
h_f = [0.0 for i in range(nlevs)]

for ilev in range(nlevs+1):
   p_h[ilev] = coefA_h[ilev] + coefB_h[ilev]*ps
   if ilev != 0: h_h[ilev] = Hs*log(p0/p_h[ilev])
   print p_h[ilev]/100.0
   #print h_h[ilev]/1000.0
print ' '
for ilev in range(nlevs):
   p_f[ilev] = 0.5*(p_h[ilev+1]+p_h[ilev])
   h_f[ilev] = Hs*log(p0/p_f[ilev])
   if TherMin < p_f[ilev] and p_f[ilev] <= TherMax: nTher = nTher + 1
   if MesoMin < p_f[ilev] and p_f[ilev] <= MesoMax: nMeso = nMeso + 1
   if StraMin < p_f[ilev] and p_f[ilev] <= StraMax: nStra = nStra + 1
   if TropMin < p_f[ilev] and p_f[ilev] <= TropMax: nTrop = nTrop + 1
   print p_f[ilev]/100.0
   #print h_f[ilev]/1000.0
print ' '


print 'nlevs    = ', nlevs, nTher+nMeso+nStra+nTrop
print 'n layers = ', nTher, nMeso, nStra, nTrop
print 'fraction = ', float(nTher)/float(nlevs)*100.0, float(nMeso)/float(nlevs)*100.0, float(nStra)/float(nlevs)*100.0, float(nTrop)/float(nlevs)*100.0


