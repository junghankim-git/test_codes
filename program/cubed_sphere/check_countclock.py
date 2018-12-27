#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


npts     = 5
chk_line = 100
#chk_line = 266581
#chk_line = 536531

infilename = './points.txt'
infile = open(infilename, 'r')


rlons = [0.0 for i in range(npts)]
rlats = [0.0 for i in range(npts)]
dlons = [0.0 for i in range(npts)]
dlats = [0.0 for i in range(npts)]

iline = 0
n2    = 0
for line in infile:
  iline = iline + 1
  sline = line.split()

  if iline == chk_line:
    for i in range(npts):
      rlons[i]  = float(sline[i])
      dlons[i]  = rlons[i]*180.0/np.pi
    for i in range(npts):
      rlats[i] = float(sline[npts+i])
      dlats[i]  = rlats[i]*180.0/np.pi

print dlons
print dlats

infile.close()


dlon0 = np.sum(dlons)/float(npts)
dlat0 = np.sum(dlats)/float(npts)
print dlon0, dlat0

ddlat = max(dlats)-min(dlats)
ddlon = max(dlons)-min(dlons)
delta = 0.7*max(ddlat, ddlon)


margin = 0.05
pltfig = plt.figure(figsize=(7,6))
pltfig.subplots_adjust(left=margin, bottom=margin, right=1.0-margin, top=1.0-margin)
axis   = pltfig.add_subplot(1,1,1)
axis.set_title('point = '+str(chk_line))


# test
#map = Basemap(projection='lcc',width=12000000,height=9000000,lat_1=dlat0-delta,lat_2=dlat0+delta,lat_0=dlat0,lon_0=dlon0,resolution='c')
#map = Basemap(projection='lcc',width=1200,height=900,lat_1=dlat0-delta,lat_2=dlat0+delta,lat_0=dlat0,lon_0=dlon0,resolution='c')
#map = Basemap(projection='cass', llcrnrlat=dlat0-delta, urcrnrlat=dlat0+delta, llcrnrlon=dlon0-delta, urcrnrlon=dlon0+delta,lat_0=dlat0,lon_0=dlon0, resolution='c')
# test

map = Basemap(projection='cyl', llcrnrlat=dlat0-delta, urcrnrlat=dlat0+delta, llcrnrlon=dlon0-delta, urcrnrlon=dlon0+delta, resolution='c')
map.drawcountries(linewidth=0.8)
map.drawcoastlines(linewidth=0.8)
map.drawmapboundary(fill_color='white')
map.drawmeridians(np.arange(0,360,delta/4.0),  labels=[False,False,False,True], fmt='%5.1f')
map.drawparallels(np.arange(-90,90,delta/4.0), labels=[True,False,False,False], fmt='%5.1f')
#map.drawmeridians(np.arange(0,360,1),  labels=[False,False,False,True], fmt='%5.1f')
#map.drawparallels(np.arange(-90,90,1), labels=[True,False,False,False], fmt='%5.1f')

x, y  = map(dlons, dlats)
map.scatter(x, y, color='r', s=15)
for i in range(npts):
  axis.text(dlons[i], dlats[i]-0.1*ddlat, str(i+1), color='k', size=12, va='center', style='italic')


plt.show()
