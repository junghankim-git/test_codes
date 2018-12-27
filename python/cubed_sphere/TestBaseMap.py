#!/usr/bin/env python
from scipy import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_toolkits.basemap import Basemap
import time
from matplotlib import cm
import os
import sys
#font
import matplotlib.font_manager as fontm


# load My Class
sys.path.append('/home/jhkim/Study/Library/Shared/Python')
from LibSpaceFillingCurve import *
from LibControlMesh import *
from LibCubedSphere import *

pltfig = plt.figure(figsize=(12,12))
axis   = pltfig.add_subplot(1,1,1)
axis.set_title('')


# seoul
lon0=127.0
lat0=37.30
# Face: 1 
#lon0=  0.0
#lat0=  0.0
# Face: 5 
#lon0=  0.0
#lat0=-90.0
# Face: 6 
#lon0=  0.0
#lat0= 90.0

map = Basemap(projection='ortho',lon_0=lon0,lat_0=lat0,resolution='l')  # seoul


# draw coastlines, country boundaries, fill continents.
#map.drawcoastlines(linewidth=0.25)
#map.drawcountries(linewidth=0.25)
map.drawcoastlines(linewidth=0.8)
#map.drawcountries(linewidth=1.0)

#map.fillcontinents(color='coral',lake_color='aqua')
#map.fillcontinents(color='white',lake_color='white')

# draw the edge of the map projection region (the projection limb)
#map.drawmapboundary(fill_color='aqua')
map.drawmapboundary(fill_color='white')


# draw lat/lon grid lines every 30 degrees.
#map.drawmeridians(np.arange(0,360,30))
#map.drawparallels(np.arange(-90,90,30))




ne     = 4

cube   = CubedSphere(4,ne,1)
nrows  = cube.nrows
ncols  = cube.ncols
cube.Ini()

iface  = cube.CubeToPlot(ival=2)
cidx   = cube.CubeToPlot(ival=6)
ab     = cube.CubeToPlot(ival=7)
lonlat = cube.CubeToPlot(ival=8)

lons   = [0.0 for i in range(4)]
lats   = [0.0 for i in range(4)]


def GetAlphaBetaLine(iface_in, npts_in, alpha_beta1, alpha_beta2):
   res_line = [[0.0 for i in range(npts_in)] for j in range(2)]
   a1 = alpha_beta1[0]
   b1 = alpha_beta1[1]
   a2 = alpha_beta2[0]
   b2 = alpha_beta2[1]

   if (a1 == a2) and (b1 != b2):
      db = abs(b2-b1)/(npts_in-1)
      for ip in range(npts_in):
         beta = min(b1,b2)+ip*db
         res_line[0][ip], res_line[1][ip] = cube.AlphaBeta2LonLat(iface_in, a1, beta)
   elif (a1 != a2) and (b1 == b2):
      da = abs(a2-a1)/(npts_in-1)
      for ip in range(npts_in):
         alpha = min(a1,a2)+ip*da
         res_line[0][ip], res_line[1][ip] = cube.AlphaBeta2LonLat(iface_in, alpha, b1)
   else:
      print 'Error : GetAlphaBetaLine...'
      quit()
   return res_line

def GetAlphaBetaSquare(iface_in, npts_in, alpha_beta):
   res_lines = [[[0.0 for i in range(npts_in)] for j in range(2)] for k in range(4)]
   for k in range(4):
      if k == 3:
         res_lines[k] = GetAlphaBetaLine(iface_in, npts_in, alpha_beta[k], alpha_beta[0])
      else:
         res_lines[k] = GetAlphaBetaLine(iface_in, npts_in, alpha_beta[k], alpha_beta[k+1])
   return res_lines


for i in range(nrows):
   for j in range(ncols):
      if iface[i][j] != -1:
#         if cidx[i][j][0]==3 and cidx[i][j][1]==4 and iface[i][j] == 6:
            for k in range(4):
               lons[k] = lonlat[i][j][k][0]/np.pi*180.0
               lats[k] = lonlat[i][j][k][1]/np.pi*180.0
            x, y = map(lons,lats)
            axis.scatter(x,y)

            npts  = 400/ne
            lines = GetAlphaBetaSquare(iface[i][j], npts, ab[i][j])
            for k in range(4):
               for ip in range(npts):
                  lines[k][0][ip] = lines[k][0][ip]/np.pi*180.0
                  lines[k][1][ip] = lines[k][1][ip]/np.pi*180.0
               x, y = map(lines[k][0], lines[k][1])
               axis.scatter(x,y,s=1)


pltfig.savefig('CS.png')

plt.show()



