#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import scipy as sp
import netCDF4 as nc
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_toolkits.basemap import Basemap
#from matplotlib import cm




class RemapMatrix_SCRIP:

  def __init__(self, filename):
    self.filename = filename

    self.xmin =   0.0
    self.xmax = 360.0
    self.ymin = -90.0
    self.ymax =  90.0

  def IniFile(self):
    # file open
    self.infile = nc.Dataset(self.filename, mode='r')

    # dimensions
    self.src_size    = len(self.infile.dimensions['src_grid_size'])
    self.dst_size    = len(self.infile.dimensions['dst_grid_size'])
    self.src_corners = len(self.infile.dimensions['src_grid_corners'])
    self.dst_corners = len(self.infile.dimensions['dst_grid_corners'])

    if self.infile.variables['src_grid_center_lon'].units == 'radians':
      self.src_lon = self.infile.variables['src_grid_center_lon'][:]*180./np.pi
      self.src_lat = self.infile.variables['src_grid_center_lat'][:]*180./np.pi
      self.dst_lon = self.infile.variables['dst_grid_center_lon'][:]*180./np.pi
      self.dst_lat = self.infile.variables['dst_grid_center_lat'][:]*180./np.pi
    else:
      self.src_lon = self.infile.variables['src_grid_center_lon'][:]
      self.src_lat = self.infile.variables['src_grid_center_lat'][:]
      self.dst_lon = self.infile.variables['dst_grid_center_lon'][:]
      self.dst_lat = self.infile.variables['dst_grid_center_lat'][:]

    self.infile.close()

  def GetSrcLonLat(self):
    return self.src_lon, self.src_lat

  def GetDstLonLat(self):
    return self.dst_lon, self.dst_lat
    
  def CloseFile(self):
    # close
    self.infile.close()





#  def IniAxis(self, axis, xmin_in, xmax_in, ymin_in, ymax_in): #, type):
#     #axis.clear()
#     axis.set_xlim(xmin_in,xmax_in)
#     axis.set_ylim(ymin_in,ymax_in)
#     #axis.set_frame_on(False)
#     #axis.set_yticklabels([])
#     #axis.set_xticklabels([])
#     #axis.set_yticks([])
#     #axis.set_xticks([])
#
#
#  def IniAxis_No(self, axis):
#     axis.clear()
#     axis.set_frame_on(False)
#     #axis.set_xlim(xmin_in,xmax_in)
#     #axis.set_ylim(ymin_in,ymax_in)
#     axis.set_yticklabels([])
#     axis.set_xticklabels([])
#     axis.set_yticks([])
#     axis.set_xticks([])
#
#
#  def CreateCanvas(self, title, type):
#    plt.title = title
#    pltfig = plt.figure(figsize=(10,20))
#    pltfig.clear()
#
#    axis1 = pltfig.add_subplot(2,1,1)
#    self.IniAxis(axis1, self.xmin, self.xmax, self.ymin, self.ymax)
#
#    axis2 = pltfig.add_subplot(2,1,2)
#    self.IniAxis(axis2, self.xmin, self.xmax, self.ymin, self.ymax)
#
#    return pltfig, axis1, axis2
#
#
#  def Draw(self, title, type=0): #, alpha, beta, lon, lat):
#
#    if type==0:
#      margin = 0.04     # 0.05 ~ 0.95 (for labels)
#      wspace = 2*margin
#      hspace = 0.01
#    elif type==1 or type==2:
#      margin = 0.06     # 0.05 ~ 0.95 (for labels)
#      wspace = 0.0
#      hspace = margin*1.2
#    else:
#      print 'unknown type..'
#      quit()
#    pltfig = plt.figure(figsize=(18,10))
#    pltfig.clear()
#    pltfig.suptitle(title,fontsize=18,y=1-margin+0.02) # default y = 0.98
#    pltfig.subplots_adjust(left=margin, bottom=margin, right=1.0-margin, top=1.0-margin, wspace=wspace, hspace=hspace)
#
#
#    print ' plot 1'
#    axis1 = pltfig.add_subplot(2,2,1)
#    self.IniAxis(axis1, self.xmin, self.xmax, self.ymin, self.ymax)
#    axis1.set_title('grid (src)')
#    self.draw_map(axis1, self.src_lons, self.src_lats, self.src_var, type, True)
#
#    print ' plot 2'
#    axis2 = pltfig.add_subplot(2,2,2)
#    self.IniAxis(axis2, self.xmin, self.xmax, self.ymin, self.ymax)
#    axis2.set_title('grid (dst)')
#    self.draw_map(axis2, self.dst_lons, self.dst_lats, self.dst_var, type, True)
#
#    print ' plot 3'
#    axis3 = pltfig.add_subplot(2,2,3)
#    self.IniAxis(axis3, self.xmin, self.xmax, self.ymin, self.ymax)
#    axis3.set_title('u10m (src)')
#    self.draw_map(axis3, self.src_lons, self.src_lats, self.src_var, type, False)
#
#    print ' plot 4'
#    axis4 = pltfig.add_subplot(2,2,4)
#    self.IniAxis(axis4, self.xmin, self.xmax, self.ymin, self.ymax)
#    axis4.set_title('u10m (dst)')
#    self.draw_map(axis4, self.dst_lons, self.dst_lats, self.dst_var, type, False)
#
#    # save fig
#    print ' saving...'
#    pltfig.savefig('./figure.png')
#    del axis1, axis2, axis3, axis4
#
#
#  def draw_map(self, axis, lons, lats, var, type, isgrid):
#
#    
#    lon0 = 127.0
#    lat0 = 37.30
#    res  = 'c'    # resolution: c<l<i<h<f
#
#    if type==0:
#      map = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=0, urcrnrlon=360, resolution=res)
#    elif type==1:
#      map = Basemap(projection='ortho',lon_0=lon0,lat_0=lat0,resolution=res)                           # seoul
#    elif type==2:
#      map = Basemap(width=12000000,height=9000000,projection='lcc',lat_1=lat0-4.0,lat_2=lat0+4.0,lat_0=lat0,lon_0=lon0,resolution=res)
# 
#    # draw coastlines, country boundaries, fill continents.
#    map.drawcountries(linewidth=0.8)
#    map.drawcoastlines(linewidth=0.8)
#    
#    # draw color for land  <= not use, because it remove field's color
#    #map.fillcontinents(color='coral',lake_color='aqua')
#    #map.fillcontinents(color='white',lake_color='white')
#    
#    # draw color for ocean
#    #map.drawmapboundary(fill_color='aqua')
#    map.drawmapboundary(fill_color='white')
#    
#    # draw lat/lon grid lines every 30 degrees. (labels=[left,right,top,down], default labelstyle = '+/-')
#    if type==1:
#      map.drawmeridians(np.arange(0,360,45),  labels=[False,False,False,False])
#      map.drawparallels(np.arange(-90,90,45), labels=[False,False,False,False])
#    else:
#      map.drawmeridians(np.arange(0,360,45),  labels=[False,False,False,False])
#      map.drawparallels(np.arange(-90,90,45), labels=[True,False,False,False])
#
#    
#    ##### Draw map #####
#    # latlon = True : use direct lon and lat
#    # tri    = True : use unstructed grid data
#    x, y  = map(lons, lats)
#    if isgrid:
#      axis.scatter(x, y, color='r', s=0.1)
#    else:
#      #clevs = arange(-30,30,10)
#      #cs    = map.colorbar(x, y, var, tri=True)
#      #cs    = map.contour(x, y, var, clevs, linewidths=0.5, colors='k', animated=True, latlon=True, tri=True)
#      #cs    = map.contour(lons, lats, var, latlon=True, tri=True)
#      #cs    = map.contour(x, y, var, tri=True)
#  
#      # Colorbar
#      im = map.pcolor(x, y, var, vmin=-20, vmax=20, shading='flat', cmap=plt.cm.jet, tri=True)
#      #cb = map.colorbar(im, location='bottom')  # <= colorbar
#      #cb.set_label('m/s')
#
#
#
#  def draw_cart(self, axis):
#
#    ##### Draw centers #####
#    x, y = self.center_lons, self.center_lats
#    axis.scatter(x,y,color='r',s=1)
#
#    ##### Draw corners #####
#    #for i in range(self.grid_size):
#    #  x, y = self.corner_lons[i][:], self.corner_lats[i][:]
#    #  axis.scatter(x,y,s=1)
#    #  axis.scatter(x,y,s=1)
#
#
# 
#
#  def Show(self):
#    plt.show()
