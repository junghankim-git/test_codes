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


class SCRIP_Input:

  def __init__(self, filename):
    self.filename = filename
    self.LoadFile()

    self.xmin =   0.0
    self.xmax = 360.0
    self.ymin = -90.0
    self.ymax =  90.0

  def LoadFile(self):
    self.infile       = nc.Dataset(self.filename, mode='r')
    self.grid_size    = len(self.infile.dimensions['grid_size'])
    self.grid_corners = len(self.infile.dimensions['grid_corners'])

    # dimensions
    size    = self.grid_size
    corners = self.grid_corners
    
    # allocate
    self.center_lats = [0.0 for i in range(size)]
    self.center_lons = [0.0 for i in range(size)]
    self.corner_lats = [[0.0 for j in range(corners)] for i in range(size)]
    self.corner_lons = [[0.0 for j in range(corners)] for i in range(size)]
    
    # variables
    self.center_lats = self.infile.variables['grid_center_lat'][:]*180.0/np.pi
    self.center_lons = self.infile.variables['grid_center_lon'][:]*180.0/np.pi
    self.corner_lats = self.infile.variables['grid_corner_lat'][:]*180.0/np.pi
    self.corner_lons = self.infile.variables['grid_corner_lon'][:]*180.0/np.pi
    
    # close
    self.infile.close()


  def IniAxis(self, axis, xmin_in, xmax_in, ymin_in, ymax_in, ismap):
     axis.clear()
     axis.set_xlim(xmin_in,xmax_in)
     axis.set_ylim(ymin_in,ymax_in)
     if ismap:
       axis.set_frame_on(False)
       axis.set_yticklabels([])
       axis.set_xticklabels([])
       axis.set_yticks([])
       axis.set_xticks([])


  def IniAxis_No(self, axis):
     axis.clear()
     axis.set_frame_on(False)
     #axis.set_xlim(xmin_in,xmax_in)
     #axis.set_ylim(ymin_in,ymax_in)
     axis.set_yticklabels([])
     axis.set_xticklabels([])
     axis.set_yticks([])
     axis.set_xticks([])


  def CreateCanvas(self, title, ismap):
    pltfig = plt.figure(figsize=(16,12))
    pltfig.clear()
    axis  = pltfig.add_subplot(1,1,1)
    self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax, ismap)
    axis.set_title(title,fontsize=25)
    return pltfig, axis


  def Draw(self, title, ismap): #, alpha, beta, lon, lat):

    size    = self.grid_size
    corners = self.grid_corners

    pltfig, axis = self.CreateCanvas(title, ismap)

    self.draw_map(axis, ismap)
    '''
    if ismap:
      self.draw_map(axis)
    else:
      self.draw_cart(axis)
    '''

    # save fig
    pltfig.savefig('./grid.png')
    del axis


  def draw_map(self, axis, ismap):

    lon0=127.0
    lat0=37.30
#    lon0=0.0
#    lat0=-90.0
    if ismap:
      map = Basemap(projection='ortho',lon_0=lon0,lat_0=lat0,resolution='l')  # seoul
    else:
      map = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=0, urcrnrlon=360, resolution='c')
 
    # draw coastlines, country boundaries, fill continents.
    map.drawcountries(linewidth=0.8)
    map.drawcoastlines(linewidth=0.8)
    
    #map.fillcontinents(color='coral',lake_color='aqua')
    #map.fillcontinents(color='white',lake_color='white')
    
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='aqua')
    map.drawmapboundary(fill_color='white')
    
    # draw lat/lon grid lines every 30 degrees.
    #map.drawmeridians(np.arange(0,360,30))
    #map.drawparallels(np.arange(-90,90,30))

    ##### Draw centers #####
    x, y = map(self.center_lons,self.center_lats)
    axis.scatter(x,y,color='r',s=1)

    ##### Draw corners #####
    #for i in range(self.grid_size):
    for i in range(10):
      x, y = map(self.corner_lons[i][:],self.corner_lats[i][:])
      axis.scatter(x,y,color='b',s=1)

    x, y = map(self.corner_lons[-1][:],self.corner_lats[-1][:])
    axis.scatter(x,y,color='b',s=1)
    #x, y = map(self.corner_lons[-2][:],self.corner_lats[-2][:])
    #axis.scatter(x,y,color='b',s=1)
    #x, y = map(self.corner_lons[-400][:],self.corner_lats[-400][:])
    #axis.scatter(x,y,color='b',s=1)


  def draw_cart(self, axis):

    ##### Draw centers #####
    x, y = self.center_lons, self.center_lats
    axis.scatter(x,y,color='r',s=1)

    ##### Draw corners #####
    #for i in range(self.grid_size):
    #  x, y = self.corner_lons[i][:], self.corner_lats[i][:]
    #  axis.scatter(x,y,s=1)
    #  axis.scatter(x,y,s=1)


 

  def Show(self):
    plt.show()
