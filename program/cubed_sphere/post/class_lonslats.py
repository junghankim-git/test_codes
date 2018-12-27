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


class lonslats_plot:

  def __init__(self, filename, has_nbrs=False):
    self.filename = filename
    self.has_nbrs = has_nbrs
    self.load_file()

    self.xmin =   0.0
    self.xmax = 360.0
    self.ymin = -90.0
    self.ymax =  90.0

  def load_file(self):
    self.infile = nc.Dataset(self.filename, mode='r')

    # dimensions
    self.nups   = len(self.infile.dimensions['nups'])
    #self.nups   = len(self.infile.dimensions['ncol'])
    nups    = self.nups
    if self.has_nbrs:
      self.nnbrs  = len(self.infile.dimensions['nnbrs'])
      self.ncell  = len(self.infile.dimensions['ncell'])
      nnbrs   = self.nnbrs
      ncell   = self.ncell

    # allocate
    self.lons = [0.0 for i in range(nups)]
    self.lats = [0.0 for i in range(nups)]
    if self.has_nbrs:
      self.nbr_lons = [[0.0 for j in range(nnbrs)] for i in range(nups)]
      self.nbr_lats = [[0.0 for j in range(nnbrs)] for i in range(nups)]
      self.cell_lons = [[0.0 for j in range(ncell)] for i in range(nups)]
      self.cell_lats = [[0.0 for j in range(ncell)] for i in range(nups)]
    
    # variables
    self.lons = self.infile.variables['lons'][:]*180.0/np.pi
    self.lats = self.infile.variables['lats'][:]*180.0/np.pi
    if self.has_nbrs:
      self.nbr_lons = self.infile.variables['nbr_lons'][:]*180.0/np.pi
      self.nbr_lats = self.infile.variables['nbr_lats'][:]*180.0/np.pi
      self.cell_lons = self.infile.variables['cell_lons'][:]*180.0/np.pi
      self.cell_lats = self.infile.variables['cell_lats'][:]*180.0/np.pi
    
    # close
    self.infile.close()



  def axis_init(self, axis, xmin_in, xmax_in, ymin_in, ymax_in, ismap):
     axis.clear()
     axis.set_xlim(xmin_in,xmax_in)
     axis.set_ylim(ymin_in,ymax_in)
     if ismap:
       axis.set_frame_on(False)
       axis.set_yticklabels([])
       axis.set_xticklabels([])
       axis.set_yticks([])
       axis.set_xticks([])


  def axis_init_no(self, axis):
     axis.clear()
     axis.set_frame_on(False)
     #axis.set_xlim(xmin_in,xmax_in)
     #axis.set_ylim(ymin_in,ymax_in)
     axis.set_yticklabels([])
     axis.set_xticklabels([])
     axis.set_yticks([])
     axis.set_xticks([])


  def create_canvas(self, title, ismap):
    margin = 0.02
    pltfig = plt.figure(figsize=(16,8))
    pltfig.subplots_adjust(left=margin, bottom=margin, right=1.0-margin,top=1.0-margin-0.04, wspace=0.1, hspace=0.01)
    pltfig.clear()
    axis  = pltfig.add_subplot(1,1,1)
    self.axis_init(axis, self.xmin, self.xmax, self.ymin, self.ymax, ismap)
    axis.set_title(title,fontsize=25)
    return pltfig, axis


  def draw(self, title, ismap): #, alpha, beta, lon, lat):

    nups    = self.nups

    pltfig, axis = self.create_canvas(title, ismap)

    self.draw_map(axis, ismap)

    # save fig
    pltfig.savefig('./'+title+'.png')
    del axis


  def draw_map(self, axis, ismap):

    lon0=127.0
    lat0=38.0
    #lat0=37.30
    #lon0=127.0
    #lat0=0.0
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

    ##### draw centers #####
    x, y = map(self.lons,self.lats)
    axis.scatter(x,y,color='r',s=1)

    if self.has_nbrs:
      iup = 1
      x, y = map(self.lons[iup],self.lats[iup])
      axis.scatter(x,y,color='k',s=3.0)
      x, y = map(self.nbr_lons[iup][:],self.nbr_lats[iup][:])
      axis.scatter(x,y,color='b',s=2.0)
      x, y = map(self.cell_lons[iup][:],self.cell_lats[iup][:])
      axis.scatter(x,y,color='g',s=2.0)

  def show(self):
    plt.show()
