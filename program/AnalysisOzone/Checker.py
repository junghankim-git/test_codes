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




class Checker:

  def __init__(self, kimfilepath, kimfilenames, varname):
    #self.varname       = 'oz'
    self.varname       = varname
    self.nfiles        = len(kimfilenames)
    self.kimfilepath   = kimfilepath
    self.kimfilenames  = kimfilenames
    self.kimfiles      = ['' for i in range(self.nfiles)]

    self.lons = []
    self.lats = []

    self.ozon = [None for i in range(self.nfiles)]

    self.LoadFile()

    self.xmin =   0.0
    self.xmax = 360.0
    self.ymin = -90.0
    self.ymax =  90.0



  def LoadFile(self):

    # read atm file
    print ' * Read atm file'
    for ifile in range(self.nfiles):
      self.kimfiles[ifile] = self.kimfilepath+'/'+self.kimfilenames[ifile]
      print '  - read file: ', self.kimfiles[ifile]
      infile = nc.Dataset(self.kimfiles[ifile], mode='r')
      if ifile == 0:
        self.lons = infile.variables['lon'][:]#*180./np.pi
        self.lats = infile.variables['lat'][:]#*180./np.pi
      self.ozon[ifile] = infile.variables[self.varname][:][:][0]
      infile.close()



  def IniAxis(self, axis, xmin_in, xmax_in, ymin_in, ymax_in): #, type):
     #axis.clear()
     axis.set_xlim(xmin_in,xmax_in)
     axis.set_ylim(ymin_in,ymax_in)
     #axis.set_frame_on(False)
     #axis.set_yticklabels([])
     #axis.set_xticklabels([])
     #axis.set_yticks([])
     #axis.set_xticks([])


  def IniAxis_No(self, axis):
     axis.clear()
     axis.set_frame_on(False)
     #axis.set_xlim(xmin_in,xmax_in)
     #axis.set_ylim(ymin_in,ymax_in)
     axis.set_yticklabels([])
     axis.set_xticklabels([])
     axis.set_yticks([])
     axis.set_xticks([])



  def Draw(self, ilev=1, type=0): #, alpha, beta, lon, lat):

    type = 0

    if type==0:
      margin = 0.04     # 0.05 ~ 0.95 (for labels)
      wspace = 2.1*margin
      hspace = 0.15
    elif type==1 or type==2:
      margin = 0.06     # 0.05 ~ 0.95 (for labels)
      wspace = 0.0
      hspace = margin*1.2
    else:
      print 'unknown type..'
      quit()

    title = self.varname+' (ilev='+str(ilev)+')'

    pltfig = plt.figure(figsize=(13,12))
    pltfig.clear()
    pltfig.suptitle(title,fontsize=18,y=1-margin+0.02) # default y = 0.98
    pltfig.subplots_adjust(left=margin, bottom=margin, right=1.0-margin, top=1.0-margin, wspace=wspace, hspace=hspace)

    ismap = False

    #self.vmax = np.max(self.ozon)/100.
    self.vmax = np.max(self.ozon[1][ilev-1][:])
    if self.varname == 'oz':
      if ilev == 1:
        self.vmax = 1.16e-07
      elif ilev == 25:
        self.vmax = 1.88e-6
      elif ilev == 50:
        self.vmax = 5.54e-6
    print 'max = ', self.vmax
     

    for iplot in range(6):
      itime = iplot*2
      title = 'itime = '+str(itime)+' hour'
      axis = pltfig.add_subplot(3,2,iplot+1)
      self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
      axis.set_title(title)
      self.draw_map(axis, self.lons, self.lats, self.ozon[itime][ilev-1][:], type, ismap)
      print itime, ilev-1, np.max(self.ozon[itime][ilev-1][:])

    print ' saving...'
    oufilename = './figure_%02d.png'%(ilev)
    pltfig.savefig(oufilename)


  def draw_map(self, axis, lons, lats, var, type, isgrid):

    
    lon0 = 127.0
    lat0 = 37.30
    res  = 'c'    # resolution: c<l<i<h<f

    if type==0:
      map = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=0, urcrnrlon=360, resolution=res)
    elif type==1:
      map = Basemap(projection='ortho',lon_0=lon0,lat_0=lat0,resolution=res)                           # seoul
    elif type==2:
      map = Basemap(width=12000000,height=9000000,projection='lcc',lat_1=lat0-4.0,lat_2=lat0+4.0,lat_0=lat0,lon_0=lon0,resolution=res)
 
    # draw coastlines, country boundaries, fill continents.
    map.drawcountries(linewidth=0.8)
    map.drawcoastlines(linewidth=0.8)
    
    # draw color for land  <= not use, because it remove field's color
    #map.fillcontinents(color='coral',lake_color='aqua')
    #map.fillcontinents(color='white',lake_color='white')
    
    # draw color for ocean
    #map.drawmapboundary(fill_color='aqua')
    map.drawmapboundary(fill_color='white')
    
    # draw lat/lon grid lines every 30 degrees. (labels=[left,right,top,down], default labelstyle = '+/-')
    if type==1:
      map.drawmeridians(np.arange(0,360,45),  labels=[False,False,False,False])
      map.drawparallels(np.arange(-90,90,45), labels=[False,False,False,False])
    else:
      map.drawmeridians(np.arange(0,360,45),  labels=[False,False,False,False])
      map.drawparallels(np.arange(-90,90,45), labels=[True,False,False,False])

    
    ##### Draw map #####
    # latlon = True : use direct lon and lat
    # tri    = True : use unstructed grid data
    x, y  = map(lons, lats)
    if isgrid:
      axis.scatter(x, y, color='r', s=0.1)
    else:
      #clevs = arange(-30,30,10)
      #cs    = map.colorbar(x, y, var, tri=True)
      #cs    = map.contour(x, y, var, clevs, linewidths=0.5, colors='k', animated=True, latlon=True, tri=True)
      #cs    = map.contour(lons, lats, var, latlon=True, tri=True)
      #cs    = map.contour(x, y, var, tri=True)
  
      # Colorbar
      im = map.pcolor(x, y, var, vmin=0, vmax=self.vmax, shading='flat', cmap=plt.cm.jet, tri=True)
      #im = map.pcolor(x, y, var, vmin=0, vmax=1, shading='flat', cmap=plt.cm.jet, tri=True)
      cb = map.colorbar(im, location='bottom', format='%.2e')  # <= colorbar
      cb.ax.tick_params(labelsize=10)
      #cb.set_label("", size=0.1)
      #cb.set_label('m/s')


 

  def Show(self):
    plt.show()
