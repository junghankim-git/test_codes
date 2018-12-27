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

  def __init__(self, kimfilenames, wavfilename, remapfilename):
    self.varname       = 'u10m'
    self.nfiles        = len(kimfilenames)
    self.kimfilenames  = kimfilenames
    self.wavfilename   = wavfilename
    self.remapfilename = remapfilename

    self.atm_size = 0
    self.wav_size = 0

    self.atm_lons = []
    self.atm_lats = []
    self.wav_lons = []
    self.wav_lats = []

    self.atm_u10m = [0. for i in range(self.nfiles)]
    self.wav_u10m = [[] for i in range(self.nfiles-1)]

    self.LoadFile()

    self.xmin =   0.0
    self.xmax = 360.0
    self.ymin = -90.0
    self.ymax =  90.0

  def LoadFile(self):

    # read atm file
    print ' * Read atm file'
    for ifile in range(self.nfiles):
      print '  - read file: ', self.kimfilenames[ifile]
      infile = nc.Dataset(self.kimfilenames[ifile], mode='r')
      self.atm_u10m[ifile] = infile.variables[self.varname][:][0]
      infile.close()

    # read wave file
    print ' * Read wav file'
    print '  - read file: ', self.wavfilename
    infile = nc.Dataset(self.wavfilename, mode='r')
    self.wav_u10m = infile.variables['wav_u10m'][:]
    infile.close()

    # read remap file
    print ' * Read remap file'
    print '  - read file: ', self.remapfilename
    infile = nc.Dataset(self.remapfilename, mode='r')
    self.atm_size = len(infile.dimensions['src_grid_size'])
    self.wav_size = len(infile.dimensions['dst_grid_size'])
    if infile.variables['src_grid_center_lon'].units == 'radians':
      self.atm_lons = infile.variables['src_grid_center_lon'][:]*180./np.pi
      self.atm_lats = infile.variables['src_grid_center_lat'][:]*180./np.pi
      self.wav_lons = infile.variables['dst_grid_center_lon'][:]*180./np.pi
      self.wav_lats = infile.variables['dst_grid_center_lat'][:]*180./np.pi
    else:
      self.atm_lons = infile.variables['src_grid_center_lon'][:]
      self.atm_lats = infile.variables['src_grid_center_lat'][:]
      self.wav_lons = infile.variables['dst_grid_center_lon'][:]
      self.wav_lats = infile.variables['dst_grid_center_lat'][:]

    infile.close()


    if len(self.atm_u10m) != len(self.wav_u10m)+1:
      print 'check times dimension ...', len(self.atm_u10m), len(self.wav_u10m)
      quit()
    if len(self.atm_lons) != len(self.atm_u10m[0]):
      print 'check atm dimension with remap...', len(self.atm_lons), len(self.atm_u10m[0])
      quit()
    if len(self.wav_lons) != len(self.wav_u10m[0]):
      print 'check wav dimension with remap...', len(self.wav_lons), len(self.wav_u10m[0])
      quit()
    print 'dimensions.. o.k.'

    '''
    print 'atm  time = 0: ', self.atm_u10m[0], ', min = ', min(self.atm_u10m[0]), ', max = ', max(self.atm_u10m[0])
    print 'atm  time = 1: ', self.atm_u10m[1], ', min = ', min(self.atm_u10m[1]), ', max = ', max(self.atm_u10m[1])
    print 'atm  time = 2: ', self.atm_u10m[2], ', min = ', min(self.atm_u10m[2]), ', max = ', max(self.atm_u10m[2])
    print 'atm  time = 3: ', self.atm_u10m[3], ', min = ', min(self.atm_u10m[3]), ', max = ', max(self.atm_u10m[3])
    print 'wave time = 0: ', self.wav_u10m[0], ', min = ', min(self.wav_u10m[0]), ', max = ', max(self.wav_u10m[0])
    print 'wave time = 1: ', self.wav_u10m[1], ', min = ', min(self.wav_u10m[1]), ', max = ', max(self.wav_u10m[1])
    print 'wave time = 2: ', self.wav_u10m[2], ', min = ', min(self.wav_u10m[2]), ', max = ', max(self.wav_u10m[2])
    print 'wave time = 3: ', self.wav_u10m[3], ', min = ', min(self.wav_u10m[3]), ', max = ', max(self.wav_u10m[3])
    '''



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


  def CreateCanvas(self, title, type):
    plt.title = title
    pltfig = plt.figure(figsize=(10,20))
    pltfig.clear()

    axis1 = pltfig.add_subplot(2,1,1)
    self.IniAxis(axis1, self.xmin, self.xmax, self.ymin, self.ymax)

    axis2 = pltfig.add_subplot(2,1,2)
    self.IniAxis(axis2, self.xmin, self.xmax, self.ymin, self.ymax)

    return pltfig, axis1, axis2


  def Draw(self, title, type=0): #, alpha, beta, lon, lat):

    itime = 1

    if type==0:
      margin = 0.04     # 0.05 ~ 0.95 (for labels)
      wspace = 2*margin
      hspace = 0.01
    elif type==1 or type==2:
      margin = 0.06     # 0.05 ~ 0.95 (for labels)
      wspace = 0.0
      hspace = margin*1.2
    else:
      print 'unknown type..'
      quit()

   ############
   # Left plot
   ############

    pltfigL = plt.figure(figsize=(13,11))
    pltfigL.clear()
    pltfigL.suptitle(title,fontsize=18,y=1-margin+0.02) # default y = 0.98
    pltfigL.subplots_adjust(left=margin, bottom=margin, right=1.0-margin, top=1.0-margin, wspace=wspace, hspace=hspace)

    for i in range(3):
      print ' plot (left): '+str(i+1)+'-3'

      iplot = 2*i
      idata = 2*i-1

      if i == 0:
        idata = 0
        ismap = True
        title = 'grid'
      else:
        ismap = False
        title = 'itime = '+str(idata+1)+' [h]'

      axis = pltfigL.add_subplot(3,2,iplot+1)
      self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
      axis.set_title(title)
      self.draw_map(axis, self.atm_lons, self.atm_lats, self.atm_u10m[idata], type, ismap)

      axis = pltfigL.add_subplot(3,2,iplot+2)
      self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
      axis.set_title(title)
      self.draw_map(axis, self.wav_lons, self.wav_lats, self.wav_u10m[idata], type, ismap)

    print ' saving...'
    #pltfigL.savefig('./figure_'+str(itime)+'.png')
    pltfigL.savefig('./figure_1.png')

   ############
   # Right plot
   ############

    pltfigR = plt.figure(figsize=(13,11))
    pltfigR.clear()
    pltfigR.suptitle(title,fontsize=18,y=1-margin+0.02) # default y = 0.98
    pltfigR.subplots_adjust(left=margin, bottom=margin, right=1.0-margin, top=1.0-margin, wspace=wspace, hspace=hspace)

    for i in range(3):
      print ' plot (right): '+str(i+1)+'-3'

      iplot = 2*i
      idata = 2*(i+3)-1

      ismap = False
      title = 'itime = '+str(idata+1)+' [h]'

      axis = pltfigR.add_subplot(3,2,iplot+1)
      self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
      axis.set_title(title)
      self.draw_map(axis, self.atm_lons, self.atm_lats, self.atm_u10m[idata], type, ismap)

      axis = pltfigR.add_subplot(3,2,iplot+2)
      self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
      axis.set_title(title)
      self.draw_map(axis, self.wav_lons, self.wav_lats, self.wav_u10m[idata], type, ismap)

    print ' saving...'
    #pltfigR.savefig('./figure_'+str(itime)+'.png')
    pltfigR.savefig('./figure_2.png')


  def Draw_flat(self, title, type=0): #, alpha, beta, lon, lat):

    itime = 1

    if type==0:
      #margin = 0.04     # 0.05 ~ 0.95 (for labels)
      #wspace = 2*margin
      #hspace = 0.01
      margin = 0.05     # 0.05 ~ 0.95 (for labels)
      wspace = 2*margin
      hspace = 0.02
    elif type==1 or type==2:
      margin = 0.06     # 0.05 ~ 0.95 (for labels)
      wspace = 0.0
      hspace = margin*1.2
    else:
      print 'unknown type..'
      quit()

    pltfig = plt.figure(figsize=(24,6))
    pltfig.clear()
    pltfig.suptitle(title,fontsize=18,y=1-margin+0.02) # default y = 0.98
    pltfig.subplots_adjust(left=margin, bottom=margin, right=1.0-margin, top=1.0-margin, wspace=wspace, hspace=hspace)

    for i in range(4):
      print ' plot: '+str(i+1)+'-4'

      ip_atm = i+1
      ip_wav = i+5

      if i == 0:
        idata  = -1
      elif i == 1:
        idata  = 2
      elif i == 2:
        idata  = 6
      elif i == 3:
        idata  = 9
      print ' time: '+str(idata)

      if i == 0:
        ismap = True
        title = 'grid'
      else:
        ismap = False
        title = 'itime = '+str(idata)+' [h]'

      axis = pltfig.add_subplot(2,4,ip_atm)
      self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
      axis.set_title(title)
      self.draw_map(axis, self.atm_lons, self.atm_lats, self.atm_u10m[idata],   type, ismap)

      axis = pltfig.add_subplot(2,4,ip_wav)
      self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
      axis.set_title(title)
      self.draw_map(axis, self.wav_lons, self.wav_lats, self.wav_u10m[idata-1], type, ismap)

    print ' saving...'
    pltfig.savefig('./figure_flat.png')



# just 4 plots
  def Draw_grid_last(self, title, type=0): #, alpha, beta, lon, lat):

    itime = 1

    if type==0:
      margin = 0.04     # 0.05 ~ 0.95 (for labels)
      wspace = 2*margin
      hspace = 0.01
    elif type==1 or type==2:
      margin = 0.06     # 0.05 ~ 0.95 (for labels)
      wspace = 0.0
      hspace = margin*1.2
    else:
      print 'unknown type..'
      quit()

   ############
   # Left plot
   ############

    pltfigL = plt.figure(figsize=(14,8))
    pltfigL.clear()
    pltfigL.suptitle(title,fontsize=18,y=1-margin+0.02) # default y = 0.98
    pltfigL.subplots_adjust(left=margin, bottom=margin, right=1.0-margin, top=1.0-margin, wspace=wspace, hspace=hspace)


    iplot = 0
    idata = 0
    ismap = True
    title = 'Grid'

    axis = pltfigL.add_subplot(2,2,1)
    self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
    axis.set_title(title)
    self.draw_map(axis, self.atm_lons, self.atm_lats, self.atm_u10m[idata], type, ismap)

    axis = pltfigL.add_subplot(2,2,2)
    self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
    axis.set_title(title)
    self.draw_map(axis, self.wav_lons, self.wav_lats, self.wav_u10m[idata], type, ismap)

    idata = -1
    ismap = False
    title = 'u10m [10h]'

    axis = pltfigL.add_subplot(2,2,3)
    self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
    axis.set_title(title)
    self.draw_map(axis, self.atm_lons, self.atm_lats, self.atm_u10m[idata], type, ismap)

    axis = pltfigL.add_subplot(2,2,4)
    self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
    axis.set_title(title)
    self.draw_map(axis, self.wav_lons, self.wav_lats, self.wav_u10m[idata], type, ismap)

    print ' saving...'
    #pltfigR.savefig('./figure_'+str(itime)+'.png')
    pltfigL.savefig('./figure_grid_last.png')


#    for i in range(3):
#      print ' plot (left): '+str(i+1)+'-3'
#
#      iplot = 2*i
#      idata = 2*i-1
#
#      if i == 0:
#        idata = 0
#        ismap = True
#        title = 'grid'
#      else:
#        ismap = False
#        title = 'itime = '+str(idata+1)+' [h]'
#
#      axis = pltfigL.add_subplot(3,2,iplot+1)
#      self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
#      axis.set_title(title)
#      self.draw_map(axis, self.atm_lons, self.atm_lats, self.atm_u10m[idata], type, ismap)
#
#      axis = pltfigL.add_subplot(3,2,iplot+2)
#      self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
#      axis.set_title(title)
#      self.draw_map(axis, self.wav_lons, self.wav_lats, self.wav_u10m[idata], type, ismap)
#
#    print ' saving...'
#    #pltfigL.savefig('./figure_'+str(itime)+'.png')
#    pltfigL.savefig('./figure_1.png')
#
#   ############
#   # Right plot
#   ############
#
#    pltfigR = plt.figure(figsize=(13,11))
#    pltfigR.clear()
#    pltfigR.suptitle(title,fontsize=18,y=1-margin+0.02) # default y = 0.98
#    pltfigR.subplots_adjust(left=margin, bottom=margin, right=1.0-margin, top=1.0-margin, wspace=wspace, hspace=hspace)
#
#    for i in range(3):
#      print ' plot (right): '+str(i+1)+'-3'
#
#      iplot = 2*i
#      idata = 2*(i+3)-1
#
#      ismap = False
#      title = 'itime = '+str(idata+1)+' [h]'
#
#      axis = pltfigR.add_subplot(3,2,iplot+1)
#      self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
#      axis.set_title(title)
#      self.draw_map(axis, self.atm_lons, self.atm_lats, self.atm_u10m[idata], type, ismap)
#
#      axis = pltfigR.add_subplot(3,2,iplot+2)
#      self.IniAxis(axis, self.xmin, self.xmax, self.ymin, self.ymax)
#      axis.set_title(title)
#      self.draw_map(axis, self.wav_lons, self.wav_lats, self.wav_u10m[idata], type, ismap)
#
#    print ' saving...'
#    #pltfigR.savefig('./figure_'+str(itime)+'.png')
#    pltfigR.savefig('./figure_2.png')



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
      im = map.pcolor(x, y, var, vmin=-20, vmax=20, shading='flat', cmap=plt.cm.jet, tri=True)
      #cb = map.colorbar(im, location='bottom')  # <= colorbar
      #cb.set_label('m/s')


 

  def Show(self):
    plt.show()
