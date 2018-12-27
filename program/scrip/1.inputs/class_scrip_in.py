#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import scipy as sp
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_toolkits.basemap import Basemap
sys.path.append('/home/jhkim/work/share/python')
from utilities import *


class class_remap_in:

  def __init__(self, filename):
    self.filename = filename
    self.read_file_scrip_in()

    self.xmin =   0.0
    self.xmax = 360.0
    self.ymin = -90.0
    self.ymax =  90.0


  def read_file_scrip_in(self):
    # file open
    self.infile       = nc.Dataset(self.filename,mode='r')

    # dimensions
    self.gsize   = len(self.infile.dimensions['grid_size'])
    self.csize   = len(self.infile.dimensions['grid_corners'])
    gsize        = self.gsize
    csize        = self.csize
    
    # allocate
    self.lats      = [0.0 for i in range(gsize)]
    self.lons      = [0.0 for i in range(gsize)]
    self.lats_nbrs = [[0.0 for j in range(csize)] for i in range(gsize)]
    self.lons_nbrs = [[0.0 for j in range(csize)] for i in range(gsize)]

    # variables
    #if self.infile.variables['gr'].units == 'radians':
    self.lats      = self.infile.variables['grid_center_lat'][:]*180./np.pi
    self.lons      = self.infile.variables['grid_center_lon'][:]*180./np.pi
    self.lats_nbrs = self.infile.variables['grid_corner_lat'][:][:]*180./np.pi
    self.lons_nbrs = self.infile.variables['grid_corner_lon'][:][:]*180./np.pi
    
    # close
    self.infile.close()


  def axis_init(self, axis, xmin_in, xmax_in, ymin_in, ymax_in): #, type):
     axis.set_xlim(xmin_in,xmax_in)
     axis.set_ylim(ymin_in,ymax_in)


  def axis_init_no(self, axis):
     axis.clear()
     axis.set_frame_on(False)
     #axis.set_xlim(xmin_in,xmax_in)
     #axis.set_ylim(ymin_in,ymax_in)
     axis.set_yticklabels([])
     axis.set_xticklabels([])
     axis.set_yticks([])
     axis.set_xticks([])


  def draw(self, title, ix=None): #, alpha, beta, lon, lat):

    '''
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
    '''
    margin = 0.04     # 0.05 ~ 0.95 (for labels)
    wspace = 2*margin
    hspace = 0.01
#
    pltfig = plt.figure(figsize=(16,5))
    pltfig.clear()
    pltfig.suptitle(title,fontsize=18,y=1-margin+0.02) # default y = 0.98
    pltfig.subplots_adjust(left=margin,bottom=margin,right=1.0-margin,top=1.0-margin,wspace=wspace,hspace=hspace)

    print ' plot 1'
    axis1 = pltfig.add_subplot(1,2,1)
    self.axis_init(axis1,self.xmin,self.xmax,self.ymin,self.ymax)
    axis1.set_title('grid')
    self.draw_map(axis1,0,ix)
    axis1 = pltfig.add_subplot(1,2,2)
    self.axis_init(axis1,self.xmin,self.xmax,self.ymin,self.ymax)
    axis1.set_title('grid')
    self.draw_map(axis1,1,ix)

    # save fig
    print ' saving...'
    pltfig.savefig('./figure.png')

#def cart_to_idx(nxx,nyy,ixx,iyy):
#    return ixx+(iyy-1)*nxx

  def draw_corners(self, title, si=None, ei=None): #, alpha, beta, lon, lat):
    if si==None: si = 0
    if ei==None: ei = self.gsize
    pltfig,axis = create_plot(1,1)
    for i in range(si-1,ei):
        x = i%240
        y = i/240
    #for x in range(240):
    #  for y in [175]:
        i = x+y*240
        print('i = {}, x = {}, y = {}'.format(i+1,x+1,y+1))
        x0 = self.lons[i]
        y0 = self.lats[i]
        xs = self.lons_nbrs[i]
        ys = self.lats_nbrs[i]
        xmin = np.amin(xs)
        xmax = np.amax(xs)
        margin = (xmax-xmin)*0.1
        xmin = xmin-margin
        xmax = xmax+margin
        ymin = np.amin(ys)
        ymax = np.amax(ys)
        margin = (ymax-ymin)*0.1
        ymin = ymin-margin
        ymax = ymax+margin
        axis.clear()
        axis_init(axis,'{} ({},{})'.format(i+1,x+1,y+1),'lon','lat',xmin,xmax,ymin,ymax)
        self.draw_corner(axis,x0,y0,xs,ys)
        pltfig.savefig('./figs/{:06d}.png'.format(i+1))


  def draw_corner(self, axis, lon0, lat0, lons, lats):
    axis.scatter(lon0,lat0,s=50,c='k',marker='o',linewidth=0)
    axis.text(lon0,lat0,'center',color='k',ha='center',va='bottom')
    nbrs = len(lons)
    colors = ['r','b','g','c','m','y']
    for i in range(nbrs):
        x0,y0 = lons[i],lats[i]
        axis.scatter(x0,y0,s=50,c=colors[i],marker='o',linewidth=0)
        #axis.text(x0,y0,'{}'.format(i+1),color=colors[i],ha='center',va='bottom')
        axis.text(x0,y0,'{}'.format(i+1),color='k',ha='center',va='bottom')

        if i==nbrs-1:
            x1,y1 = lons[0],lats[0]
        else:
            x1,y1 = lons[i+1],lats[i+1]
        if x0!=x1 or y0!=y1:
            width  = 0.000001
            head_w = np.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))*0.001
            head_l = np.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))*0.02
            axis.arrow(x0,y0,x1-x0,y1-y0,fc='y',ec='y',alpha=0.5,width=0.000000,length_includes_head=True)#,head_width=head_w,head_length=head_l)


  def draw_map(self, axis, type, ix=None): #, lons, lats):

    res  = 'c'    # resolution: c<l<i<h<f
    if ix==None:
        lon0 = 127.0
        lat0 = 37.30
    else:
        lon0 = self.lons[ix-1]
        lat0 = self.lats[ix-1]

    if type==0:
        map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360,resolution=res)
    elif type==1:
        map = Basemap(projection='ortho',lon_0=lon0,lat_0=lat0,resolution=res)                           # seoul
    elif type==2:
        map = Basemap(width=12000000,height=9000000,projection='lcc',lat_1=lat0-4.0,lat_2=lat0+4.0,lat_0=lat0,lon_0=lon0,resolution=res)
 
    # draw coastlines,country boundaries,fill continents.
    map.drawcountries(linewidth=0.8)
    map.drawcoastlines(linewidth=0.8)
    
    # draw color for land  <= not use,because it remove field's color
    #map.fillcontinents(color='coral',lake_color='aqua')
    #map.fillcontinents(color='white',lake_color='white')
    
    # draw color for ocean
    #map.drawmapboundary(fill_color='aqua')
    map.drawmapboundary(fill_color='white')
    
    # draw lat/lon grid lines every 30 degrees. (labels=[left,right,top,down],default labelstyle = '+/-')
    if type==1:
        map.drawmeridians(np.arange(0,360,45), labels=[False,False,False,False])
        map.drawparallels(np.arange(-90,90,45),labels=[False,False,False,False])
    else:
        map.drawmeridians(np.arange(0,360,45), labels=[False,False,False,True])
        map.drawparallels(np.arange(-90,90,45),labels=[True,False,False,False])

    ##### Draw map #####
    if ix==None:
        x,y = map(self.lons,self.lats)
        axis.scatter(x,y,color='r',s=0.1)
    else:
        x,y = map(self.lons[ix-1],self.lats[ix-1])
        axis.scatter(x,y,color='r',s=0.1)
        x,y = map(self.lons_nbrs[ix-1],self.lats_nbrs[ix-1])
        axis.scatter(x,y,color='b',s=0.1)



  def draw_cart(self, axis):

    ##### Draw centers #####
    x, y = self.center_lons, self.center_lats
    axis.scatter(x,y,color='r',s=1)

    ##### Draw corners #####
    #for i in range(self.gsize):
    #  x, y = self.corner_lons[i][:], self.corner_lats[i][:]
    #  axis.scatter(x,y,s=1)
    #  axis.scatter(x,y,s=1)


 

  def show(self):
    plt.show()
