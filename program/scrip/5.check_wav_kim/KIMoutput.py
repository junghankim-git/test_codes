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





class KIMoutput:

  def __init__(self, filename):
    self.filename = filename

    self.xmin =   0.0
    self.xmax = 360.0
    self.ymin = -90.0
    self.ymax =  90.0

  def OpenFile(self):
    # file open
    self.infile = nc.Dataset(self.filename, mode='r')

    # dimensions
    self.nhoriz = len(self.infile.dimensions['ncol'])
    self.nlevs  = len(self.infile.dimensions['lev'])

    # coordinate (lon, lat)
    self.lons, self.lats = self.GetLonLat()
    

  def GetLonLat(self):
    return self.infile.variables['lon'][:], self.infile.variables['lat'][:]

  def GetVariable(self, name):
    return self.infile.variables[name][:]
    
  def CloseFile(self):
    # close
    self.infile.close()




