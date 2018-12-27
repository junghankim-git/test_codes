#!/usr/bin/env python
import os
import sys
import time
from numpy import *
from scipy import *
import netCDF4 as nc
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_toolkits.basemap import Basemap
#from matplotlib import cm


pltfig = plt.figure(figsize=(10,8))
pltfig.clear()
axis  = pltfig.add_subplot(1,1,1)

map = Basemap(projection='ortho',lon_0=127.0,lat_0=37.3,resolution='l')

map.drawcoastlines(linewidth=0.8)
map.drawcountries(linewidth=0.8)

map.drawgreatcircle(0,0,90,2,linewidth=2,color='k')

plt.show()
