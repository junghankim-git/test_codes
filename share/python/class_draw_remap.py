#!/usr/bin/env python
import os
import sys
import netCDF4 as nc
sys.path.append('/home/jhkim/work/share/python')
from class_draw_unstruct import *


class class_draw_remap:



  def __init__(self, src_ll_file, src_latname, src_lonname, dst_ll_file, dst_latname, dst_lonname):

    self.src_ll_file = src_ll_file
    self.src_latname = src_latname
    self.src_lonname = src_lonname
    self.dst_ll_file = dst_ll_file
    self.dst_latname = dst_latname
    self.dst_lonname = dst_lonname

    self.src_nsize = 0
    self.src_lats  = []
    self.src_lons  = []

    self.dst_nsize = 0
    self.dst_lats  = []
    self.dst_lons  = []

    ngrids = 2

    print ' * draw_remap%draw_latlon: Initialize ...'
    self.drawer = class_draw_unstruct(ngrids)



  def set_variables(self, src_file, dst_file, varnames, ndims):

    nvars = len(varnames)
    if len(ndims) != nvars:
      print 'check dimension of varnames and ndims... ', nvars, len(ndims)

    self.nvars     = nvars
    self.varnames  = varnames
    self.src_vars  = [[] for i in range(nvars)]
    self.dst_vars  = [[] for i in range(nvars)]

    print ' * draw_remap: read variables...'
    self.src_nsize = self.drawer.read_variables_in_file(src_file,varnames,self.src_vars)
    self.dst_nsize = self.drawer.read_variables_in_file(dst_file,varnames,self.dst_vars)

    print ' * draw_remap: read lat-lon information...'
    self.src_lats,self.src_lons = \
               self.drawer.read_latlon_in_file(self.src_nsize,self.src_ll_file,self.src_latname,self.src_lonname)
    self.dst_lats,self.dst_lons = \
               self.drawer.read_latlon_in_file(self.dst_nsize,self.dst_ll_file,self.dst_latname,self.dst_lonname)




  def draw(self, title, maptype=0):

    print '  - draw_remap%draw_latlon: add grid (source)...'
    self.drawer.add_grid(self.src_lons, self.src_lats)
    print '  - draw_remap%draw_latlon: add grid (destination)...'
    self.drawer.add_grid(self.dst_lons, self.dst_lats)
    self.drawer.check_grid()

    for iv in range(self.nvars):
      print '  - draw_remap%draw_latlon: add variable :', self.varnames[iv]
      self.drawer.add_variable(self.varnames[iv], [self.src_vars[iv], self.dst_vars[iv]])

      '''
      # for hycom accuracy
      mins = [52905/2,-23.*2]; maxs = [103469.67968750000*1.3,22.164796829223633*2]
      for i in range(self.dst_nsize):
        if self.dst_vars[iv][i]<mins[iv]:
          print 'check min ',iv,i+1,i%240+1,i/240+1,self.dst_vars[iv][i]
        if self.dst_vars[iv][i]>maxs[iv]:
          print 'check max ',iv,i+1,i%240+1,i/240+1,self.dst_vars[iv][i]
      '''

    #for i in range(len(self.dst_lats)):
    #  if self.dst_vars[0][i]<-500.0:
    #    print i+1, self.dst_lons[i], self.dst_lats[i]

    print '  - draw_remap%draw_latlon: draw...'
    self.drawer.draw(title=title,maptype=maptype,figname='./figs/test.png',ongrid=True)
    self.drawer.show()



