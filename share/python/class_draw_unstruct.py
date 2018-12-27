#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
SHARE_DIR='/home/jhkim/work/share/python'
sys.path.append(SHARE_DIR)
from utilities import *

#------------------------------------------------------------------------------
# Class canvas utils
#------------------------------------------------------------------------------
class class_draw_unstruct:

    def __init__(self, ngrids):
        self.ne030np4file_ur = '/data/KIM3.2/cs_grid/cs_grid_ne030np4.nc'
        self.ne060np4file_ur = '/data/KIM3.2/cs_grid/cs_grid_ne060np4.nc'
        self.ne120np4file_ur = '/data/KIM3.2/cs_grid/cs_grid_ne120np4.nc'
        self.ne240np4file_ur = '/data/KIM3.2/cs_grid/cs_grid_ne240np4.nc'
        self.ne030np4file_r  = '/data/KIM3.2/cs_grid/cs_grid_ne030np4_rotated.nc'
        self.ne060np4file_r  = '/data/KIM3.2/cs_grid/cs_grid_ne060np4_rotated.nc'
        self.ne120np4file_r  = '/data/KIM3.2/cs_grid/cs_grid_ne120np4_rotated.nc'
        self.ne240np4file_r  = '/data/KIM3.2/cs_grid/cs_grid_ne240np4_rotated.nc'
        self.ne045np3file_r  = '/data/KIM3.2/cs_grid/cs_grid_ne045np3_rotated.nc'
        self.ne090np3file_r  = '/data/KIM3.2/cs_grid/cs_grid_ne090np3_rotated.nc'
        self.ne180np3file_r  = '/data/KIM3.2/cs_grid/cs_grid_ne180np3_rotated.nc'
        self.ne360np3file_r  = '/data/KIM3.2/cs_grid/cs_grid_ne360np3_rotated.nc'

        if ngrids < 1 or ngrids > 5:
            print 'Error: maximum ngrids is 4'
            quit()
        self.ngrids   = ngrids
        self.igrid    = -1
        self.size     = [0  for i in range(ngrids)] # [ngrids]
        self.lons     = [[] for i in range(ngrids)] # [ncols x ngrids]
        self.lats     = [[] for i in range(ngrids)] # [ncols x ngrids]
        self.nvars    = 0
        self.ndims    = []   # [nvars]
        self.varnames = []   # [nvars]
        self.vmin     = []   # [nvars]
        self.vmax     = []   # [nvars]
        self.vars     = [[] for i in range(ngrids)] #[nvars x ngrids]

        self.xmin     =   0.0
        self.xmax     = 360.0
        self.ymin     = -90.0
        self.ymax     =  90.0
        self.lon0     = 127.0
        self.lat0     =  37.3
        self.def_res  = 'c'    # resolution: c<l<i<h<f



    def add_grid(self, lons, lats):
        ig         = self.igrid + 1
        self.igrid = ig
        if self.igrid > self.ngrids-1:
            print 'Error: check number of grids..'
            quit()
        if len(lons) != len(lats):
            print 'Error: check the size of lons and lats...'
        self.size[ig] = len(lons)
        self.lons[ig] = lons
        self.lats[ig] = lats
  
  
  
    def check_grid(self):
        if self.igrid != self.ngrids-1:
            print 'Error: check number of grids..', self.igrid, self.ngrids
            quit()
  
  
  
  
    def add_variable(self, varname, vars):
        # check ngrids
        if len(vars) != self.ngrids:
            print 'Error: check number of vars in add_variable...'
            quit()
        # check ndims
        if self.ngrids == 1:
            ndims = self.get_ndims(vars[0])
        elif self.ngrids > 1:
            ndims = [0 for i in range(self.ngrids)]
            for ig in range(self.ngrids):
                ndims[ig] = self.get_ndims(vars[ig])
            if self.ngrids != ndims.count(ndims[0]):
                print 'Error: check ndims ', ndims
                quit()
            ndims = ndims[0]
        self.nvars = self.nvars + 1
        self.ndims.append(ndims)
        self.varnames.append(varname)
        self.vmin.append(min(vars[0]))
        self.vmax.append(max(vars[0]))
        if varname=='div_v10m':
          self.vmin[self.nvars-1] = -8e-5
          self.vmax[self.nvars-1] =  8e-5
        elif varname=='qcimps':
          self.vmin[self.nvars-1] = 0
          self.vmax[self.nvars-1] = 5e-5
        elif varname=='dx':
          self.vmin[self.nvars-1] = 10000.0
          self.vmax[self.nvars-1] = 42000.0
  
        # add variable
        for ig in range(self.ngrids):
            ncol  = self.get_dimsize(ndims, vars[ig])
            if ncol != self.size[ig]:
                print 'Error: check dimension size of variable'
                quit()
            self.vars[ig].append(vars[ig])



    def get_ndims(self, var):
        ndim = 0
        if isinstance(var, list):
            ndim = 1
        else:
            return ndim
        if isinstance(var[0], list):
            ndim = 2
        else:
            return ndim
        if isinstance(var[0][0], list):
            ndim = 3
        else:
            return ndim
        if isinstance(var[0][0][0], list):
            ndim = 4
        else:
            return ndim
        if isinstance(var[0][0][0][0], list):
            ndim = 5
        else:
            return ndim
        return ndim



    def get_dimsize(self, ndims, var):
        if ndims == 0:
            print 'Error: variable is scalar..'
            quit()
        elif ndims == 1:
            ncol = len(var[:]) 
        elif ndims == 2:
            ncol = len(var[0][:]) 
        elif ndims == 3:
            ncol = len(var[0][0][:]) 
        elif ndims == 4:
            ncol = len(var[0][0][0][:]) 
        elif ndims == 5:
            ncol = len(var[0][0][0][0][:]) 
        return ncol



    def axis_initialize(self, axis, xmin_in, xmax_in, ymin_in, ymax_in):
        axis.set_xlim(xmin_in,xmax_in)
        axis.set_ylim(ymin_in,ymax_in)



    def axis_initialize_no(self, axis):
        axis.clear()
        axis.set_frame_on(False)
        #axis.set_xlim(xmin_in,xmax_in)
        #axis.set_ylim(ymin_in,ymax_in)
        axis.set_yticklabels([])
        axis.set_xticklabels([])
        axis.set_yticks([])
        axis.set_xticks([])



    def draw(self,title,maptype,figname=None,ongrid=False,onbar=False,res=0):
        self.type = maptype
        ngrids    = self.ngrids
        nvars     = self.nvars
        if ngrids > 2:
            print 'draw api support until 2 grid...'
            quit()
        if ongrid:
            ncanv  = self.nvars+1
        else:
            ncanv  = self.nvars
        #pltfig = make_pltfig(ngrids,ncanv,5.0,0.7,0.35)
        pltfig = make_pltfig(ngrids,ncanv,7.0,0.6,0.40)
        if ongrid:
            for ig in range(ngrids):
                axis = pltfig.add_subplot(ncanv, ngrids, ig+1)
                axis.set_title('grid')
                self.draw_grid(axis, ig, res)
        # Fields
        for ig in range(ngrids):
            for iv in range(nvars):
                if ongrid:
                    axis = pltfig.add_subplot(ncanv, ngrids, ngrids+iv*ngrids+ig+1)
                else:
                    axis = pltfig.add_subplot(ncanv, ngrids, ngrids+(iv-1)*ngrids+ig+1)
                axis.set_title(self.varnames[iv])
                self.draw_field(axis, ig, iv, onbar, res)
        #plt.tight_layout()
        if figname != None:
            pltfig.savefig(figname)



    def draw_templete(self, axis, res=0):
        # c<l<i<h<f
        if   res==1:
            resol = 'c'
        elif res==2:
            resol = 'l'
        elif res==3:
            resol = 'i'
        elif res==4:
            resol = 'h'
        elif res>=5:
            resol = 'f'
        else:
            resol = self.def_res
        if self.type==0:
            map = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=0, \
                                            urcrnrlon=360, resolution=resol)
        elif self.type==1:
            map = Basemap(projection='ortho',lon_0=self.lon0,lat_0=self.lat0,resolution=resol)
        elif self.type==2:
            map = Basemap(width=12000000,height=9000000,projection='lcc',lat_1=self.lat0-4.0,\
                         lat_2=self.lat0+4.0,lat_0=self.lat0,lon_0=self.lon0,resolution=resol)
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
        if self.type==1:
            map.drawmeridians(np.arange(0,360,45),  labels=[False,False,False,False])
            map.drawparallels(np.arange(-90,90,45), labels=[False,False,False,False])
        else:
            map.drawmeridians(np.arange(0,360,45),  labels=[False,False,False,True])
            map.drawparallels(np.arange(-90,90,45), labels=[True,False,False,False])
        return map



    def draw_grid(self, axis, ig, res=0):
        map   = self.draw_templete(axis,res)
        x, y  = map(self.lons[ig], self.lats[ig])
        axis.scatter(x, y, color='r', s=0.1)



    def draw_field(self, axis, ig, ivar, onbar=False, res=0):
        map  = self.draw_templete(axis,res)
        if False:
           x, y = map(self.lons[ig], self.lats[ig])  # lats, lons must be numpy array
           im   = map.pcolor(x,y,self.vars[ig][ivar],vmin=self.vmin[ivar],vmax=self.vmax[ivar],\
                             shading='flat',cmap=plt.cm.jet,tri=True)
        else:
           im   = map.contourf(self.lons[ig],self.lats[ig],self.vars[ig][ivar],vmin=self.vmin[ivar],vmax=self.vmax[ivar],\
                               cmap=plt.cm.jet,tri=True)
        if onbar:
            #cb = map.colorbar(im, location='bottom')
            cb = map.colorbar(im,location='right',size='1.8%')
            #cb.set_label('?')
        axis.text(0  ,95,'min = %7.3e'%(min(self.vars[ig][ivar])),fontsize=12,ha='left', va='center')
        axis.text(360,95,'max = %7.3e'%(max(self.vars[ig][ivar])),fontsize=12,ha='right',va='center')



    def show(self):
        plt.show()



    def read_variable_in_file(self, filename, varname):
        size = 0
        # file open
        infile = nc.Dataset(filename, mode='r')
        ndiv = varname.count(':')
        #
        if (ndiv > 4):
            print 'Error: ndiv must less than 5.'
            quit()
        i2 = 0
        i3 = 0
        i4 = 0
        i5 = 0
        varn = varname
        ii = [-1 for i in range(ndiv)]
        np1 = -1
        nm1 = -1
        for id in range(ndiv):
            np1 = varname.index(':',nm1+1)
            if   id == 0 and id == ndiv-1:
                varn = varname[nm1+1:np1]
                ii[ndiv-1] = int(varname[np1+1:])
            elif id == 0:
                varn = varname[nm1+1:np1]
            elif id == ndiv-1:
                ii[ndiv-2] = int(varname[nm1+1:np1])
                ii[ndiv-1] = int(varname[np1+1:])
            else:
                ii[id-1] = int(varname[nm1+1:np1])
            nm1 = np1
        if   ndiv == 1:
            i2 = ii[0]-1
        elif ndiv == 2:
            i2 = ii[0]-1
            i3 = ii[1]-1
        elif ndiv == 3:
            i2 = ii[0]-1
            i3 = ii[1]-1
            i4 = ii[2]-1
        elif ndiv == 4:
            i2 = ii[0]-1
            i3 = ii[1]-1
            i4 = ii[2]-1
            i5 = ii[3]-1
        tmp = infile.variables[varn][:].tolist()
        if len(tmp) == 1:
            tmp = tmp[0][:]
        ndims = self.get_ndims(tmp)
        if   ndims == 1:
            var = tmp[:]
            print '  - read: ', varn
        elif ndims == 2:
            var = tmp[i2][:]
            print '  - read: ', varn, i2
        elif ndims == 3:
            var = tmp[i3][i2][:]
            print '  - read: ', varn, i2, i3
        elif ndims == 4:
            var = tmp[i4][i3][i2][:]
            print '  - read: ', varn, i2, i3, i4
        elif ndims == 5:
            var = tmp[i5][i4][i3][i2][:]
            print '  - read: ', varn, i2, i3, i4, i5
        size = len(var)
        del tmp    
        print '  - min, max: ', min(var), max(var)

        # close
        infile.close()
        return size, var



    def read_variables_in_file(self, filename, varnames, vars):
        nvars    = len(varnames)
        if len(vars) != nvars:
            print 'Error: check dimensions of arguments... in read_variables_in_file.. ', nvars, len(vars)
            quit()
        sizes = [0 for i in range(nvars)]
        # file open
        infile = nc.Dataset(filename, mode='r')
        for ivar in range(nvars):
            sizes[ivar], vars[ivar] = self.read_variable_in_file(filename, varnames[ivar])
            #print vars[ivar][0], vars[ivar][10]
        if sizes.count(sizes[0]) != nvars:
            print 'Error: check the horizontal size of variables...'
            quit()
        # close
        infile.close()
        return sizes[0]



    def read_variables_in_file_org(self, filename, varnames, vars):
        nvars    = len(varnames)
        if len(vars) != nvars:
            print 'Error: check dimensions of arguments... in read_variables_in_file.. ', nvars, len(vars)
            quit()
        sizes = [0 for i in range(nvars)]
        # file open
        infile = nc.Dataset(filename, mode='r')
        for ivar in range(nvars):
            ndiv = varnames[ivar].count(':')
            if (ndiv > 4):
                print 'Error: ndiv must less than 5.'
                quit()
            i2 = 0
            i3 = 0
            i4 = 0
            i5 = 0
            varname = varnames[ivar]
            ii = [-1 for i in range(ndiv)]
            np1 = -1
            nm1 = -1
            for id in range(ndiv):
                np1 = varnames[ivar].index(':',nm1+1)
                if   id == 0 and id == ndiv-1:
                    varname = varnames[ivar][nm1+1:np1]
                    ii[ndiv-1] = int(varnames[ivar][np1+1:])
                elif id == 0:
                    varname = varnames[ivar][nm1+1:np1]
                elif id == ndiv-1:
                    ii[ndiv-2] = int(varnames[ivar][nm1+1:np1])
                    ii[ndiv-1] = int(varnames[ivar][np1+1:])
                else:
                    ii[id-1] = int(varnames[ivar][nm1+1:np1])
                nm1 = np1
            #varnames[ivar] = varname
            if   ndiv == 1:
                i2 = ii[0]-1
            elif ndiv == 2:
                i2 = ii[0]-1
                i3 = ii[1]-1
            elif ndiv == 3:
                i2 = ii[0]-1
                i3 = ii[1]-1
                i4 = ii[2]-1
            elif ndiv == 4:
                i2 = ii[0]-1
                i3 = ii[1]-1
                i4 = ii[2]-1
                i5 = ii[3]-1
            #tmp = infile.variables[varnames[ivar]][:].tolist()
            tmp = infile.variables[varname][:].tolist()
            if len(tmp) == 1:
                tmp = tmp[0][:]
            ndims = self.get_ndims(tmp)
            if   ndims == 1:
                vars[ivar] = tmp[:]
                print '  - read: ', varname
            elif ndims == 2:
                vars[ivar] = tmp[i2][:]
                print '  - read: ', varname, i2
            elif ndims == 3:
                vars[ivar] = tmp[i3][i2][:]
                print '  - read: ', varname, i2, i3
            elif ndims == 4:
                vars[ivar] = tmp[i4][i3][i2][:]
                print '  - read: ', varname, i2, i3, i4
            elif ndims == 5:
                vars[ivar] = tmp[i5][i4][i3][i2][:]
                print '  - read: ', varname, i2, i3, i4, i5
            sizes[ivar] = len(vars[ivar])
            del tmp    
            print '  - min, max: ', min(vars[ivar]), max(vars[ivar])
        if sizes.count(sizes[0]) != nvars:
            print 'Error: check the horizontal size of variables...'
            quit()
        # close
        infile.close()
        return sizes[0]



    def read_latlon_in_file(self, nqp, ncols, filename, latname, lonname):
        # file open
        infile = nc.Dataset(filename, mode='r')
        # processing None
        if latname==None: latname='lats'
        if lonname==None: lonname='lons'
        if infile.variables.keys().count(latname)<1 or infile.variables.keys().count(lonname)<1:
            infile.close()
            if nqp==3:
              if   ncols ==   48602:
                  infile = nc.Dataset(self.ne045np3file_r, mode='r')
              elif ncols ==  194402:
                  infile = nc.Dataset(self.ne090np3file_r, mode='r')
              elif ncols ==  777602:
                  infile = nc.Dataset(self.ne180np3file_r, mode='r')
              elif ncols == 3110402:
                  infile = nc.Dataset(self.ne360np3file_r, mode='r')
              else:
                  print 'Error: It is an unsupported size, ncols = ', ncols
                  quit()
            elif nqp==4:
              if   ncols ==   48602:
                  infile = nc.Dataset(self.ne030np4file_r, mode='r')
              elif ncols ==  194402:
                  infile = nc.Dataset(self.ne060np4file_r, mode='r')
              elif ncols ==  777602:
                  infile = nc.Dataset(self.ne120np4file_r, mode='r')
              elif ncols == 3110402:
                  infile = nc.Dataset(self.ne240np4file_r, mode='r')
              else:
                  print 'Error: It is an unsupported size, ncols = ', ncols
                  quit()
            else:
              print 'check np = ',nqp
              quit()
            print ' * warn: can not find lat, lon dimension in file, it read from the default llfile...'
            latlons = infile.variables['latlons'][:]
            lats = np.array([0.0 for i in range(ncols)])
            lons = np.array([0.0 for i in range(ncols)])
            for i in range(ncols):
                lats[i] = latlons[i][0]
                lons[i] = latlons[i][1]
        else:
            # variables
            lats = infile.variables[latname][:]#.tolist() ! lats, lons must be numpy array
            lons = infile.variables[lonname][:]#.tolist() ! lats, lons must be numpy array
        nsize = len(lats)
        if nsize != ncols:
            print 'Error: the size of the coordinate is different from the size of the variables', nsize, ncols
            quit()
        if max(lons) < 2.0*np.pi+1e-6: # if radian
            for i in range(nsize):
                lats[i] = lats[i]*180./np.pi
                lons[i] = lons[i]*180./np.pi
        # close
        infile.close()
        return lats, lons

#------------------------------------------------------------------------------
#
#
#
#
#
#
#------------------------------------------------------------------------------
# Class class_draw_interface
#------------------------------------------------------------------------------

class class_draw_interface:
      def __init__(self,infiles,llfiles,nqp):
          self.times    = [0.0 for i in range(7)]
          self.times[0] = time.time()
          self.times[1] = time.time()
          self.ngrids   = len(infiles)
          self.infiles  = infiles
          self.nqp      = nqp
          self.llfiles  = ['' for i in range(self.ngrids)]
          for ig in range(self.ngrids):
              if not os.path.isfile(self.infiles[ig]):
                  print 'Error: check file name.... ', self.infiles[ig]
                  quit()
          for ig in range(self.ngrids):
              if llfiles[ig]==None:
                  self.llfiles[ig] = infiles[ig]
              else:
                  self.llfiles[ig] = llfiles[ig]
          self.ncols = [ 0 for i in range(self.ngrids)]
          self.lats  = [[] for i in range(self.ngrids)]
          self.lons  = [[] for i in range(self.ngrids)]
          print ' * log: initialize ...'
          self.canvas = class_draw_unstruct(self.ngrids)
          self.times[1] = time.time()-self.times[1]



      def set_variables(self, varnames):
          self.times[2] = time.time()
          self.varnames = varnames
          self.nvars    = len(varnames)
          self.vars     = [[[] for i in range(self.nvars)] for i in range(self.ngrids)]
          print ' * log: read variables...'
          for ig in range(self.ngrids):
              self.ncols[ig] = self.canvas.read_variables_in_file(self.infiles[ig],self.varnames,self.vars[ig])
          self.times[2] = time.time()-self.times[2]



      def set_coordinates(self, latnames, lonnames):
          self.times[3] = time.time()
          print ' * log: read lat-lon information...'
          for ig in range(self.ngrids):
              self.lats[ig],self.lons[ig] = self.canvas.read_latlon_in_file(self.nqp,self.ncols[ig],self.llfiles[ig],latnames[ig],lonnames[ig])
          self.times[3] = time.time()-self.times[3]



      def draw(self, title, maptype=0, figname=None, ongrid=False, onbar=False, res=0):
          self.times[4] = time.time()
          print ' * log: add grids...'
          for ig in range(self.ngrids):
              self.canvas.add_grid(self.lons[ig], self.lats[ig])
          self.canvas.check_grid()
          self.times[4] = time.time()-self.times[4]
          self.times[5] = time.time()
          for iv in range(self.nvars):
              print ' * log: add variables: ', self.varnames[iv]
              vars = []
              for ig in range(self.ngrids):
                  vars.append(self.vars[ig][iv][:])
              self.canvas.add_variable(self.varnames[iv],vars)
              del vars
          self.times[5] = time.time()-self.times[5]
          self.times[6] = time.time()
          print ' * log: draw...'
          self.canvas.draw(title=title,maptype=maptype,figname=figname,ongrid=ongrid,onbar=onbar,res=res)
          self.times[6] = time.time()-self.times[6]
          self.times[0] = time.time()-self.times[0]
          print ' * log: elapse time = %5.2f sec'%(self.times[0])
          print '  - initialize time = %5.2f sec'%(self.times[1])
          print '  - read vars  time = %5.2f sec'%(self.times[2])
          print '  - read ll.   time = %5.2f sec'%(self.times[3])
          print '  - add grid   time = %5.2f sec'%(self.times[4])
          print '  - add vars.  time = %5.2f sec'%(self.times[5])
          print '  - draw       time = %5.2f sec'%(self.times[6])
          if figname==None:
              self.canvas.show()
          else:
              print ' * log: saving file...: ', figname
