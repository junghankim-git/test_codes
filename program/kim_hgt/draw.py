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

    def __init__(self, ngrids, new=False, level=1):
        self.ne030file_ur = '/data/KIM2.3/cs_grid/cs_grid_ne030np4.nc'
        self.ne060file_ur = '/data/KIM2.3/cs_grid/cs_grid_ne060np4.nc'
        self.ne120file_ur = '/data/KIM2.3/cs_grid/cs_grid_ne120np4.nc'
        self.ne240file_ur = '/data/KIM2.3/cs_grid/cs_grid_ne240np4.nc'
        self.ne030file_r  = '/data/KIM2.3/cs_grid/cs_grid_ne030np4_rotated.nc'
        self.ne060file_r  = '/data/KIM2.3/cs_grid/cs_grid_ne060np4_rotated.nc'
        self.ne120file_r  = '/data/KIM2.3/cs_grid/cs_grid_ne120np4_rotated.nc'
        self.ne240file_r  = '/data/KIM2.3/cs_grid/cs_grid_ne240np4_rotated.nc'

        self.new      = new
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
        vmin = np.amin(vars[1])
        vmax = np.amax(vars[1])
        self.vmin.append(vmin)
        self.vmax.append(vmax)
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
        if ngrids > 3:
            print 'draw api support until 2 grid...'
            quit()
        if ongrid:
            ncanv  = self.nvars+1
        else:
            ncanv  = self.nvars
        #pltfig = make_pltfig(ngrids,ncanv,5.0,0.7,0.35)
        #pltfig = make_pltfig(ngrids,ncanv,5.0,0.6,0.40)
        pltfig = make_pltfig(ngrids,ncanv,base=5.0,mrg_wgt=0.5,yscale=0.45)
        if ongrid:
            for ig in range(ngrids):
                axis = pltfig.add_subplot(ncanv, ngrids, ig+1)
                axis.set_title('grid')
                self.draw_grid(axis, ig, res)
        # Fields
        if not self.new:
            names = ['(dyn,0step)','(phy,1step)','(1step-0step)']
        else:
            names = ['(phy,0step)','(phy,1step)','(1step-0step)']
        for ig in range(ngrids):
            for iv in range(nvars):
                if ongrid:
                    axis = pltfig.add_subplot(ncanv, ngrids, ngrids+iv*ngrids+ig+1)
                else:
                    axis = pltfig.add_subplot(ncanv, ngrids, ngrids+(iv-1)*ngrids+ig+1)
                axis.set_title(self.varnames[iv]+names[ig])
                #print('draw: var {},{}'.format(ig,iv))
                self.draw_field(axis, ig, iv, onbar, res)
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
            map = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=0, urcrnrlon=360, resolution=resol)
        elif self.type==1:
            map = Basemap(projection='ortho',lon_0=self.lon0,lat_0=self.lat0,resolution=resol)
        elif self.type==2:
            map = Basemap(width=12000000,height=9000000,projection='lcc',lat_1=self.lat0-4.0,lat_2=self.lat0+4.0,lat_0=self.lat0,lon_0=self.lon0,resolution=resol)
        map.drawcountries(linewidth=0.8)
        map.drawcoastlines(linewidth=0.8)
        map.drawmapboundary(fill_color='white')
        if self.type==1:
            map.drawmeridians(np.arange(0,360,45),  labels=[False,False,False,False])
            map.drawparallels(np.arange(-90,90,45), labels=[False,False,False,False])
        else:
            map.drawmeridians(np.arange(0,360,45),  labels=[False,False,False,True])
            map.drawparallels(np.arange(-90,90,45), labels=[True,False,False,False])
        return map



    def draw_grid(self, axis, ig, res=0):
        map  = self.draw_templete(axis,res)
        x, y = map(self.lons[ig], self.lats[ig])
        axis.scatter(x, y, color='r', s=0.1)



    def draw_field(self, axis, ig, ivar, onbar=False, res=0):
        map  = self.draw_templete(axis,res)
        x, y = map(self.lons[ig], self.lats[ig])  # lats, lons must be numpy array
        if ig!=2:
            vmin = self.vmin[ivar]
            vmax = self.vmax[ivar]
        else:
            #vmin = np.amin(self.vars[ig][ivar])*0.5
            #vmax = np.amax(self.vars[ig][ivar])*0.5
            vmin =-1000.0
            vmax = 1000.0
            if   ivar==0:
                vmin = -0.5
                vmax =  0.5
            elif ivar==1:
                vmin = -150.
                vmax =  150.
            elif ivar==2:
                vmin = -800.
                vmax =  800.
            elif ivar==3:
                vmin = -1500.
                vmax =  1500.
        #print 'vmin = {}, vmax = {}'.format(vmin,vmax)
        #print 'vmin = {}, vmax = {}'.format(np.amin(self.vars[ig][ivar]),np.amax(self.vars[ig][ivar]))
        im = map.pcolor(x,y,self.vars[ig][ivar],vmin=vmin,vmax=vmax,shading='flat',cmap=plt.cm.jet,tri=True)
        if onbar:
            cb = map.colorbar(im,location='right',size='1.8%')
        axis.text(0  ,95,'min = %7.3e'%(np.amin(self.vars[ig][ivar])),fontsize=12,ha='left', va='center')
        axis.text(360,95,'max = %7.3e'%(np.amax(self.vars[ig][ivar])),fontsize=12,ha='right',va='center')



    def show(self):
        plt.show()



    def read_variables_in_file(self, filename, varnames, vars):
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



    def read_latlon_in_file(self, ncols, filename, latname, lonname):
        # file open
        infile = nc.Dataset(filename, mode='r')
        # processing None
        if latname==None: latname='lats'
        if lonname==None: lonname='lons'
        if infile.variables.keys().count(latname)<1 or infile.variables.keys().count(lonname)<1:
            infile.close()
            if   ncols ==   48602:
                infile = nc.Dataset(self.ne030file_r, mode='r')
            elif ncols ==  194402:
                infile = nc.Dataset(self.ne060file_r, mode='r')
            elif ncols ==  777602:
                infile = nc.Dataset(self.ne120file_r, mode='r')
            elif ncols == 3110402:
                infile = nc.Dataset(self.ne240file_r, mode='r')
            else:
                print 'Error: It is an unsupported size, ncols = ', ncols
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



#new = False
new = True
level = [1,25,49,50]
#level = [50]
nlevs = len(level)

if new:
    def_path = '/data/jhkim/TestBed/KIM/Output/3.0.10/ne30/gnu/10h/SW_G/101_new'
    infiles = ['/scratch/jhkim/TestBed/Data/3.0.09.04/10h/SW_G/101/UP-20110725120000-000000.nc',def_path+'/UP-20110725130000-000001.nc']
else:
    def_path = '/data/jhkim/TestBed/KIM/Output/3.0.10/ne30/gnu/10h/SW_G/101_org'
    infiles = [def_path+'/UP-20110725120000-000000.nc',def_path+'/UP-20110725130000-000001.nc']

if new:
   fig_name = 'phy.png'
else:
   fig_name = 'dyn.png'


ngrids = 3
varnames = ['hgt:'+str(lev) for lev in level]
nvars  = len(level)
print('ngrids = {}, nvars = {}'.format(ngrids,nvars))
vars = [[[] for j in range(nvars)] for i in range(ngrids)]
canvas = class_draw_unstruct(ngrids,new,1)


print('read variables')
print(' - read file :',infiles)
for iv in range(nvars):
   ncols = canvas.read_variables_in_file(infiles[0],varnames,vars[0])
   ncols = canvas.read_variables_in_file(infiles[1],varnames,vars[1])
   vars[2][iv] = [0.0 for i in range(ncols)]
   for i in range(ncols):
      vars[2][iv][i] = vars[1][iv][i]-vars[0][iv][i]
   print('vmin = {}, vmax = {}'.format(np.amin(vars[2][iv]),np.amax(vars[2][iv])))
lats,lons = canvas.read_latlon_in_file(ncols,infiles[0],'lats','lons')
canvas.add_grid(lons,lats)
canvas.add_grid(lons,lats)
canvas.add_grid(lons,lats)
canvas.check_grid()
for iv in range(nvars):
   gvars = []
   for ig in range(ngrids):
      gvars.append(vars[ig][iv][:])
   canvas.add_variable(varnames[iv],gvars)
   del gvars

#canvas.add_variable('hgt',vars)
print('draw')
canvas.draw(title='test',maptype=0,figname=fig_name,onbar=True)
print('save fig: '+fig_name)
canvas.show()


