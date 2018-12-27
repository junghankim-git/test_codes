#!/usr/bin/env python
import os
import sys
import time
from numpy import *
from scipy import *
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_toolkits.basemap import Basemap



class ocn_grid:

  def __init__(self):
    self.load_file()
    self.set_plot_options()


  def load_file(self):
    self.nprocs  = 4
    self.nx      = 16
    self.ny      = 12
    self.gsize   = self.nx*self.ny

    nprocs       = self.nprocs
    nx           = self.nx
    ny           = self.ny
    print '* nprocs = ',nprocs
    print '* nx     = ',nx
    print '* ny     = ',ny
    
    # cart index,1d index,mask
    self.idx_2d = [[[i+1,j+1] for j in range(ny)] for i in range(nx)]
    self.idx_1d = [[j*nx+i+1 for j in range(ny)] for i in range(nx)]
    self.mask, self.mask_l, self.mask_s = self.gen_lsmask(nx,ny)


    # cart index,1d index,mask
    self.gsize_l, self.gsize_s, self.lsize, self.lsize_l, self.lsize_s = \
                                          self.define_sizes(nx,ny,nprocs,self.mask)
    self.rank, self.rank_l, self.rank_s = self.define_ranks(nx,ny,nprocs,self.mask)
    self.lsize_m, self.rank_m = self.define_my_domain(nx,ny,nprocs)





  def define_sizes(self,nx,ny,nprocs,mask):

    gsize   = nx*ny
    gsize_l = 0
    gsize_s = 0
    for i in range(nx):
      for j in range(ny):
        if   mask[i][j]==0:
          gsize_s = gsize_s+1
        elif mask[i][j]==1:
          gsize_l = gsize_l+1
    lsize   = [-1 for j in range(nprocs)]
    lsize_l = [-1 for j in range(nprocs)]
    lsize_s = [-1 for j in range(nprocs)]

    for iproc in range(nprocs):
      lsize[iproc],   tmp, tmp = self.decompose_1d(1,gsize  ,nprocs,iproc)
      lsize_l[iproc], tmp, tmp = self.decompose_1d(1,gsize_l,nprocs,iproc)
      lsize_s[iproc], tmp, tmp = self.decompose_1d(1,gsize_s,nprocs,iproc)

    return gsize_l, gsize_s, lsize, lsize_l, lsize_s





  def define_ranks(self,nx,ny,nprocs,mask):

    rank    = [[-1 for j in range(ny)] for i in range(nx)]
    rank_l  = [[-1 for j in range(ny)] for i in range(nx)]
    rank_s  = [[-1 for j in range(ny)] for i in range(nx)]

    iproc   = 0
    iproc_l = 0
    iproc_s = 0
    for j in range(ny):
      for i in range(nx):
        if   mask[i][j]==0:
          rank_s[i][j] = iproc_s
          iproc_s = iproc_s+1
          if iproc_s==nprocs: iproc_s=0
          rank[i][j] = rank_s[i][j]
        elif mask[i][j]==1:
          rank_l[i][j] = iproc_l
          iproc_l = iproc_l+1
          if iproc_l==nprocs: iproc_l=0
          rank[i][j] = rank_l[i][j]

    return rank, rank_l, rank_s



  def define_my_domain(self,nx,ny,nprocs):

    gsize   = nx*ny
    lsize_m = [-1 for j in range(nprocs)]
    rank_m  = [[-1 for j in range(ny)] for i in range(nx)]

    for iproc in range(nprocs):
      lsize_m[iproc], ista, iend = self.decompose_1d(1,gsize  ,nprocs,iproc)
      for i in range(ista,iend+1):
        jj = (i-1)/nx
        ii = i%nx-1
        if ii==-1: ii = nx-1
        rank_m[ii][jj] = iproc

    return lsize_m, rank_m





  def decompose_1d(self,n1,n2,nprocs,rank):

    domain = n2-n1+1
    res    = domain/nprocs
    extra  = domain%nprocs
    ista   = res*rank+n1+min(rank, extra)
    iend   = ista+res-1
    if rank<=extra-1:
      iend = iend+1
    res    = iend-ista+1

    return res, ista, iend



  def gen_lsmask(self,nx,ny):

    mask   = [[-1 for j in range(ny)] for i in range(nx)]
    mask_l = [[-1 for j in range(ny)] for i in range(nx)]
    mask_s = [[-1 for j in range(ny)] for i in range(nx)]

    fx = 0.4
    fy = 0.4
    max_land_xwidth = int(float(nx)*fx)
    max_land_ywidth = int(float(ny)*fy)
    #print max_land_xwidth,max_land_ywidth

    # left, down
    cnt = 0
    for j in range(max_land_ywidth):
      for i in range(max_land_xwidth-cnt):
        mask[i][j] = 1
      cnt = cnt+1

    # right, up
    cnt = 0
    for j in range(ny-1,max_land_ywidth,-1):
      for i in range(nx-max_land_xwidth+cnt,nx):
        mask[i][j] = 1
      cnt = cnt+1

    for i in range(nx):
      for j in range(ny):
        if mask[i][j]!=1: mask[i][j]=0
        if mask[i][j]==0:
          mask_s[i][j] = 0
        elif mask[i][j]==1:
          mask_l[i][j] = 1

    return mask, mask_l, mask_s



  def set_plot_options(self):
    nprocs = self.nprocs
    nx  = self.nx
    ny  = self.ny
    np  = 4
 
    # text size
    self.canv_size = 1.0
    if self.canv_size == 1.0:
      self.tsize  = 100./float(ny)
      self.utsize = 100./float(ny)*0.8
    else:
      self.tsize  = 160./float(ny)
      self.utsize = 160./float(ny)*0.8
     
    # view processor
    self.vproc = 1
    
    # text colors
    self.colors  = []
    self.ocolors = []
    cstep  = ['ff','aa','55','00']
    ocstep = ['99','77','44','11']
    for i in range(4):
       for j in range(3):
          self.colors.append('#'+cstep[i]+cstep[4-j-1]+cstep[4-j-1])
          self.colors.append('#'+cstep[4-j-1]+cstep[i]+cstep[4-j-1])
          self.colors.append('#'+cstep[4-j-1]+cstep[4-j-1]+cstep[i])
          self.colors.append('#'+cstep[i]+cstep[i]+cstep[4-j-1])
          self.colors.append('#'+cstep[4-j-1]+cstep[i]+cstep[i])
          self.colors.append('#'+cstep[i]+cstep[4-j-1]+cstep[i])
          self.ocolors.append('#'+ocstep[i]+ocstep[4-j-1]+ocstep[4-j-1])
          self.ocolors.append('#'+ocstep[4-j-1]+ocstep[i]+ocstep[4-j-1])
          self.ocolors.append('#'+ocstep[4-j-1]+ocstep[4-j-1]+ocstep[i])
          self.ocolors.append('#'+ocstep[i]+ocstep[i]+ocstep[4-j-1])
          self.ocolors.append('#'+ocstep[4-j-1]+ocstep[i]+ocstep[i])
          self.ocolors.append('#'+ocstep[i]+ocstep[4-j-1]+ocstep[i])
    
    # Text Postion (center)
    self.itu,self.jtu = +0.45,+0.75  # up text
    self.itd,self.jtd = +0.45,+0.45  # down text
     
    # Domain Size
    margin    = 0.2
    self.xmin = -0.1-margin
    self.xmax = nx+margin
    self.ymin = -0.1-margin
    self.ymax = ny+margin

    self.figs_dir = './figs'
  
  
  
  
  def axis_init(self,axis,xmin_in,xmax_in,ymin_in,ymax_in):
     axis.clear()
     axis.set_frame_on(False)
     axis.set_xlim(xmin_in,xmax_in)
     axis.set_ylim(ymin_in,ymax_in)
     axis.set_yticklabels([])
     axis.set_xticklabels([])
     axis.set_yticks([])
     axis.set_xticks([])


  def axis_init_no(self,axis):
     axis.clear()
     axis.set_frame_on(False)
     axis.set_yticklabels([])
     axis.set_xticklabels([])
     axis.set_yticks([])
     axis.set_xticks([])
  
  
  def coord2box(self,ii,jj,color):
     dx = 1.0
     dl = 0.9*dx
     si = ii
     ei = ii+dl
     sj = jj
     ej = jj+dl
  
     verts = [[si,sj],[si,ej],[ei,ej],[ei,sj],[si,sj]]
     codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY]
     path  = Path(verts,codes)
     return patches.PathPatch(path,facecolor=color,lw=1.0)


  def create_canvas(self,title):
    margin = 0.02
    pltfig = plt.figure(figsize=(self.nx*self.canv_size/2.0,self.ny*self.canv_size/2.0))
    pltfig.subplots_adjust(left=margin,bottom=margin,right=1.0-margin,top=1.0-margin-0.04,wspace=0.1,hspace=0.01)
    pltfig.clear()
    axis  = pltfig.add_subplot(1,1,1)
    #self.axis_init_no(axis)
    self.axis_init(axis,self.xmin,self.xmax,self.ymin,self.ymax)
    axis.set_title(title,fontsize=25)
    return pltfig,axis




  def grid_structure(self,axis,mask,rank,iscolor=False,isvproc=False):

    for i in range(self.nx):
      for j in range(self.ny):

        if mask[i][j]>=0:
          iproc = rank[i][j]
          if iscolor and mask[i][j]!=-1:
            color = self.colors[iproc]
          else:
            color = 'none'
  
          if isvproc:
             if iproc == self.vproc:
               patch = self.coord2box(i,j,color)
               axis.add_patch(patch)
          else:
             patch = self.coord2box(i,j,color)
             axis.add_patch(patch)




  def write_index_in_box(self,axis,mask,rank,index,isup=False,iscolor=False,isvproc=False):
    if isup:
      it,jt = self.itu,self.jtu
      tsize = self.utsize
      if iscolor:
        tcolor='w'
      else:
        tcolor='k'
    else:
      it,jt = self.itd,self.jtd
      tsize = self.tsize
      if iscolor:
        tcolor='k'
      else:
        tcolor='r'

    for i in range(self.nx):
      for j in range(self.ny):
        if mask[i][j]>=0:
          if isvproc:
            iproc = rank[i][j]
            if iproc == self.vproc:   # draw only viw proc
              axis.text(i+it,j+jt,str(index[i][j]),color=tcolor,size=tsize,va='center',ha='center')
          else:
            axis.text(i+it,j+jt,str(index[i][j]),color=tcolor,size=tsize,va='center',ha='center')



  def draw_all(self):
    self.draw_grid_structure('01_grid')
    '''
    self.draw_general_map('02_general_map')
    self.draw_ww3_ocean('03_ww3_ocean')
    self.draw_ww3_land('04_ww3_land')
    self.draw_ww3_all('05_ww3_all')
    self.draw_ww3_map('06_ww3_map')

    self.draw_ww3_ocean_map('07_convert',True)
    self.draw_ww3_map('08_inv_convert',True)
    '''


  def draw_grid_structure(self,title,isvproc=False):
    pltfig,axis = self.create_canvas('')
    iscolor = False
    self.grid_structure(axis,self.mask,self.rank,iscolor=iscolor,isvproc=isvproc)
    self.write_index_in_box(axis,self.mask,self.rank,self.idx_2d,isup=True,iscolor=iscolor,isvproc=isvproc)
    self.write_index_in_box(axis,self.mask,self.rank,self.idx_1d,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis



  def draw_general_map(self,title,isvproc=False):
    pltfig,axis = self.create_canvas('')
    iscolor = True
    self.grid_structure(axis,self.mask,self.rank_m,iscolor=iscolor,isvproc=isvproc)
    self.write_index_in_box(axis,self.mask,self.rank_m,self.idx_1d,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis



  def draw_ww3_ocean(self,title,isvproc=False):
    pltfig,axis = self.create_canvas('')
    iscolor = True
    self.grid_structure(axis,self.mask_s,self.rank,iscolor=iscolor,isvproc=isvproc)
    self.write_index_in_box(axis,self.mask_s,self.rank,self.rank_s,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis



  def draw_ww3_land(self,title,isvproc=False):
    pltfig,axis = self.create_canvas('')
    iscolor = True
    self.grid_structure(axis,self.mask_l,self.rank,iscolor=iscolor,isvproc=isvproc)
    self.write_index_in_box(axis,self.mask_l,self.rank,self.rank_l,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis



  def draw_ww3_all(self,title,isvproc=False):
    pltfig,axis = self.create_canvas('')
    iscolor = True
    self.grid_structure(axis,self.mask,self.rank,iscolor=iscolor,isvproc=isvproc)
    self.write_index_in_box(axis,self.mask,self.rank,self.rank,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis



  def draw_ww3_ocean_map(self,title,isvproc=False):
    pltfig,axis = self.create_canvas('')
    iscolor = True
    self.grid_structure(axis,self.mask_s,self.rank_s,iscolor=iscolor,isvproc=isvproc)
    self.write_index_in_box(axis,self.mask_s,self.rank_s,self.idx_1d,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis



  def draw_ww3_map(self,title,isvproc=False):
    pltfig,axis = self.create_canvas('')
    iscolor = True
    self.grid_structure(axis,self.mask,self.rank,iscolor=iscolor,isvproc=isvproc)
    self.write_index_in_box(axis,self.mask,self.rank,self.idx_1d,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis


  def show(self):
    plt.show()


