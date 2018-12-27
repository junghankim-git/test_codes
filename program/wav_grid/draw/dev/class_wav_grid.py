import sys

sys.path.append('/home/jhkim/work/share/python')
from class_plot_cart import *
from utilities       import decompose_1d



class wav_grid:

  def __init__(self):

    self.nprocs  = 4
    self.nx      = 15
    self.ny      = 9
    self.gsize   = self.nx*self.ny

    nprocs       = self.nprocs
    nx           = self.nx
    ny           = self.ny

    self.figs_dir = './figs'
    self.cart     = plot_cart(nx,ny,nprocs)
    
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
      lsize[iproc],   tmp, tmp = decompose_1d(1,gsize  ,nprocs,iproc)
      lsize_l[iproc], tmp, tmp = decompose_1d(1,gsize_l,nprocs,iproc)
      lsize_s[iproc], tmp, tmp = decompose_1d(1,gsize_s,nprocs,iproc)

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



  def draw_all(self):
    self.draw_grid_structure('01_grid')
    self.draw_general_map('02_general_map')
    self.draw_ww3_ocean('03_ww3_ocean')
    self.draw_ww3_land('04_ww3_land')
    self.draw_ww3_all('05_ww3_all')
    self.draw_ww3_map('06_ww3_map')

    self.draw_ww3_ocean_map('07_convert',True)
    self.draw_ww3_map('08_inv_convert',True)


  def draw_grid_structure(self,title,isvproc=False):
    pltfig,axis = self.cart.create_canvas('')
    iscolor = False
    self.cart.grid_structure(axis,self.mask,self.rank,iscolor=iscolor,isvproc=isvproc)
    self.cart.write_index_in_box(axis,self.mask,self.rank,self.idx_2d,isup=True,iscolor=iscolor,isvproc=isvproc)
    self.cart.write_index_in_box(axis,self.mask,self.rank,self.idx_1d,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis



  def draw_general_map(self,title,isvproc=False):
    pltfig,axis = self.cart.create_canvas('')
    iscolor = True
    self.cart.grid_structure(axis,self.mask,self.rank_m,iscolor=iscolor,isvproc=isvproc)
    self.cart.write_index_in_box(axis,self.mask,self.rank_m,self.idx_1d,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis



  def draw_ww3_ocean(self,title,isvproc=False):
    pltfig,axis = self.cart.create_canvas('')
    iscolor = True
    self.cart.grid_structure(axis,self.mask_s,self.rank,iscolor=iscolor,isvproc=isvproc)
    self.cart.write_index_in_box(axis,self.mask_s,self.rank,self.rank_s,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis



  def draw_ww3_land(self,title,isvproc=False):
    pltfig,axis = self.cart.create_canvas('')
    iscolor = True
    self.cart.grid_structure(axis,self.mask_l,self.rank,iscolor=iscolor,isvproc=isvproc)
    self.cart.write_index_in_box(axis,self.mask_l,self.rank,self.rank_l,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis



  def draw_ww3_all(self,title,isvproc=False):
    pltfig,axis = self.cart.create_canvas('')
    iscolor = True
    self.cart.grid_structure(axis,self.mask,self.rank,iscolor=iscolor,isvproc=isvproc)
    self.cart.write_index_in_box(axis,self.mask,self.rank,self.rank,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis



  def draw_ww3_ocean_map(self,title,isvproc=False):
    pltfig,axis = self.cart.create_canvas('')
    iscolor = True
    self.cart.grid_structure(axis,self.mask_s,self.rank_s,iscolor=iscolor,isvproc=isvproc)
    self.cart.write_index_in_box(axis,self.mask_s,self.rank_s,self.idx_1d,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis



  def draw_ww3_map(self,title,isvproc=False):
    pltfig,axis = self.cart.create_canvas('')
    iscolor = True
    self.cart.grid_structure(axis,self.mask,self.rank,iscolor=iscolor,isvproc=isvproc)
    self.cart.write_index_in_box(axis,self.mask,self.rank,self.idx_1d,isup=False,iscolor=iscolor,isvproc=isvproc)
    pltfig.savefig(self.figs_dir+'/'+title)
    del axis


  def show(self):
    self.cart.show()


