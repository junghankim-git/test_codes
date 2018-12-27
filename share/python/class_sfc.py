import os
import sys
import numpy as np
sys.path.append('/home/jhkim/work/share/python')
from control_mesh import *
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

moveto = [[-1,0],[1,0],[0,-1],[0,1]]

class SFC:
   def __init__(self, n):
      self.n            = n
      self.n_base_fac   = 3
      self.base_facs    = [2, 3, 5]
      self.factorable   = False
      self.factors      = []
      self.nfactors     = 0
      self.nprim        = 4
      self.major        = [[[] for iprim in range(self.nprim)]  \
                               for ifac in range(self.n_base_fac)]
      self.joiner       = [[[] for iprim in range(self.nprim)]  \
                               for ifac in range(self.n_base_fac)]
      self.ndir         = 4
      self.move         = [[-1,0], [1,0], [0,-1], [0,1]]
      self.result_joiner = [[-1 for i in range(n)] for j in range(n)]

      self.do_factorize(n)
      self.primitive_info()

      n2 = 0
      n3 = 0
      n5 = 0
      for i in self.factors:
         if i == 2: n2=n2+1
         if i == 3: n3=n3+1
         if i == 5: n5=n5+1
      self.name = str(n2)+'hilvert x '+str(n3)+'peano x '+str(n5)+'cinco'



   def do_factorize(self, num):
      # Look for factors of 2 & 3 & 5
      for ifac in self.base_facs:
         tmp = num
         found = True
         while found:
            found = False
            tmpd  = tmp/ifac
            if tmpd*ifac == tmp:
               (self.factors).append(ifac)
               self.nfactors = self.nfactors+1
               tmp   = tmpd
               found = True
      if self.nfactors > 0: self.factorable = True



   def primitive_info(self):
      n_base_fac = self.n_base_fac
      nprim = self.nprim
      ndir  = self.ndir

      il, ir, id, iu, tmp1, tmp2, tmp3, tmp4 = get_direction_index()
      xaxis, yaxis, xyaxis, mxyaxis, orgaxis = get_axis_index()

      self.iprim    = [ir, il, id, iu]
      # [axis=0: x-axis], [axis=1: y-axis], [axis=2: x.y-axis],
      # [axis=3: (-x).y or x.(-y)-axis], [axis=4: (0,0)]
      axis = [orgaxis, xyaxis, mxyaxis]

      # primitive direction
      self.major[0][0]  = [[il,iu],[il,id]]
      self.major[1][0]  = [[il,ir,il],[il,iu,iu],[il,id,id]]
      self.major[2][0]  = [[il,ir,il,ir,il],[il,iu,iu,iu,iu],  \
                           [il,id,id,id,id],[il,iu,iu,ir,il],[il,id,id,ir,il]]
      self.joiner[0][0] = [[iu,-1],[il,id]]
      self.joiner[1][0] = [[iu,ir,-1],[il,iu,il],[il,id,id]]
      self.joiner[2][0] = [[iu,ir,iu,ir,-1],[il,iu,il,iu,il],  \
                           [il,id,ir,id,id],[iu,il,iu,ir,il],[il,id,id,id,il]]
      for ifac in range(n_base_fac):
         for i in range(nprim-1):
            self.major[ifac][i+1]  = symmetry_dir_mesh(self.base_facs[ifac],axis[i], \
                                                                 self.major[ifac][0])
            self.joiner[ifac][i+1] = symmetry_dir_mesh(self.base_facs[ifac],axis[i], \
                                                                self.joiner[ifac][0])



   def get_first_end(self, n, imajor):
      if imajor == 0:
         ifirst, jfirst = n-1, n-1
         iend,  jend  = 0, n-1
      elif imajor == 1:
         ifirst, jfirst = 0, 0
         iend,  jend  = n-1, 0
      elif imajor == 2:
         ifirst, jfirst = n-1, n-1
         iend,  jend  = n-1, 0
      elif imajor == 3:
         ifirst, jfirst = 0, 0
         iend,  jend  = 0, n-1
      return ifirst, jfirst, iend, jend



   def generate_sfc(self):
      if self.nfactors < 1: 
         print 'n was not factorable...'

      #imajor  = 0
      #ijoiner = 0
      imajor  = 1
      ijoiner = 1
      #imajor  = 1
      #ijoiner = 0
      #imajor  = 2
      #ijoiner = 3
      #imajor  = 3
      #ijoiner = 3

      idx     = 1
      ilevel  = 0
      istart, jstart = 0, 0

      nlevels = self.nfactors+1
      size    = [1]
      for i in range(self.nfactors):
          size.append(self.factors[self.nfactors-i-1])

      nsize   = [1 for i in range(nlevels)]
      for ilev in range(1,nlevels):
          nsize[ilev] = size[ilev]*nsize[ilev-1]

      majors  = [[] for ilev in range(nlevels)]
      joiners = [[] for ilev in range(nlevels)]
      majors[0]  = [[imajor]]
      joiners[0] = [[ijoiner]]
      for ilev in range(1,nlevels):
          majors[ilev]  = [[0 for j in range(nsize[ilev])] for i in range(nsize[ilev])]
          joiners[ilev] = [[0 for j in range(nsize[ilev])] for i in range(nsize[ilev])]

      for ilev in range(1,nlevels):
          for i in range(nsize[ilev-1]):
              for j in range(nsize[ilev-1]):
                 imajor  = majors[ilev-1][i][j]
                 ijoiner = joiners[ilev-1][i][j]
                 self.get_major_joiner(ilev,size[ilev],i,j,imajor,ijoiner, \
                                                  majors[ilev],joiners[ilev])

      mesh = self.indexing_sfc(self.n,self.factors[-1],imajor,joiners[nlevels-1])
      self.result_joiner = joiners[nlevels-1]
      return mesh



   def get_joiner(self):
      return self.result_joiner



   def indexing_sfc(self, n, nsfc, imajor, joiners):
      ifirst, jfirst, iend, jend = self.get_first_end(nsfc,imajor)
      mesh = [[0 for j in range(self.n)] for i in range(self.n)]
      idx = 1
      ip, jp = ifirst, jfirst
      for i in range(n*n):
         mesh[ip][jp] = idx
         dir = joiners[ip][jp]
         iinc, jinc = self.move[dir]
         ip  = ip + iinc
         jp  = jp + jinc
         idx = idx + 1

      return mesh



   def get_major_joiner(self, nlev, nsfc, ipos, jpos, imajor, ijoiner, major, joiner):

      ibase, jbase, ilast, jlast = self.get_first_end(nsfc,imajor)
 
      if nsfc == 2:
         ifac = 0
      elif nsfc == 3:
         ifac = 1
      elif nsfc == 5:
         ifac = 2

      for i in range(nsfc):
         for j in range(nsfc):
            ii = nsfc*ipos+i
            jj = nsfc*jpos+j
            major[ii][jj]  = self.major[ifac][imajor][i][j]
            if i == ilast and j == jlast:
               joiner[ii][jj] = ijoiner
            else:
               joiner[ii][jj] = self.joiner[ifac][imajor][i][j]

      return major, joiner



   def get_sfc_line(self):

      n      = self.n
      joiner = self.result_joiner
      dl     = 0.45
      start  = [0,0]
      verts  = []
      codes  = []
      i = start[0]
      j = start[1]
      verts.append([i+dl,j+dl])
      codes.append(Path.MOVETO)
      dir = [0,0]
      for ie in range((n*n)-1):
         dir = moveto[joiner[i][j]]
         i   = i + dir[0]
         j   = j + dir[1]
         verts.append([i+dl,j+dl])
         codes.append(Path.LINETO)
      path  = Path(verts,codes)
      patch = patches.PathPatch(path,lw=1.5,edgecolor='k',facecolor='none',alpha=0.5)
      return patch
      #axis.add_patch(patch)
      #del verts,codes




