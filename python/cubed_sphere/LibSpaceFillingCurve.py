import os
import sys
from numpy import *
sys.path.append('/home/jhkim/Study/Library/Shared/Python')
from LibControlMesh import *
class SFC:
   def __init__(self, n):
      self.n            = n
      self.nbaseFacs    = 3
      self.baseFacs     = [2, 3, 5]
      self.isFactorable = False
      self.factors      = []
      self.nfactors     = 0
      self.nprim        = 4
      self.major        = [[[] for iprim in range(self.nprim)] for ifacSFC in range(self.nbaseFacs)]
      self.joiner       = [[[] for iprim in range(self.nprim)] for ifacSFC in range(self.nbaseFacs)]
      self.ndir         = 4
      self.move         = [[-1,0], [1,0], [0,-1], [0,1]]
      self.ResultJoiner = [[-1 for i in range(n)] for j in range(n)]

      self.DoFactor(n)
      self.GenPrimitiveInfo()

      n2 = 0
      n3 = 0
      n5 = 0
      for i in self.factors:
         if i == 2: n2=n2+1
         if i == 3: n3=n3+1
         if i == 5: n5=n5+1
      self.name = str(n2)+'hilvert x '+str(n3)+'peano x '+str(n5)+'cinco'

   def DoFactor(self, num):
      # Look for factors of 2 & 3 & 5
      for ifac in self.baseFacs:
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
      if self.nfactors > 0: self.isFactorable = True

   def GenPrimitiveInfo(self):
      nbaseFacs = self.nbaseFacs
      nprim = self.nprim
      ndir  = self.ndir


      mymesh = Mesh()
      il, ir, id, iu, tmp1, tmp2, tmp3, tmp4 = mymesh.GetDirectionIndex()
      xaxis, yaxis, xyaxis, mxyaxis, orgaxis = mymesh.GetAxisIndex()

      self.iprim    = [ir, il, id, iu]
      # [axis=0: x-axis], [axis=1: y-axis], [axis=2: x.y-axis], [axis=3: (-x).y or x.(-y)-axis], [axis=4: (0,0)]
      axis = [orgaxis, xyaxis, mxyaxis]


      # - Primitive Direction
      self.major[0][0]  = [[il,iu],[il,id]]
      self.major[1][0]  = [[il,ir,il],[il,iu,iu],[il,id,id]]
      self.major[2][0]  = [[il,ir,il,ir,il],[il,iu,iu,iu,iu],[il,id,id,id,id],[il,iu,iu,ir,il],[il,id,id,ir,il]]
      self.joiner[0][0] = [[iu,-1],[il,id]]
      self.joiner[1][0] = [[iu,ir,-1],[il,iu,il],[il,id,id]]
      self.joiner[2][0] = [[iu,ir,iu,ir,-1],[il,iu,il,iu,il],[il,id,ir,id,id],[iu,il,iu,ir,il],[il,id,id,id,il]]
      for ifacSFC in range(nbaseFacs):
         for i in range(nprim-1):
            self.major[ifacSFC][i+1]  = mymesh.SymmetryDirMesh(self.baseFacs[ifacSFC], axis[i], self.major[ifacSFC][0])
            self.joiner[ifacSFC][i+1] = mymesh.SymmetryDirMesh(self.baseFacs[ifacSFC], axis[i], self.joiner[ifacSFC][0])


   def GetFirstEnd(self, n, imajor):
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


   def GenSFC(self):
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


      mymesh = Mesh()
      for ilev in range(1,nlevels):
          for i in range(nsize[ilev-1]):
              for j in range(nsize[ilev-1]):
                 imajor  = majors[ilev-1][i][j]
                 ijoiner = joiners[ilev-1][i][j]
                 self.GetMajorJoiner(ilev, size[ilev], i, j, imajor, ijoiner, majors[ilev], joiners[ilev])

      #mymesh.PrintMesh2D(nsize[nlevels-1], nsize[nlevels-1], majors[nlevels-1], 'major', True)
      #mymesh.PrintMesh2D(nsize[nlevels-1], nsize[nlevels-1], joiners[nlevels-1], 'joiner', True)

      mesh = self.IndexingSFC(self.n, self.factors[-1], imajor, joiners[nlevels-1])
      self.ResultJoiner = joiners[nlevels-1]
      #mymesh.PrintMesh2D(nsize[nlevels-1], nsize[nlevels-1], self.ResultJoiner, 'result joiner', True)
      #njoiner = mymesh.SymmetryDirMesh(nsize[nlevels-1], nsize[nlevels-1], 0, joiners[nlevels-1])
      #mymesh.PrintMesh2D(nsize[nlevels-1], nsize[nlevels-1], njoiner, 'result joiner', True)
      return mesh


   def GetJoiner(self):
      return self.ResultJoiner


   def IndexingSFC(self, n, nsfc, imajor, joiners):
      ifirst, jfirst, iend, jend = self.GetFirstEnd(nsfc, imajor)
      mesh    = [[0 for j in range(self.n)] for i in range(self.n)]

      idx = 1

      ip, jp = ifirst, jfirst

      for i in range(n*n):
         mesh[ip][jp] = idx
         dir = joiners[ip][jp]
         iinc, jinc = self.move[dir]
         ip    = ip + iinc
         jp    = jp + jinc
         idx = idx + 1

      return mesh

   def GetMajorJoiner(self, nlev, nsfc, ipos, jpos, imajor, ijoiner, major, joiner):

      ibase, jbase, ilast, jlast = self.GetFirstEnd(nsfc, imajor)
 
      if nsfc == 2:
         ifacSFC = 0
      elif nsfc == 3:
         ifacSFC = 1
      elif nsfc == 5:
         ifacSFC = 2

      for i in range(nsfc):
         for j in range(nsfc):
            ii = nsfc*ipos+i
            jj = nsfc*jpos+j
            major[ii][jj]  = self.major[ifacSFC][imajor][i][j]
            if i == ilast and j == jlast:
               joiner[ii][jj] = ijoiner
            else:
               joiner[ii][jj] = self.joiner[ifacSFC][imajor][i][j]

      return major, joiner
