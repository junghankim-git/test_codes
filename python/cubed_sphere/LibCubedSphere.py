import os
import sys
from numpy import *
sys.path.append('/home/jhkim/Study/Library/Shared/Python')
from LibSpaceFillingCurve import *
from LibControlMesh import *

mymesh = Mesh()
ndir = 8
ww, ee, ss, nn, ws, es, wn, en = mymesh.GetDirectionIndex()



class Neighborhood:

   def __init__(self):
      self.cidx  = [-1,-1] # cartesian index
      self.gidx  = -1      # global index
      self.iface = -1      # face number
      self.iproc = -1      # processor number


class Element:
   
   def __init__(self):
      self.nbr    = [ Neighborhood() for i in range(ndir) ]
      self.iface  = -1
      self.cidx   = [-1,-1]
      self.gidx   = -1
      self.isfc   = -1
      self.joiner = -1
      self.iproc  = -1
      self.lidx   = -1
      self.nlidx  = -1
      self.isOut  = False
    
      self.ab     = [[0.0, 0.0] for i in range(4)]
      self.lonlat = [[0.0, 0.0] for i in range(4)]

class CubedSphere:
   def __init__(self, np, ne, nprocs):
      self.nface  = 6
      self.np     = np
      self.ne     = ne
      self.nelem  = self.nface*ne*ne
      self.nprocs = nprocs
      self.ncols  = 3*ne
      self.nrows  = 4*ne
      self.nelemd = [0 for i in range(nprocs)]
      self.nelemo = [0 for i in range(nprocs)]
      if np  == 2:
         self.qp = array([-1.0, 1.0])
      elif np  == 3:
         self.qp = array([-1.0, 0.0, 1.0])
      elif np  == 4:
         self.qp = array([-1.0, -0.5, 0.5, 1.0])
      self.elem  = [[[Element() for k in range(ne)] for j in range(ne)] for i in range(self.nface)]

   def Ini(self):
      print self.elem[0][0][0].__dict__.keys()
      self.FaceNumber()
      self.Cartesian()
      self.GlobalIndex()
      self.Coordinate()
      self.SFCIndex()
      self.DomainDecomposition()
      self.LocalIndex()
      self.Neighborhood()
      self.NewLocalIndex()


   def FaceNumber(self):
      print ' # START: Get face number         ...'
      ne = self.ne
      for iface in range(self.nface):
         for i in range(ne):
            for j in range(ne):
               self.elem[iface][i][j].iface = iface+1
      print ' # END  : Get face number         ...\n'


   def Cartesian(self):
      print ' # START: Get cartesian number    ...'
      ne = self.ne
      for iface in range(self.nface):
         for i in range(ne):
            for j in range(ne):
               self.elem[iface][i][j].cidx = [i+1, j+1]
      print ' # END  : Get cartesian number    ...\n'


   def GlobalIndex(self):
      print ' # START: Get Global index        ...'
      ne = self.ne
      for iface in range(self.nface):
         for i in range(ne):
            for j in range(ne):
               self.elem[iface][i][j].gidx = iface*ne*ne + ne*j + i + 1
      print ' # END  : Get Global index        ...\n'



   # Coordinate
   def AlphaBeta2xy(self, alpha_in, beta_in):
      res_x, res_y = tan(alpha_in), tan(beta_in)
      return res_x, res_y


   def xy2LonLat(self, iface_in, x_in, y_in):
      #print iface_in, x_in, y_in
      if iface_in >= 1 and iface_in <= 4:
         res_lon = arctan(x_in)+pi/2.0*(iface_in-1)
         res_lat = arctan(y_in*cos(res_lon-pi/2.0*(iface_in-1)))
      elif iface_in == 5:
         if x_in == 0.0 and y_in > 0.0:
            res_lon = 0.0
            res_lat = -pi/2.0+arctan(y_in)
         elif x_in == 0.0 and y_in < 0.0:
            res_lon = pi
            res_lat = -pi/2.0-arctan(y_in)
         elif x_in > 0.0 and y_in == 0.0:
            res_lon = pi/2.0
            res_lat = -pi/2.0+arctan(x_in)
         elif x_in < 0.0 and y_in == 0.0:
            res_lon = 3.0*pi/2.0
            res_lat = -pi/2.0-arctan(x_in)
         elif x_in == 0.0 and y_in == 0.0:
            res_lon = 0.0
            res_lat = -pi/2.0
         elif x_in != 0.0 and y_in > 0.0:
            res_lon = arctan(x_in/y_in)
            res_lat = arctan(-sin(res_lon)/x_in)
         elif x_in != 0.0 and y_in < 0.0:
            res_lon = pi+arctan(x_in/y_in)
            res_lat = arctan(-sin(res_lon)/x_in)
      elif iface_in == 6:
         if x_in == 0.0 and y_in > 0.0:
            res_lon = pi
            res_lat = pi/2.0-arctan(y_in)
         elif x_in == 0.0 and y_in < 0.0:
            res_lon = 0.0
            res_lat = pi/2.0+arctan(y_in)
         elif x_in > 0.0 and y_in == 0.0:
            res_lon = pi/2.0
            res_lat = pi/2.0-arctan(x_in)
         elif x_in < 0.0 and y_in == 0.0:
            res_lon = 3.0*pi/2.0
            res_lat = pi/2.0+arctan(x_in)
         elif x_in == 0.0 and y_in == 0.0:
            res_lon = 0.0
            res_lat = pi/2.0
         elif x_in != 0.0 and y_in > 0.0:
            res_lon = pi-arctan(x_in/y_in)
            res_lat = arctan(sin(res_lon)/x_in)
         elif x_in != 0.0 and y_in < 0.0:
            res_lon = arctan(-x_in/y_in)
            res_lat = arctan(sin(res_lon)/x_in)
      else:
         print 'check face number...', iface_in
         quit()
      return res_lon, res_lat


   def AlphaBeta2LonLat(self, iface_in, alpha_in, beta_in):
      int_x, int_y     = self.AlphaBeta2xy(alpha_in, beta_in)
      res_lon, res_lat = self.xy2LonLat(iface_in, int_x, int_y)
      return res_lon, res_lat


   def Coordinate(self):
      print ' # START: Get Coordinate          ...'
      ne = self.ne
      #da, db = 90.0/float(ne), 90.0/float(ne)
      da, db = pi/2.0/float(ne), pi/2.0/float(ne)
      for iface in range(self.nface):
         for i in range(ne):
            for j in range(ne):
               a1 = -pi/4.0+da*float(i)
               a2 = -pi/4.0+da*float(i+1)
               b1 = -pi/4.0+da*float(j)
               b2 = -pi/4.0+da*float(j+1)
               self.elem[iface][i][j].ab  = [[a1,b1],[a2,b1],[a2,b2],[a1,b2]]
               lon11, lat11 = self.AlphaBeta2LonLat(iface+1, a1, b1)
               lon21, lat21 = self.AlphaBeta2LonLat(iface+1, a2, b1)
               lon22, lat22 = self.AlphaBeta2LonLat(iface+1, a2, b2)
               lon12, lat12 = self.AlphaBeta2LonLat(iface+1, a1, b2)
               self.elem[iface][i][j].lonlat = [[lon11,lat11],[lon21,lat21],[lon22,lat22],[lon12,lat12]]
               #print self.elem[iface][i][j].ab 
      print ' # END  : Get Coordinate          ...\n'


   def SFCIndex(self):
      print ' # START: Get SFC    index        ...'
      ne = self.ne
      ## space filling curve
      b_sfc    = array([[0 for i in range(ne)] for j in range(ne)])
      isfc     = SFC(self.ne)
      b_sfc    = isfc.GenSFC()
      #- joiner
      b_joiner = isfc.GetJoiner()
      xjoiner  = mymesh.SymmetryDirMesh(self.ne, 0, b_joiner)
      yjoiner  = mymesh.SymmetryDirMesh(self.ne, 1, b_joiner)
      xyjoiner = mymesh.SymmetryDirMesh(self.ne, 2, b_joiner)
      ojoiner  = mymesh.SymmetryDirMesh(self.ne, 4, b_joiner)
      xyxjoiner= mymesh.SymmetryDirMesh(self.ne, 0, xyjoiner)
      for iface in range(self.nface):
         for i in range(ne):
            for j in range(ne):
               if iface == 0:
                  self.elem[iface][i][j].isfc   = 0       + b_sfc[i][(ne-1)-j]
                  self.elem[iface][i][j].joiner = xjoiner[i][j]
               elif iface == 1:
                  self.elem[iface][i][j].isfc   =  ne*ne  + b_sfc[i][(ne-1)-j]
                  self.elem[iface][i][j].joiner = xjoiner[i][j]
               elif iface == 2:
                  self.elem[iface][i][j].isfc   = 5*ne*ne + b_sfc[i][j]
                  self.elem[iface][i][j].joiner = b_joiner[i][j]
               elif iface == 3:
                  self.elem[iface][i][j].isfc   = 3*ne*ne + b_sfc[(ne-1)-j][i]
                  self.elem[iface][i][j].joiner = xyxjoiner[i][j]
               elif iface == 4:
                  self.elem[iface][i][j].isfc   = 4*ne*ne + b_sfc[i][j]
                  self.elem[iface][i][j].joiner = b_joiner[i][j]
               elif iface == 5:
                  self.elem[iface][i][j].isfc   = 2*ne*ne + b_sfc[(ne-1)-i][(ne-1)-j]
                  self.elem[iface][i][j].joiner = ojoiner[i][j]
      print ' # END  : Get SFC    index        ...\n'

   def DomainDecomposition(self):
      print ' # START: Get Process number      ...'
      nelemd = self.nelem/self.nprocs
      extra  = self.nelem%self.nprocs
      for iproc in range(self.nprocs):
         if iproc < extra:
            self.nelemd[iproc] = nelemd+1
         else:
            self.nelemd[iproc] = nelemd
      np1    = extra*(nelemd+1)
      ne = self.ne
      for iface in range(self.nface):
         for i in range(ne):
            for j in range(ne):
               isfc = self.elem[iface][i][j].isfc - 1
               if isfc <= np1:
                  self.elem[iface][i][j].iproc = isfc/(nelemd+1)+1
               else:
                  self.elem[iface][i][j].iproc = (isfc-np1)/nelemd+1+extra
      print ' # END  : Get Process number      ...\n'


   def LocalIndex(self):
      print ' # START: Get Local index         ...'
      ilocal = array([0 for i in range(self.nprocs)])
      ne = self.ne
      for iface in range(self.nface):
         for j in range(ne):
            for i in range(ne):
               ic = self.elem[iface][i][j].iproc-1
               ilocal[ic] = ilocal[ic] + 1
               self.elem[iface][i][j].lidx = ilocal[ic]
      print ' # END  : Get Local index         ...\n'


   def Neighborhood(self):
      print ' # START: Get Neighborhood Index  ...'
      # west, east, south, north
      ne = self.ne
      for iface in range(self.nface):
        for i in range(ne):
          for j in range(ne):
            left_f  = iface
            left_i  = i-1
            left_j  = j
            right_f = iface
            right_i = i+1
            right_j = j
            up_f    = iface
            up_i    = i
            up_j    = j+1
            down_f  = iface
            down_i  = i
            down_j  = j-1
            if iface == 0:
               if i == 0:
                  left_f  = 3
                  left_i  = ne-1
               elif i == ne-1:
                  right_f = 1
                  right_i = 0
               if j == 0:
                  down_f  = 4
                  down_j  = ne-1
               elif j == ne-1:
                  up_f    = 5
                  up_j    = 0
            elif iface == 1:
               if i == 0:
                  left_f  = 0
                  left_i  = ne-1
               elif i == ne-1:
                  right_f = 2
                  right_i = 0
               if j == 0:
                  down_f  = 4
                  down_i  = ne-1
                  down_j  = (ne-1)-i
               elif j == ne-1:
                  up_f    = 5
                  up_i    = ne-1
                  up_j    = i
            elif iface == 2:
               if i == 0:
                  left_f  = 1
                  left_i  = ne-1
               elif i == ne-1:
                  right_f = 3
                  right_i = 0
               if j == 0:
                  down_f  = 4
                  down_i  = (ne-1)-i
                  down_j  = 0
               elif j == ne-1:
                  up_f    = 5
                  up_i    = (ne-1)-i
                  up_j    = ne-1
            elif iface == 3:
               if i == 0:
                  left_f  = 2
                  left_i  = ne-1
               elif i == ne-1:
                  right_f = 0
                  right_i = 0
               if j == 0:
                  down_f  = 4
                  down_i  = 0
                  down_j  = i
               elif j == ne-1:
                  up_f    = 5
                  up_i    = 0
                  up_j    = (ne-1)-i
            elif iface == 4:
               if i == 0:
                  left_f  = 3
                  left_i  = j
                  left_j  = 0
               elif i == ne-1:
                  right_f = 1
                  right_i = (ne-1)-j
                  right_j = 0
               if j == 0:
                  down_f  = 2
                  down_j  = ne-1
               elif j == ne-1:
                  up_f    = 0
                  up_j    = 0
            elif iface == 5:
               if i == 0:
                  left_f  = 3
                  left_i  = (ne-1)-j
                  left_j  = ne-1
               elif i == ne-1:
                  right_f = 1
                  right_i = j
                  right_j = ne-1
               if j == 0:
                  down_f  = 0
                  down_j  = ne-1
               elif j == ne-1:
                  up_f    = 4
                  up_j    = 0
            # left
            self.elem[iface][i][j].nbr[ww].cidx  = [left_i+1, left_j+1]
            self.elem[iface][i][j].nbr[ww].gidx  = self.elem[left_f][left_i][left_j].gidx
            self.elem[iface][i][j].nbr[ww].iface = self.elem[left_f][left_i][left_j].iface
            self.elem[iface][i][j].nbr[ww].iproc = self.elem[left_f][left_i][left_j].iproc
            # right
            self.elem[iface][i][j].nbr[ee].cidx  = [right_i+1, right_j+1]
            self.elem[iface][i][j].nbr[ee].gidx  = self.elem[right_f][right_i][right_j].gidx
            self.elem[iface][i][j].nbr[ee].iface = self.elem[right_f][right_i][right_j].iface
            self.elem[iface][i][j].nbr[ee].iproc = self.elem[right_f][right_i][right_j].iproc
            # down
            self.elem[iface][i][j].nbr[ss].cidx  = [down_i+1, down_j+1]
            self.elem[iface][i][j].nbr[ss].gidx  = self.elem[down_f][down_i][down_j].gidx
            self.elem[iface][i][j].nbr[ss].iface = self.elem[down_f][down_i][down_j].iface
            self.elem[iface][i][j].nbr[ss].iproc = self.elem[down_f][down_i][down_j].iproc
            # up
            #print iface, i, j, ':', up_i, up_j
            self.elem[iface][i][j].nbr[nn].cidx  = [up_i+1, up_j+1]
            self.elem[iface][i][j].nbr[nn].gidx  = self.elem[up_f][up_i][up_j].gidx
            self.elem[iface][i][j].nbr[nn].iface = self.elem[up_f][up_i][up_j].iface
            self.elem[iface][i][j].nbr[nn].iproc = self.elem[up_f][up_i][up_j].iproc

      # west-south, east-south, west-north, east-north
      for iface in range(self.nface):
        for i in range(ne):
          for j in range(ne):
            left_f  = self.elem[iface][i][j].nbr[ww].iface-1
            left_i  = self.elem[iface][i][j].nbr[ww].cidx[0]-1
            left_j  = self.elem[iface][i][j].nbr[ww].cidx[1]-1
            right_f = self.elem[iface][i][j].nbr[ee].iface-1
            right_i = self.elem[iface][i][j].nbr[ee].cidx[0]-1
            right_j = self.elem[iface][i][j].nbr[ee].cidx[1]-1
            down_f  = self.elem[iface][i][j].nbr[ss].iface-1
            down_i  = self.elem[iface][i][j].nbr[ss].cidx[0]-1
            down_j  = self.elem[iface][i][j].nbr[ss].cidx[1]-1
            up_f    = self.elem[iface][i][j].nbr[nn].iface-1
            up_i    = self.elem[iface][i][j].nbr[nn].cidx[0]-1
            up_j    = self.elem[iface][i][j].nbr[nn].cidx[1]-1

            # ws
            leftdown_f = self.elem[left_f][left_i][left_j].nbr[ss].iface-1
            leftdown_i = self.elem[left_f][left_i][left_j].nbr[ss].cidx[0]-1
            leftdown_j = self.elem[left_f][left_i][left_j].nbr[ss].cidx[1]-1
            # es
            rightdown_f = self.elem[right_f][right_i][right_j].nbr[ss].iface-1
            rightdown_i = self.elem[right_f][right_i][right_j].nbr[ss].cidx[0]-1
            rightdown_j = self.elem[right_f][right_i][right_j].nbr[ss].cidx[1]-1
            # wn
            leftup_f = self.elem[left_f][left_i][left_j].nbr[nn].iface-1
            leftup_i = self.elem[left_f][left_i][left_j].nbr[nn].cidx[0]-1
            leftup_j = self.elem[left_f][left_i][left_j].nbr[nn].cidx[1]-1
            # en
            rightup_f = self.elem[right_f][right_i][right_j].nbr[nn].iface-1
            rightup_i = self.elem[right_f][right_i][right_j].nbr[nn].cidx[0]-1
            rightup_j = self.elem[right_f][right_i][right_j].nbr[nn].cidx[1]-1
            if (iface == 4) or (iface == 5):
               if i == 0:
                  # ws
                  leftdown_f = self.elem[down_f][down_i][down_j].nbr[ww].iface-1
                  leftdown_i = self.elem[down_f][down_i][down_j].nbr[ww].cidx[0]-1
                  leftdown_j = self.elem[down_f][down_i][down_j].nbr[ww].cidx[1]-1
                  # wn
                  leftup_f = self.elem[up_f][up_i][up_j].nbr[ww].iface-1
                  leftup_i = self.elem[up_f][up_i][up_j].nbr[ww].cidx[0]-1
                  leftup_j = self.elem[up_f][up_i][up_j].nbr[ww].cidx[1]-1
               if i == ne-1:
                  # es
                  rightdown_f = self.elem[down_f][down_i][down_j].nbr[ee].iface-1
                  rightdown_i = self.elem[down_f][down_i][down_j].nbr[ee].cidx[0]-1
                  rightdown_j = self.elem[down_f][down_i][down_j].nbr[ee].cidx[1]-1
                  # en
                  rightup_f = self.elem[up_f][up_i][up_j].nbr[ee].iface-1
                  rightup_i = self.elem[up_f][up_i][up_j].nbr[ee].cidx[0]-1
                  rightup_j = self.elem[up_f][up_i][up_j].nbr[ee].cidx[1]-1
            #print iface, i, j, ' ws: ', leftdown_f, leftdown_i, leftdown_j
            #print iface, i, j, ' es: ', rightdown_f, rightdown_i, rightdown_j
            #print iface, i, j, ' wn: ', leftup_f, leftup_i, leftup_j
            #print iface, i, j, ' en: ', rightup_f, rightup_i, rightup_j
            # set ws
            #if iface == 1 and i == 0: print leftdown_f, leftdown_i, leftdown_j, self.elem[leftdown_f][leftdown_i][leftdown_j].iproc
            self.elem[iface][i][j].nbr[ws].cidx  = [leftdown_i+1, leftdown_j+1]
            self.elem[iface][i][j].nbr[ws].gidx  = self.elem[leftdown_f][leftdown_i][leftdown_j].gidx
            self.elem[iface][i][j].nbr[ws].iface = self.elem[leftdown_f][leftdown_i][leftdown_j].iface
            self.elem[iface][i][j].nbr[ws].iproc = self.elem[leftdown_f][leftdown_i][leftdown_j].iproc
            # set es
            self.elem[iface][i][j].nbr[es].cidx  = [rightdown_i+1, rightdown_j+1]
            self.elem[iface][i][j].nbr[es].gidx  = self.elem[rightdown_f][rightdown_i][rightdown_j].gidx
            self.elem[iface][i][j].nbr[es].iface = self.elem[rightdown_f][rightdown_i][rightdown_j].iface
            self.elem[iface][i][j].nbr[es].iproc = self.elem[rightdown_f][rightdown_i][rightdown_j].iproc
            # set wn
            self.elem[iface][i][j].nbr[wn].cidx  = [leftup_i+1, leftup_j+1]
            self.elem[iface][i][j].nbr[wn].gidx  = self.elem[leftup_f][leftup_i][leftup_j].gidx
            self.elem[iface][i][j].nbr[wn].iface = self.elem[leftup_f][leftup_i][leftup_j].iface
            self.elem[iface][i][j].nbr[wn].iproc = self.elem[leftup_f][leftup_i][leftup_j].iproc
            # set en
            self.elem[iface][i][j].nbr[en].cidx  = [rightup_i+1, rightup_j+1]
            self.elem[iface][i][j].nbr[en].gidx  = self.elem[rightup_f][rightup_i][rightup_j].gidx
            self.elem[iface][i][j].nbr[en].iface = self.elem[rightup_f][rightup_i][rightup_j].iface
            self.elem[iface][i][j].nbr[en].iproc = self.elem[rightup_f][rightup_i][rightup_j].iproc

            if i == 0  and j == 0:
               self.elem[iface][i][j].nbr[ws].cidx  = [-1, -1]
               self.elem[iface][i][j].nbr[ws].gidx  = -1
               self.elem[iface][i][j].nbr[ws].iface = -1
               self.elem[iface][i][j].nbr[ws].iproc = -1
            if i == ne and j == 0:
               self.elem[iface][i][j].nbr[es].cidx  = [-1, -1]
               self.elem[iface][i][j].nbr[es].gidx  = -1
               self.elem[iface][i][j].nbr[es].iface = -1
               self.elem[iface][i][j].nbr[es].iproc = -1
            if i == 0  and j == ne:
               self.elem[iface][i][j].nbr[wn].cidx  = [-1, -1]
               self.elem[iface][i][j].nbr[wn].gidx  = -1
               self.elem[iface][i][j].nbr[wn].iface = -1
               self.elem[iface][i][j].nbr[wn].iproc = -1
            if i == ne and j == ne:
               self.elem[iface][i][j].nbr[en].cidx  = [-1, -1]
               self.elem[iface][i][j].nbr[en].gidx  = -1
               self.elem[iface][i][j].nbr[en].iface = -1
               self.elem[iface][i][j].nbr[en].iproc = -1
      print ' # END  : Get Neighborhood Index  ...\n'


   def NewLocalIndex(self):
      print ' # START: Get New Local index     ...'
      ne = self.ne
      ielem = [0 for i in range(self.nprocs)]
      for iproc in range(self.nprocs):
         for iface in range(self.nface):
            for j in range(ne):
               for i in range(ne):
                  for l in range(ndir):
                     iiproc = self.elem[iface][i][j].nbr[l].iproc
                     if iiproc == -1: continue
                     if (self.elem[iface][i][j].iproc == iproc+1) and (iiproc != iproc+1):
                        #if iproc == 3 and j == 0:
                        #   print '('+str(iface)+','+str(i)+','+str(j)+'): nbr:'+str(l)+'  '+str(iiproc)
                        ielem[iproc] = ielem[iproc] + 1
                        self.elem[iface][i][j].nlidx = ielem[iproc]
                        self.elem[iface][i][j].isOut = True
                        break
         self.nelemo[iproc] = ielem[iproc]
         for iface in range(self.nface):
            for j in range(ne):
               for i in range(ne):
                  if (not self.elem[iface][i][j].isOut) and (self.elem[iface][i][j].iproc == iproc+1):
                     ielem[iproc] = ielem[iproc] + 1
                     self.elem[iface][i][j].nlidx = ielem[iproc]

      #print ' nelemd =', self.nelemd
      #print ' nelemo =', self.nelemo
      print ' # END  : Get New Local index     ...\n'


   def GetInfoOuterElement(self):
      return max(self.nelemd), std(self.nelemd), mean(self.nelemo), std(self.nelemo), max(self.nelemo), min(self.nelemo)



   def CanvasToCubeAxis(self, iface, ic, jc):
      ne = self.ne
      base = [[2*ne, ne], [2*ne, 2*ne], [2*ne, 3*ne], [2*ne, 0], [3*ne, ne], [ne, ne]]
      i = jc - base[iface][1]
      j = base[iface][0] - ic -1
      return i, j

   def CubeToCanvasAxis(self, iface, i, j):
      ne = self.ne
      base = [[2*ne, ne], [2*ne, 2*ne], [2*ne, 3*ne], [2*ne, 0], [3*ne, ne], [ne, ne]]
      ic = base[iface][0] - j - 1
      jc = base[iface][1] + i
      return ic, jc


   def CubeToCanvas(self, ival):
      # for multi subplots - old version
      #           0         1       2         3      4      5        6       7        8        9
      # ival = ['isOut', 'nlidx', 'iface', 'iproc', 'nbr', 'gidx', 'cidx', 'lidx', 'joiner', 'isfc']
      if ival==0:
         canvas = array([[False for j in range(self.nrows)] for i in range(self.ncols)])
      elif ival==4:
         canvas = array([[[-1 for k in range(ndir)] for j in range(self.nrows)] for i in range(self.ncols)])
      elif ival==6:
         canvas = array([[[-1,-1] for j in range(self.nrows)] for i in range(self.ncols)])
      elif ival==7 or ival==8:
         canvas = array([[[-1 for k in range(4)] for j in range(self.ncols)] for i in range(self.nrows)])
      else:
         canvas = array([[-1 for j in range(self.nrows)] for i in range(self.ncols)])

      ne = self.ne
      base = [[2*ne, ne], [2*ne, 2*ne], [2*ne, 3*ne], [2*ne, 0], [3*ne, ne], [ne, ne]]

      for iface in range(self.nface):
         for i in range(self.ne):
            for j in range(self.ne):
               ic = base[iface][0] - j - 1
               jc = base[iface][1] + i
               #print self.elem[iface][i][j].__dict__.values()[ival]
               if ival == 4:
                  for k in range(8):
                     canvas[ic][jc][k] = self.elem[iface][i][j].__dict__.values()[ival][k].iproc
               else:
                  canvas[ic][jc] = self.elem[iface][i][j].__dict__.values()[ival]
      return canvas


   def CubeToPlot(self, ival):
      # for 1 plot
      #==old           0         1       2         3      4      5        6       7        8        9
      #==old ival = ['isOut', 'nlidx', 'iface', 'iproc', 'nbr', 'gidx', 'cidx', 'lidx', 'joiner', 'isfc']
      #           0         1       2         3      4      5        6       7      8        9        10       11
      # ival = ['isOut', 'nlidx', 'iface', 'iproc', 'nbr', 'gidx', 'cidx', 'ab', 'lonlat', 'lidx', 'joiner', 'isfc']
      if ival==0:
         canvas = array([[False for j in range(self.ncols)] for i in range(self.nrows)])
      elif ival==4:
         canvas = array([[[-1 for k in range(ndir)] for j in range(self.ncols)] for i in range(self.nrows)])
      elif ival==6:
         canvas = array([[[-1,-1] for j in range(self.ncols)] for i in range(self.nrows)])
      elif ival==7 or ival==8:
         canvas = array([[[[-1.0,-1.0] for k in range(4)] for j in range(self.ncols)] for i in range(self.nrows)])
      else:
         canvas = array([[-1 for j in range(self.ncols)] for i in range(self.nrows)])

      ne = self.ne
      base = [[ne, ne], [2*ne, ne], [3*ne, ne], [0, ne], [ne, 0], [ne, 2*ne]]

      for iface in range(self.nface):
         for i in range(self.ne):
            for j in range(self.ne):
               ic = base[iface][0] + i
               jc = base[iface][1] + j
               if ival == 4:
                  for k in range(8):
                     canvas[ic][jc][k] = self.elem[iface][i][j].__dict__.values()[ival][k].iproc
               else:
                  #print self.elem[iface][i][j].__dict__.keys()[ival]
                  canvas[ic][jc] = self.elem[iface][i][j].__dict__.values()[ival]
                  #print canvas[ic][jc]
                  #print self.elem[iface][i][j].ab
      return canvas



