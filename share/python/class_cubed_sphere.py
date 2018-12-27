#!/usr/bin/env python
import os
import sys
import time
import numpy
from scipy import *
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_toolkits.basemap import Basemap
SHARE_DIR='/home/jhkim/work/share/python'
sys.path.append(SHARE_DIR)
from class_draw_box import *


moveto = [[-1,0],[1,0],[0,-1],[0,1]]



class neighbor:

   def __init__(self):
      self.ci    = -1      # cartesian index
      self.cj    = -1      # cartesian index
      self.grb   = -1      # global index
      self.iface = -1      # face number
      self.rank  = -1      # processor number




class cubed_sphere:

   def __init__(self, filename,simple=False):
      self.filename = filename
      self.load_file(simple)

      # quadrature point position 
      if self.np==4:
         self.qp = [-1.00,-0.5,0.5,1.00]
      elif self.np==3:
         self.qp = [-1.00,0.0,1.00]
      else:
         print 'check np...'

      self.plot = class_draw_box(self.np,self.qp,self.pnx,self.pny)
      self.generate_flags(simple)
      self.dpi = 300 # default: 300



   def load_file(self,simple):

      infile = nc.Dataset(self.filename, mode='r')
      self.nface  = len(infile.dimensions['nface'])
      self.np     = len(infile.dimensions['np'])
      self.ne     = len(infile.dimensions['ne'])
      self.nprocs = len(infile.dimensions['nprocs'])
      self.pnx    = len(infile.dimensions['pnx'])
      self.pny    = len(infile.dimensions['pny'])
      if not simple:
         self.nnbrs = len(infile.dimensions['nnbrs'])
         self.n_nbr = len(infile.dimensions['n_nbr'])
      nface  = self.nface
      np     = self.np
      ne     = self.ne
      nprocs = self.nprocs
      pnx    = self.pnx
      pny    = self.pny
      if not simple:
         nnbrs = self.nnbrs
         n_nbr = self.n_nbr
      print ' '
      print '* nface     = ', nface
      print '* np        = ', np
      print '* ne        = ', ne
      print '* nprocs    = ', nprocs
      print '* pnx       = ', pnx
      print '* pny       = ', pny
      print ' '
  
      # attributes
      if not simple:
         rotation = infile.rotation
         if rotation==1:
            self.rotation = True
         else:
            self.rotation = False
     
      # view processor
      self.vproc = 1
      
      self.allocate_variables_numpy(simple)
      #self.allocate_variables_list(simple)

      iface_m = infile.variables['iface'][:].tolist()
      ci_m    = infile.variables['ie'][:].tolist()
      cj_m    = infile.variables['je'][:].tolist()
      grb_m   = infile.variables['global'][:].tolist()
      rank_m  = infile.variables['rank'][:].tolist()
      iup_m   = infile.variables['iup'][:].tolist()
      isu_m   = infile.variables['isunique'][:].tolist()
      if not simple:
         isfc_m    = infile.variables['isfc'][:].tolist()
         lcl_m     = infile.variables['local'][:].tolist()
         joiner_m  = infile.variables['joiner'][:].tolist()
         nbrs_m    = infile.variables['nbrs'][:].tolist()
         alpha_m   = infile.variables['alpha'][:].tolist()
         beta_m    = infile.variables['beta'][:].tolist()
         lon_m     = infile.variables['lon'][:].tolist()
         lat_m     = infile.variables['lat'][:].tolist()
      else:
         v1_m    = infile.variables['var1'][:].tolist()
         v2_m    = infile.variables['var2'][:].tolist()
         v3_m    = infile.variables['var3'][:].tolist()
         v4_m    = infile.variables['var4'][:].tolist()
         v5_m    = infile.variables['var5'][:].tolist()
         v6_m    = infile.variables['var6'][:].tolist()

      # fortran index to c(python) index
      for i in range(pnx):
         for j in range(pny):
            self.iface[i][j]  = iface_m[j][i]
            self.ci[i][j]     = ci_m[j][i]
            self.cj[i][j]     = cj_m[j][i]
            self.cij[i][j][0] = self.ci[i][j]
            self.cij[i][j][1] = self.cj[i][j]
            self.grb[i][j]    = grb_m[j][i]
            self.rank[i][j]   = rank_m[j][i]
            for x in range(np):
               for y in range(np):
                  self.iup[i][j][x][y] = iup_m[j][i][y][x]
                  self.isu[i][j][x][y] = isu_m[j][i][y][x]
            if not simple:
               self.isfc[i][j]   = isfc_m[j][i]
               self.lcl[i][j]    = lcl_m[j][i]
               self.joiner[i][j] = joiner_m[j][i]
               for k in range(nnbrs):
                  self.nbrs[i][j][k].ci    = nbrs_m[j][i][k][0]
                  self.nbrs[i][j][k].cj    = nbrs_m[j][i][k][1]
                  self.nbrs[i][j][k].grb   = nbrs_m[j][i][k][2]
                  self.nbrs[i][j][k].iface = nbrs_m[j][i][k][3]
                  self.nbrs[i][j][k].rank  = nbrs_m[j][i][k][4]
               self.alpha[i][j] = alpha_m[j][i]
               self.beta[i][j]  = beta_m[j][i]
               self.lon[i][j]   = lon_m[j][i]
               self.lat[i][j]   = lat_m[j][i]
            else:
               for x in range(np):
                  for y in range(np):
                     self.value[0][i][j][x][y] = v1_m[j][i][y][x]
                     self.value[1][i][j][x][y] = v2_m[j][i][y][x]
                     self.value[2][i][j][x][y] = v3_m[j][i][y][x]
                     self.value[3][i][j][x][y] = v4_m[j][i][y][x]
                     self.value[4][i][j][x][y] = v5_m[j][i][y][x]
                     self.value[5][i][j][x][y] = v6_m[j][i][y][x]
  
      del iface_m,ci_m,cj_m,grb_m,rank_m,iup_m,isu_m
      if not simple:
         del isfc_m,lcl_m,joiner_m,alpha_m,beta_m,lon_m,lat_m
      else:
         del v1_m,v2_m,v3_m,v4_m,v5_m,v6_m

      infile.close()



   def allocate_variables_list(self,simple):

      np     = self.np
      pnx    = self.pnx
      pny    = self.pny
      if not simple:
         nnbrs = self.nnbrs
         n_nbr = self.n_nbr

      self.iface = [[0 for j in range(pny)] for i in range(pnx)]
      self.ci    = [[0 for j in range(pny)] for i in range(pnx)]
      self.cj    = [[0 for j in range(pny)] for i in range(pnx)]
      self.cij   = [[[0 for k in range(2)] for j in range(pny)] for i in range(pnx)]
      self.grb   = [[0 for j in range(pny)] for i in range(pnx)]
      self.rank  = [[0 for j in range(pny)] for i in range(pnx)]
      self.iup   = [[[[-1 for l in range(np)] for k in range(np)] \
                         for j in range(pny)] for i in range(pnx)]
      self.isu   = [[[[-1 for l in range(np)] for k in range(np)] \
                         for j in range(pny)] for i in range(pnx)]
  
      if not simple:
         self.isfc   = [[0 for j in range(pny)] for i in range(pnx)]
         self.lcl    = [[0 for j in range(pny)] for i in range(pnx)]
         self.joiner = [[0 for j in range(pny)] for i in range(pnx)]
         self.nbrs   = [[[neighbor() for k in range(nnbrs)] for j in range(pny)] \
                                                             for i in range(pnx)]
         self.alpha = [[[0.0 for k in range(4)] for j in range(pny)] for i in range(pnx)]
         self.beta  = [[[0.0 for k in range(4)] for j in range(pny)] for i in range(pnx)]
         self.lon   = [[[0.0 for k in range(4)] for j in range(pny)] for i in range(pnx)]
         self.lat   = [[[0.0 for k in range(4)] for j in range(pny)] for i in range(pnx)]
      else:
         self.value = [[[[[0.0 for l in range(np)] for k in range(np)]  \
                      for j in range(pny)] for i in range(pnx)] for o in range(6)]



   def allocate_variables_numpy(self,simple):

      np     = self.np
      pnx    = self.pnx
      pny    = self.pny
      if not simple:
         nnbrs = self.nnbrs
         n_nbr = self.n_nbr

      self.iface = numpy.empty(shape=(pnx,pny)).tolist()
      self.ci    = numpy.empty(shape=(pnx,pny)).tolist()
      self.cj    = numpy.empty(shape=(pnx,pny)).tolist()
      self.cij   = numpy.empty(shape=(pnx,pny,2)).tolist()
      self.grb   = numpy.empty(shape=(pnx,pny)).tolist()
      self.rank  = numpy.empty(shape=(pnx,pny)).tolist()
      self.iup   = numpy.empty(shape=(pnx,pny,np,np)).tolist()
      self.isu   = numpy.empty(shape=(pnx,pny,np,np)).tolist()
  
      if not simple:
         self.isfc   = numpy.empty(shape=(pnx,pny)).tolist()
         self.lcl    = numpy.empty(shape=(pnx,pny)).tolist()
         self.joiner = numpy.empty(shape=(pnx,pny)).tolist()
         #self.nbrs   = numpy.full(shape=(pnx,pny,nnbrs),fill_value=neighbor(), \
         #                                                dtype=object).tolist()
         self.nbrs   = numpy.empty(shape=(pnx,pny,nnbrs),dtype=object).tolist()
         for ie in range(pnx):
            for je in range(pny):
               for k in range(nnbrs):
                  self.nbrs[ie][je][k] = neighbor()
         self.alpha  = numpy.empty(shape=(pnx,pny,4)).tolist()
         self.beta   = numpy.empty(shape=(pnx,pny,4)).tolist()
         self.lon    = numpy.empty(shape=(pnx,pny,4)).tolist()
         self.lat    = numpy.empty(shape=(pnx,pny,4)).tolist()
      else:
         self.value = numpy.empty(shape=(6,pnx,pny,np,np)).tolist()



   def generate_flags(self,simple=False):

      # flags
      pnx = self.pnx
      pny = self.pny
      grbs = [1,2,3,4,5,7,21,22]
      # flag1 (draw box, e.g. face)
      self.flag1 = [[False for j in range(pny)] for i in range(pnx)]
      self.flag2 = [[False for j in range(pny)] for i in range(pnx)]
      self.flag3 = [[False for j in range(pny)] for i in range(pnx)]
      self.flag4 = [[False for j in range(pny)] for i in range(pnx)]
      self.fflag = [[False for j in range(pny)] for i in range(pnx)]
      self.color = [['#ffffff' for j in range(pny)] for i in range(pnx)]
  
      for ie in range(self.pnx):
         for je in range(self.pny):
            if self.iface[ie][je] != -1:
               self.flag1[ie][je] = True
            if grbs == -1 or grbs.count(self.grb[ie][je])!=0:
               self.flag2[ie][je] = True
  
            iproc = self.rank[ie][je]
            if iproc == self.vproc:
               self.flag3[ie][je] = True
            else:                           # draw neighbors
               if not simple:
                  for inbr in arange(self.nnbrs):
                     iproc_nbr = self.nbrs[ie][je][inbr].rank
                     if iproc_nbr == self.vproc:
                        self.flag4[ie][je] = True
      self.color = self.plot.get_box_color(self.rank)



####################################################
###                Functions                     ###
####################################################

   def xy_to_lon(self, x, y):
      if x==0.0 and y==0.0:
         theta = 0.0
      elif x>0.0 and y>=0.0 or x==0.0 and y>0.0:
         theta = arctan(y/x)
      elif x<0.0 and y>0.0:
         theta = pi+arctan(y/x)
      elif x<0.0 and y<=0.0:
         theta = pi+arctan(y/x)
      elif x>= 0.0 and y<0.0:
         theta = 2.0*pi+arctan(y/x)
      else:
         print 'check xy_to_lon...'
      if theta >= 2.0*pi: theta = theta-2.0*pi

      return theta
 
 

   def lonlat_to_cartesian(self, r, lon, lat):
      #if lon<0.0 or lon>2.0*pi:
      #  print 'check the range of lon... in lonlat_to_cartesian',lon
      x = r*cos(lat)*cos(lon)
      y = r*cos(lat)*sin(lon)
      z = r*sin(lat)
      return x, y, z
 

 
   def cartesian_to_lonlat(self, x, y, z):
      r   = sqrt(x*x+y*y+z*z)
      lon = self.xy_to_lon(x, y)
      lat = arcsin(z/r)
      return r, lon, lat
   
   
   def galphabeta_to_xy(self, a, b):
      x, y = tan(a), tan(b)
      return x, y
 
   
   def xy_to_lonlat(self, iface, x, y, rotation=False):
      if iface >= 1 and iface <= 4:
         lon = arctan(x)+pi/2.0*(iface-1)
         lat = arctan(y*cos(lon-pi/2.0*(iface-1)))
      elif iface == 5:
         if x == 0.0 and y > 0.0:
            lon = 0.0
            lat = -pi/2.0+arctan(y)
         elif x == 0.0 and y < 0.0:
            lon = pi
            lat = -pi/2.0-arctan(y)
         elif x > 0.0 and y == 0.0:
            lon = pi/2.0
            lat = -pi/2.0+arctan(x)
         elif x < 0.0 and y == 0.0:
            lon = 3.0*pi/2.0
            lat = -pi/2.0-arctan(x)
         elif x == 0.0 and y == 0.0:
            lon = 0.0
            lat = -pi/2.0
         elif x != 0.0 and y > 0.0:
            lon = arctan(x/y)
            lat = arctan(-sin(lon)/x)
         elif x != 0.0 and y < 0.0:
            lon = pi+arctan(x/y)
            lat = arctan(-sin(lon)/x)
      elif iface == 6:
         if x == 0.0 and y > 0.0:
            lon = pi
            lat = pi/2.0-arctan(y)
         elif x == 0.0 and y < 0.0:
            lon = 0.0
            lat = pi/2.0+arctan(y)
         elif x > 0.0 and y == 0.0:
            lon = pi/2.0
            lat = pi/2.0-arctan(x)
         elif x < 0.0 and y == 0.0:
            lon = 3.0*pi/2.0
            lat = pi/2.0+arctan(x)
         elif x == 0.0 and y == 0.0:
            lon = 0.0
            lat = pi/2.0
         elif x != 0.0 and y > 0.0:
            lon = pi-arctan(x/y)
            lat = arctan(sin(lon)/x)
         elif x != 0.0 and y < 0.0:
            lon = arctan(-x/y)
            lat = arctan(sin(lon)/x)
      else:
         print 'check face number...', iface
         quit()
 
      if rotation:
         lon0 = 127.0*pi/180.0
         lat0 =  38.0*pi/180.0
         xx, yy, zz = self.lonlat_to_cartesian(1.0,lon,lat)
         rxx = cos(lat0)*cos(lon0)*xx-sin(lon0)*yy-sin(lat0)*cos(lon0)*zz
         ryy = cos(lat0)*sin(lon0)*xx+cos(lon0)*yy-sin(lat0)*sin(lon0)*zz
         rzz = sin(lat0)*xx+cos(lat0)*zz
         rrr, lon, lat = self.cartesian_to_lonlat(rxx,ryy,rzz)
 
      return lon, lat
   
   
   def alphabeta_to_lonlat(self, iface, a, b):
      x, y = self.galphabeta_to_xy(a,b)
      lon, lat = self.xy_to_lonlat(iface,x,y,self.rotation)
      return lon, lat
   
   
   
   def two_alphabeta_to_lonlat_line(self, iface, a1, b1, a2, b2, npts):

      lons = [0.0 for i in range(npts)]
      lats = [0.0 for i in range(npts)]
   
      if (a1 == a2) and (b1 != b2):
         db = abs(b2-b1)/(npts-1)
         for ip in range(npts):
            b = min(b1,b2)+ip*db
            lons[ip],lats[ip] = self.alphabeta_to_lonlat(iface,a1,b)
      elif (a1 != a2) and (b1 == b2):
         da = abs(a2-a1)/(npts-1)
         for ip in range(npts):
            a = min(a1,a2)+ip*da
            lons[ip],lats[ip] = self.alphabeta_to_lonlat(iface,a,b1)
      else:
         print 'Error : two_alphabeta_to_lonlat_line...'
         quit()
      return lons, lats
   
   

   def get_alphabeta_square(self, iface, npts, a4, b4):
      lons = [[0.0 for j in range(npts)] for i in range(4)]
      lats = [[0.0 for j in range(npts)] for i in range(4)]
      for k in range(4):
         if k == 3:
            lons[k],lats[k] = self.two_alphabeta_to_lonlat_line(iface,a4[k],b4[k], \
                                                              a4[0],b4[0],npts)
         else:
            lons[k],lats[k] = self.two_alphabeta_to_lonlat_line(iface,a4[k],b4[k], \
                                                              a4[k+1],b4[k+1],npts)
      return lons, lats
 
 
 
 
   def show(self):
      self.plot.show()
   
   
   
   
   def set_view(self, axis, ne, np, nprocs, vproc, optp, optvp):
      # basic position
      base    = [ 3.0, 2.0]
      tpos    = [-1.0,-1.0]
      # size
      tsize   = 25.0 #100./ne
      # postion
      tpos[0] = ne*(base[0]-0.1)
      tpos[1] = ne*(base[1]+0.75)
      # text
      axis.text(tpos[0],tpos[1],'np = '+str(np),color='k',size=tsize, \
                        va='center',style='italic')#,fontname='sans-serif'
      tpos[1] = ne*(base[1]+0.60)
      axis.text(tpos[0],tpos[1],'ne = '+str(ne),color='k',size=tsize, \
                        va='center',style='italic')#,fontname='sans-serif'
      if optp:
         tpos[1] = ne*(base[1]+0.45)
         axis.text(tpos[0],tpos[1],'nprocs = '+str(nprocs),color='k',size=tsize, \
                        va='center',style='italic')#,fontname='sans-serif'
      if optvp:
         tpos[1] = ne*(base[1]+0.3)
         axis.text(tpos[0],tpos[1],'iproc = '+str(vproc),color='k',size=tsize, \
                        va='center',style='italic')#,fontname='sans-serif'
 

 
   def draw_face_number(self, axis):
      dx    = 2
      nface = self.nface
      ne    = self.ne
      base  = [[1.0,1.0],[2.0,1.0],[3.0,1.0],[0.0,1.0],[1.0,0.0],[1.0,2.0]]
      tpos  = [[-1,-1] for i in range(nface)]
      for iface in range(nface):
         #tpos[iface][0] = ne*(base[iface][0]+0.5)
         #tpos[iface][1] = ne*(base[iface][1]+0.5)
         tpos[iface][0] = ne*(base[iface][0]+0.45)
         tpos[iface][1] = ne*(base[iface][1]+0.45)
      for iface in range(nface):
         #axis.text(tpos[iface][0],tpos[iface][1],str(iface+1),color='k',size=150.0,alpha=0.3,weight='semibold',fontname='sans-serif',style='italic',ha='center',va='center')
         axis.text(tpos[iface][0],tpos[iface][1],str(iface+1),color='b',size=130.0, \
                   alpha=0.3,fontname='sans-serif',style='italic',ha='center',va='center')
 
   
   
   def draw_sfc_line(self, axis):
      ne    = self.ne
      nface = self.nface
     
      dl    = 0.45
      start = [[ne,2*ne-1],    [2*ne,2*ne-1],  [3*ne,ne],  [0,2*ne-1], \
               [ne,0],         [2*ne-1,3*ne-1]]
      end   = [[2*ne-1,2*ne-1],[3*ne-1,2*ne-1],[4*ne-1,ne],[0,ne],     \
               [2*ne-1,0],     [ne,3*ne-1]]
      for iface in range(nface):
         verts = []
         codes = []
         i = start[iface][0]
         j = start[iface][1]
         verts.append([i+dl,j+dl])
         codes.append(Path.MOVETO)
         dir = [0,0]
         for ie in range((ne*ne)-1):
            dir = moveto[self.joiner[i][j]-1]
            i   = i + dir[0]
            j   = j + dir[1]
            verts.append([i+dl,j+dl])
            codes.append(Path.LINETO)
         path  = Path(verts,codes)
         patch = patches.PathPatch(path,lw=1.5,edgecolor='b',facecolor='none')
         axis.add_patch(patch)
         del verts,codes
  
      # between faces
      start   = [[ne,    2*ne-1],[2*ne,  2*ne-1],[2*ne-1,3*ne-1],[0,2*ne-1], \
                 [ne,   0],[3*ne, ne]]
      end     = [[2*ne-1,2*ne-1],[3*ne-1,2*ne-1],[ne,    3*ne-1],[0,ne],    \
                 [2*ne-1,0],[4*ne-1,ne]]
      for i in range(5):
         #verts[i] = [[end[i][0]+dl,end[i][1]+dl],[start[i+1][0]+dl,start[i+1][1]+dl]]
         #codes[i] = [Path.MOVETO,Path.LINETO]
         xfac = 0.5
         yfac = 0.5
         if i==0: xcurv,ycurv = +0,+0
         if i==1: xcurv,ycurv = +1,+1
         if i==2: xcurv,ycurv = -1,+1
         if i==3: xcurv,ycurv = -1,-1
         if i==4: xcurv,ycurv = +1,-1
         xcurv,ycurv = xcurv*ne*3.0/8.0,ycurv*ne*3.0/8.0
         center = [abs(start[i+1][0]+end[i][0])*xfac+xcurv, \
                   abs(start[i+1][1]+end[i][1])*yfac+ycurv]
         verts  = [[end[i][0]+dl,end[i][1]+dl],[center[0]+dl,center[1]+dl], \
                   [start[i+1][0]+dl,start[i+1][1]+dl]]
         codes  = [Path.MOVETO,Path.CURVE3,Path.CURVE3]
         path   = Path(verts,codes)
         patch  = patches.PathPatch(path,lw=1.5,edgecolor='b',facecolor='none')
         axis.add_patch(patch)



   def get_quadrature_points_lonlat(self, iface, a4, b4, np, qp):

      a = [0.0 for i in range(np)]
      b = [0.0 for i in range(np)]
      lons = []
      lats = []
      an = a4[1]
      a0 = a4[0]
      bn = b4[2]
      b0 = b4[1]
      for i in range(np):
         a[i] = -(qp[np-1]-qp[i])*(an-a0)/(qp[np-1]-qp[0])+an
         b[i] = -(qp[np-1]-qp[i])*(bn-b0)/(qp[np-1]-qp[0])+bn
 
      for i in range(np):
         for j in range(np):
            lon, lat = self.alphabeta_to_lonlat(iface,a[i],b[j])
            lon = lon/pi*180.0
            lat = lat/pi*180.0
            lons.append(lon)
            lats.append(lat)

      return lons, lats
 


   #########################################
   # high level APIs
   #########################################
 
   def draw_cube_structure(self, title):
      pltfig, axis = self.plot.create_canvas(title)
  
      self.plot.draw_box(axis,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.grb,upper=False,iscolor=False,flag1=self.flag1)
      self.plot.draw_points_in_box(axis,var=self.iup,kind=self.isu,flag1=self.flag1)
      self.draw_face_number(axis)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 
 
   def draw_cube_cart_glb_index(self, title):
      pltfig, axis = self.plot.create_canvas(title)
  
      self.plot.draw_box(axis,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.cij,upper=True,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.grb,upper=False,flag1=self.flag1)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 
 
   def draw_cube_glb_sfc_index(self, title):
      pltfig, axis = self.plot.create_canvas(title)
  
      self.plot.draw_box(axis,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.grb,upper=True,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.isfc,upper=False,flag1=self.flag1)
      self.draw_sfc_line(axis)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 
 
   def draw_sfc_rank_dist(self, title):
      pltfig, axis = self.plot.create_canvas(title)
  
      self.plot.draw_box(axis,color=self.color,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.isfc,upper=True,iscolor=True,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.rank,upper=False,iscolor=True,flag1=self.flag1)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 
 
   def draw_sfc_lcl_index(self, title):
      pltfig, axis = self.plot.create_canvas(title)
  
      self.plot.draw_box(axis,color=self.color,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.isfc,upper=True,iscolor=True,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.lcl,upper=False,iscolor=True,flag1=self.flag1)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 
 
   def draw_domain(self, title):
      pltfig, axis = self.plot.create_canvas('')
  
      self.plot.draw_box(axis,color=self.color, \
                                    flag1=self.flag1,flag3=self.flag3)
                                    #flag1=self.flag1,flag3=self.flag3,flag4=self.flag4)
      self.plot.write_values_in_box(axis,self.lcl,upper=False,iscolor=True, \
                                    flag1=self.flag1,flag3=self.flag3)
                                    #flag1=self.flag1,flag3=self.flag3,flag4=self.flag4)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 

   # domain decomposition
   def draw_local_domain(self, title):
      pltfig, axis = self.plot.create_canvas('')
  
      self.plot.draw_box(axis,color=self.color,\
                                    flag1=self.flag1,flag3=self.flag3,flag4=self.fflag)
      self.plot.write_values_in_box(axis,self.grb,upper=False,iscolor=True, \
                                    flag1=self.flag1,flag3=self.flag3,flag4=self.fflag)
      #self.plot.draw_points_in_box(axis,var=self.iup,kind=self.isu,flag1=self.flag1)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 
 
   def draw_cubed_sphere(self, title, draw_gll=False):
 
      ne    = self.ne
      nface = self.nface
      pnx   = self.pnx
      pny   = self.pny
  
      pltfig, axis = self.plot.create_canvas('')
  
      # seoul (127E, 38.0N)
      #lon0 =  127.0
      #lat0 =   0.0
      #lon0 = 127.0
      #lat0 =  38.0
      lon0 =  95.0
      lat0 =  40.0
      # resolution: c<l<i<h<f
      map = Basemap(projection='ortho',lon_0=lon0,lat_0=lat0,resolution='c')  # seoul
   
      # draw coastlines, country boundaries, fill continents.
      map.drawcoastlines(linewidth=0.8)
      #map.drawcountries(linewidth=0.8)
      
      #map.fillcontinents(color='coral',lake_color='aqua')
      map.fillcontinents(color='black',alpha=0.1,lake_color='white')
      #map.fillcontinents(color='white',lake_color='white')
      
      # draw the edge of the map projection region (the projection limb)
      #map.drawmapboundary(fill_color='aqua')
      map.drawmapboundary(fill_color='white')
      
      # draw lat/lon grid lines every 30 degrees.
      #map.drawmeridians(np.arange(0,360,30))
      #map.drawparallels(np.arange(-90,90,30))
   
      lons = [0.0 for i in range(4)]
      lats = [0.0 for i in range(4)]
   
      # draw cubed-sphere line
      for i in range(pnx):
         for j in range(pny):
            if self.iface[i][j] != -1:
               # check edge
               edge = [False for k in range(4)]
               if j==0:
                  edge[0] = True
               else:
                  if self.iface[i][j-1]!=self.iface[i][j]: edge[0] = True
               if i==pnx-1:
                  edge[1] = True
               else:
                  if self.iface[i+1][j]!=self.iface[i][j]: edge[1] = True
               if j==pny-1:
                  edge[2] = True
               else:
                  if self.iface[i][j+1]!=self.iface[i][j]: edge[2] = True
               if i==0:
                  edge[3] = True
               else:
                  if self.iface[i-1][j]!=self.iface[i][j]: edge[3] = True
     
               # draw element line
               #for k in range(4):
               #   lons[k] = self.lon[i][j][k]/pi*180.0
               #   lats[k] = self.lat[i][j][k]/pi*180.0
               '''
               x, y = map(lons,lats)
               axis.scatter(x,y,s=1)
               map.drawgreatcircle(lons[0],lats[0],lons[1],lats[1],linewidth=2,color='k')
               map.drawgreatcircle(lons[1],lats[1],lons[2],lats[2],linewidth=2,color='k')
               map.drawgreatcircle(lons[2],lats[2],lons[3],lats[3],linewidth=2,color='k')
               map.drawgreatcircle(lons[3],lats[3],lons[0],lats[0],linewidth=2,color='k')
               '''
     
               npts  = 500/ne
               lons, lats = self.get_alphabeta_square(self.iface[i][j],npts, \
                                                 self.alpha[i][j],self.beta[i][j])
               for k in range(4):
                  for ip in range(npts):
                     lons[k][ip] = lons[k][ip]/pi*180.0
                     lats[k][ip] = lats[k][ip]/pi*180.0
                  x,y = map(lons[k],lats[k])
                  if edge[k]:
                     axis.scatter(x,y,s=8.0,c='r',marker='.',linewidth=0.0)
                  else:
                     axis.scatter(x,y,s=4.0,c='b',marker='.',linewidth=0.0)


      # draw quadrature point
      if draw_gll:

         lons = []
         lats = []
         for i in range(pnx):
            for j in range(pny):
               if self.iface[i][j] != -1:
                  lons_np,lats_np = self.get_quadrature_points_lonlat(self.iface[i][j], \
                                    self.alpha[i][j],self.beta[i][j],self.np,self.qp)
                  lons.extend(lons_np)
                  lats.extend(lats_np)
         x,y = map(lons,lats)
         axis.scatter(x,y,s=50.0,c='k',marker='.',linewidth=0.0)
         del lons, lats


      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis

 
 
   # domain decomposition
   def draw_one_elem(self, title):
 
      pltfig, axis = self.plot.create_canvas('')
      self.plot.draw_box(axis,color=self.color,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.grb,upper=False,iscolor=True,flag1=self.flag1)
      isu = 0
      self.isu[1][0][3][0] = isu
      self.isu[1][0][3][1] = isu
      self.isu[1][0][3][2] = isu
      self.isu[1][0][3][3] = isu
      self.isu[1][0][0][3] = isu
      self.isu[1][0][1][3] = isu
      self.isu[1][0][2][3] = isu
      self.isu[1][0][3][3] = isu
      self.isu[0][1][0][3] = isu
      self.isu[0][1][1][3] = isu
      self.isu[0][1][2][3] = isu
      self.isu[0][1][3][3] = isu
      self.isu[3][1][0][0] = isu
      self.isu[3][1][1][0] = isu
      self.isu[3][1][2][0] = isu
      self.isu[3][1][3][0] = isu
      self.isu[1][2][0][0] = isu
      self.isu[1][2][1][0] = isu
      self.isu[1][2][2][0] = isu
      self.isu[1][2][3][0] = isu
      self.plot.draw_points_in_box(axis,var=self.iup,kind=self.isu,bigger=True,flag1=self.flag1)

      self.draw_info_legend(axis,True)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis


   def draw_info_legend(self, axis, extend=False):

      np = self.np
      ne = self.ne
      base    = [ 3.0, 3.0]
      tpos    = [-1.0,-1.0]
      tsize   = 20.
      dy      = 0.15

      #fontstyle = 'normal'
      #fontstyle = 'italic'
      fontstyle = 'oblique'

      tpos[0] = ne*base[0]
      tpos[1] = ne*base[1]-dy

      axis.text(tpos[0],tpos[1],'np     = {}'.format(self.np),color='k', \
                                size=tsize,va='center',style=fontstyle)

      tpos[1] = tpos[1]-dy
      axis.text(tpos[0],tpos[1],'ne     = {}'.format(self.ne),color='k', \
                                size=tsize,va='center',style=fontstyle)

      tpos[1] = tpos[1]-dy
      axis.text(tpos[0],tpos[1],'nprocs = {}'.format(self.nprocs),color='k', \
                                size=tsize,va='center',style=fontstyle)


      if extend:
         tpos[1] = tpos[1]-dy
         axis.text(tpos[0]+0.3,tpos[1],': EP',color='k', \
                                   size=tsize,va='center',style=fontstyle)
         circle = plt.Circle((tpos[0],tpos[1]),radius=0.035,color='k')
         axis.add_patch(circle)
         circle = plt.Circle((tpos[0]+0.2,tpos[1]),radius=0.035,color='b')
         axis.add_patch(circle)
   
         tpos[1] = tpos[1]-dy
         axis.text(tpos[0]+0.3,tpos[1],': UP',color='k', \
                                   size=tsize,va='center',style=fontstyle)
         circle = plt.Circle((tpos[0]+0.2,tpos[1]),radius=0.035,color='b')
         axis.add_patch(circle)

 
 
 
   #########################################
   # paper (2018)
   #########################################
 
   def create_canvas_2_2(self, title):

      fig_xsize = self.plot.canv_ysize/1.5
      fig_ysize = self.plot.canv_ysize/1.5
      pltfig = plt.figure(figsize=(fig_xsize,fig_ysize))
      margin = 0.01
      pltfig.subplots_adjust(left=margin,bottom=margin,right=1.0-margin, \
                             top=1.0-margin-0.04,wspace=0.1,hspace=0.01)
  
      pltfig.clear()
      axis  = pltfig.add_subplot(1,1,1)
  
      ### domain size
      margin = 0.2
      xmin = -0.1-margin+2
      xmax = 2+3+margin
      ymin = -0.1-margin+2
      ymax = 2+3+margin
      self.plot.axis_init(axis,xmin,xmax,ymin,ymax)
      axis.set_title(title,fontsize=25)
      return pltfig, axis
 

 
   # sphere
   def draw_paper_figure_1_1(self, title):
      self.draw_cubed_sphere(title)
 
 
   # cube (element and face number)
   def draw_paper_figure_1_2(self, title):
      pltfig, axis = self.plot.create_canvas('')
  
      self.plot.draw_box(axis,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.grb,upper=False,iscolor=True,flag1=self.flag1)
      self.draw_face_number(axis)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 
   # unique point (element and unique number)
   def draw_paper_figure_2_1(self, title):
      pltfig, axis = self.plot.create_canvas('')
  
      self.plot.draw_box(axis,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.grb,upper=False,iscolor=True,flag1=self.flag1)
      self.plot.draw_points_in_box(axis,var=self.iup,kind=self.isu,flag1=self.flag1)
      self.draw_info_legend(axis)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 
   # unique point (element and unique number)
   def draw_paper_figure_2_2(self, title):
      pltfig, axis = self.create_canvas_2_2('')
  
      self.plot.draw_box(axis,flag1=self.flag1,flag2=self.flag2)
      self.plot.write_values_in_box(axis,self.grb,upper=False,iscolor=True, \
                                   flag1=self.flag1,flag2=self.flag2)
      self.plot.draw_points_in_box(axis,var=self.iup,kind=self.isu, \
                                       flag1=self.flag1,flag2=self.flag2)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 
   # space filling curve (sfc index and unique number)
   def draw_paper_figure_3_1(self, title):
      pltfig, axis = self.plot.create_canvas('')
  
      self.plot.draw_box(axis,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.grb,upper=False,iscolor=True,flag1=self.flag1)
      self.plot.draw_points_in_box(axis,var=self.iup,kind=self.isu,flag1=self.flag1)
      self.draw_sfc_line(axis)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 

   # domain decomposition
   def draw_paper_figure_3_2(self, title):
      pltfig, axis = self.plot.create_canvas('')
  
      self.plot.draw_box(axis,color=self.color,flag1=self.flag1)
      self.plot.write_values_in_box(axis,self.grb,upper=False,iscolor=True,flag1=self.flag1)
      self.plot.draw_points_in_box(axis,var=self.iup,kind=self.isu,flag1=self.flag1)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 
 
   def draw_unique(self, title, opt, kind=False):
      pltfig, axis = self.plot.create_canvas('')
  
      if kind:
         flag3 = self.flag3
      else:
         flag3 = False
     
      self.plot.draw_box(axis,color=self.color,flag1=self.flag1,flag3=flag3)

      self.plot.write_values_in_box(axis,self.grb,upper=False,iscolor=True, \
                                   flag1=self.flag1,flag3=flag3)

      if opt>0:
         self.plot.draw_points_in_box(axis,bigger=True,var=self.value[opt-1], \
                           kind=kind,size=0.7,flag1=self.flag1,flag3=flag3)
      else:
         self.plot.draw_points_in_box(axis,bigger=True,var=self.iup, \
                           kind=kind,size=0.7,flag1=self.flag1,flag3=flag3)

      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 
 
   def draw_unique_org(self, title, opt=1, extend=False):
      pltfig, axis = self.plot.create_canvas('')
  
      if extend:
         flag3 = self.flag3
      else:
         flag3 = False
     
      self.plot.draw_box(axis,color=self.color,flag1=self.flag1,flag3=flag3)
      self.plot.write_values_in_box(axis,self.grb,upper=False,iscolor=True, \
                                   flag1=self.flag1,flag3=flag3)
      
      isu = [[[[self.isu[ip][jp][i][j] for j in range(self.np)] for i in range(self.np)]  \
                                   for jp in range(self.pny)] for ip in range(self.pnx)]
      for ip in range(self.pnx):
         for jp in range(self.pny):
            for i in range(self.np):
               for j in range(self.np):
                  isu[ip][jp][i][j] = self.isu[ip][jp][i][j]
                  if isu[ip][jp][i][j] > 1 and (not extend):
                     isu[ip][jp][i][j] = 0


      bigger = True
      if extend:
         self.plot.draw_points_in_box(axis,bigger=bigger, \
                  var=self.value[opt-1],kind=isu,size=0.7,flag1=self.flag1,flag3=flag3)
      else:
         self.plot.draw_points_in_box(axis,bigger=bigger, \
                  var=self.value[opt-1],kind=False,size=0.7,flag1=self.flag1,flag3=flag3)
  
      pltfig.savefig(self.plot.figs_dir+'/'+title,dpi=self.dpi)
      del axis
 
 
 
   def check_result(self):

      issame = self.compare_two_values(self.value[1],self.value[3])
      if issame:
         print 'LUP: Two values are same!'
      else:
         print 'LUP: Two values are different!'
      print('----')
  
      issame = self.compare_two_values(self.value[1],self.value[5])
      if issame:
         print 'EUP: Two values are same!'
      else:
         print 'EUP: Two values are different!'
      print('----')
 
 
 
   def compare_two_values(self, v1, v2):
 
      small = 1e-13

      issame = True
      for ie in range(self.pnx):
         for je in range(self.pny):
            if self.iface[ie][je]!=-1:
               for i in range(self.np):
                  for j in range(self.np):
                     if self.iface[ie][je]!=-1 and self.isu[ie][je][i][j]>0:
                        #if v1[ie][je][i][j]!=v2[ie][je][i][j]:
                        if abs(v1[ie][je][i][j]-v2[ie][je][i][j])>small:
                           issame = False
                           print 'ie={:2d}, je={:2d}, i={:2d}, j={:2d}:: v1={:20.17f}, v2={:20.17f}'.format \
                            (ie+1,je+1,i+1,j+1,v1[ie][je][i][j],v2[ie][je][i][j])
                           break
      return issame
