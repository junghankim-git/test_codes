#!/usr/bin/env python
from numpy import *
from scipy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_toolkits.basemap import Basemap
import time
from matplotlib import cm
import os
import sys
#font
import matplotlib.font_manager as fontm


# load My Class
#SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
#sys.path.append(SHARE_DIR)
from LibSpaceFillingCurve import *
from LibControlMesh import *
from LibCubedSphere import *

os.system('rm -rf ./FigCube')
os.system('mkdir ./FigCube')

ndir = 8
mymesh = Mesh()
ww, ee, ss, nn, ws, es, wn, en = mymesh.GetDirectionIndex()
moveto                         = [[-1,0], [1,0], [0,-1], [0,1]]




####################################################
###                Functions                     ###
####################################################

def IniAxis(axis_in,xmin_in,xmax_in,ymin_in,ymax_in):
   axis_in.clear()
   axis_in.set_frame_on(False)
   axis_in.set_xlim(xmin_in,xmax_in)
   axis_in.set_ylim(ymin_in,ymax_in)
   axis_in.set_yticklabels([])
   axis_in.set_xticklabels([])
   axis_in.set_yticks([])
   axis_in.set_xticks([])


def CoordToBox(i_in, j_in, color_in):
   dx = 1.0
   dl = 0.9*dx
   si = i_in
   ei = i_in + dl
   sj = j_in
   ej = j_in + dl

   #verts = [(si,sj), (si,ej), (ei,ej), (ei,sj), (si,sj)]
   verts = [[si,sj], [si,ej], [ei,ej], [ei,sj], [si,sj]]
   codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
   path = Path(verts, codes)

   return patches.PathPatch(path, facecolor=color_in, lw=1.5)



def GetAlphaBetaLine(cube_in, iface_in, npts_in, alpha_beta1, alpha_beta2):
   res_line = [[0.0 for i in range(npts_in)] for j in range(2)]
   a1 = alpha_beta1[0]
   b1 = alpha_beta1[1]
   a2 = alpha_beta2[0]
   b2 = alpha_beta2[1]

   if (a1 == a2) and (b1 != b2):
      db = abs(b2-b1)/(npts_in-1)
      for ip in range(npts_in):
         beta = min(b1,b2)+ip*db
         res_line[0][ip], res_line[1][ip] = cube_in.AlphaBeta2LonLat(iface_in, a1, beta)
   elif (a1 != a2) and (b1 == b2):
      da = abs(a2-a1)/(npts_in-1)
      for ip in range(npts_in):
         alpha = min(a1,a2)+ip*da
         res_line[0][ip], res_line[1][ip] = cube_in.AlphaBeta2LonLat(iface_in, alpha, b1)
   else:
      print 'Error : GetAlphaBetaLine...'
      quit()
   return res_line


def GetAlphaBetaSquare(cube_in, iface_in, npts_in, alpha_beta):
   res_lines = [[[0.0 for i in range(npts_in)] for j in range(2)] for k in range(4)]
   for k in range(4):
      if k == 3:
         res_lines[k] = GetAlphaBetaLine(cube_in, iface_in, npts_in, alpha_beta[k], alpha_beta[0])
      else:
         res_lines[k] = GetAlphaBetaLine(cube_in, iface_in, npts_in, alpha_beta[k], alpha_beta[k+1])
   return res_lines




def FaceNumber(axis_in, ne_in):
   ## start font
   #fname = []
   #for flist in fontm.fontManager.afmlist:
   #   fname.append(flist.name)
   ## end font
   dx = 2
   nface = 6
   base  = [[1.0,1.0], [2.0,1.0], [3.0,1.0], [0.0,1.0], [1.0,0.0], [1.0,2.0]]
   tpos  = [[-1,-1] for i in range(nface)]
   #tpos  = ne_in*(base+0.5)
   for iface in range(nface):
      tpos[iface][0]  = ne_in*(base[iface][0]+0.5)
      tpos[iface][1]  = ne_in*(base[iface][1]+0.5)
   for iface in range(nface):
      axis_in.text(tpos[iface][0], tpos[iface][1], str(iface+1), color='k', size=150.0, alpha=0.3, weight='semibold', fontname='sans-serif', style='italic', ha='center', va='center')
      axis_in.text(tpos[iface][0], tpos[iface][1], str(iface+1), color='b', size=150.0, alpha=0.3, fontname='sans-serif', style='italic', ha='center', va='center')


def ViewSet(axis_in, ne_in, np_in, nprocs_in, vproc_in, optp_in, optvp_in):
   # basic position
   base    = [ 3.0, 2.0]
   tpos    = [-1.0,-1.0]
   # size
   tsize   = 25.0 #100./ne_in
   # postion
   tpos[0] = ne_in*(base[0]-0.1)
   tpos[1] = ne_in*(base[1]+0.75)
   # text
   axis_in.text(tpos[0], tpos[1], 'np = '+str(np_in), color='k', size=tsize, va='center', style='italic')#, fontname='sans-serif'
   tpos[1] = ne_in*(base[1]+0.60)
   axis_in.text(tpos[0], tpos[1], 'ne = '+str(ne_in), color='k', size=tsize, va='center', style='italic')#, fontname='sans-serif'
   if optp_in:
      tpos[1] = ne_in*(base[1]+0.45)
      axis_in.text(tpos[0], tpos[1], 'nprocs = '+str(nprocs_in), color='k', size=tsize, va='center', style='italic')#, fontname='sans-serif'
   if optvp_in:
      tpos[1] = ne_in*(base[1]+0.3)
      axis_in.text(tpos[0], tpos[1], 'iproc = '+str(vproc_in), color='k', size=tsize, va='center', style='italic')#, fontname='sans-serif'


def DrawLineSFC(axis_in, iface_in, ne_in, joiner_in):
   dl = 0.45
   start = [[ne_in,2*ne_in-1],     [2*ne_in,2*ne_in-1],   [3*ne_in,ne_in],   [0,2*ne_in-1], [ne_in,0],     [2*ne_in-1,3*ne_in-1]]
   end   = [[2*ne_in-1,2*ne_in-1], [3*ne_in-1,2*ne_in-1], [4*ne_in-1,ne_in], [0,ne_in],     [2*ne_in-1,0], [ne_in,3*ne_in-1]]
   verts = []
   codes = []
   i = start[iface_in][0]
   j = start[iface_in][1]
   verts.append([i+dl,j+dl])
   codes.append(Path.MOVETO)
   dir = [0,0]
   for ie in range((ne_in*ne_in)-1):
      dir = moveto[joiner_in[i][j]]
      i   = i + dir[0]
      j   = j + dir[1]
      verts.append([i+dl,j+dl])
      codes.append(Path.LINETO)
#   codes[-1] = Path.CLOSEPOLY
   path = Path(verts, codes)
   patch = patches.PathPatch(path, lw=1.5, edgecolor='b', facecolor='none')
   axis_in.add_patch(patch)
   #for i in range(6):
   axis_in.scatter(start[iface_in][0]+dl, start[iface_in][1]+dl, s=800.0/ne_in, marker ='>', c ='r')
   axis_in.scatter(end[iface_in][0]+dl, end[iface_in][1]+dl, s=800.0/ne_in, marker =(5,1), c ='g')


def DrawLineSFC_FacePath(axis_in, ne_in):
   dl = 0.45
   #                1 face                   2 face                3 face              4 face         5 face           6 face        
   #start   = [[ne_in,    2*ne_in-1], [2*ne_in,  2*ne_in-1], [3*ne_in,  ne_in], [0,2*ne_in-1], [ne_in,0],     [2*ne_in-1,3*ne_in-1]]
   #end     = [[2*ne_in-1,2*ne_in-1], [3*ne_in-1,2*ne_in-1], [4*ne_in-1,ne_in], [0,ne_in],     [2*ne_in-1,0], [ne_in,3*ne_in-1]]
   #                1 face                   2 face                6 face              4 face         5 face           3 face        
   start   = [[ne_in,    2*ne_in-1], [2*ne_in,  2*ne_in-1], [2*ne_in-1,3*ne_in-1], [0,2*ne_in-1], [ne_in,    0], [3*ne_in,  ne_in]]
   end     = [[2*ne_in-1,2*ne_in-1], [3*ne_in-1,2*ne_in-1], [ne_in,    3*ne_in-1], [0,ne_in],     [2*ne_in-1,0], [4*ne_in-1,ne_in]]
   for i in range(5):
      #verts[i] = [[end[i][0]+dl, end[i][1]+dl], [start[i+1][0]+dl, start[i+1][1]+dl]]
      #codes[i] = [Path.MOVETO, Path.LINETO]
      xfac = 0.5
      yfac = 0.5
      if i==0: xcurv, ycurv = +0, +0
      if i==1: xcurv, ycurv = +1, +1
      if i==2: xcurv, ycurv = -1, +1
      if i==3: xcurv, ycurv = -1, -1
      if i==4: xcurv, ycurv = +1, -1
      xcurv, ycurv = xcurv*ne_in*3.0/8.0, ycurv*ne_in*3.0/8.0
      center   = [abs(start[i+1][0]+end[i][0])*xfac+xcurv, abs(start[i+1][1]+end[i][1])*yfac+ycurv]
      verts = [[end[i][0]+dl, end[i][1]+dl], [center[0]+dl, center[1]+dl], [start[i+1][0]+dl, start[i+1][1]+dl]]
      codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
      path = Path(verts, codes)
      patch = patches.PathPatch(path, lw=1.5, edgecolor='b', facecolor='none')
      axis_in.add_patch(patch)
      #axis_in.scatter(center[0]+dl, center[1]+dl, s=9, marker ='o', c ='r')
    


'''
   verts12 = [[end[0][0]+dl, end[0][1]+dl], [start[1][0]+dl, start[1][1]+dl]]
   codes12 = [Path.MOVETO, Path.LINETO]
   path = Path(verts12, codes12)
   patch = patches.PathPatch(path, lw=1.5, edgecolor='b', facecolor='none')
   axis_in.add_patch(patch)
'''


# Main Draw Algorithm
def DrawCubeStructure(axis_in, pos, color, qx, qy):
   ii, jj = pos[0], pos[1]
   patch = CoordToBox(ii,jj,color)
   axis_in.add_patch(patch)
   axis_in.scatter(qx+ii, qy+jj, s=7, marker ='o', c ='r')



def DrawCubedSphere(axis_in, cube_in, nrows_in, ncols_in, iface_in, ab_in, lonlat_in):
   #axis_in.set_title('')
   lon0=127.0
   lat0=37.30
   map = Basemap(projection='ortho',lon_0=lon0,lat_0=lat0,resolution='l')  # seoul

   # draw coastlines, country boundaries, fill continents.
   map.drawcoastlines(linewidth=0.8)
   map.drawcountries(linewidth=0.8)
   
   #map.fillcontinents(color='coral',lake_color='aqua')
   #map.fillcontinents(color='white',lake_color='white')
   
   # draw the edge of the map projection region (the projection limb)
   #map.drawmapboundary(fill_color='aqua')
   map.drawmapboundary(fill_color='white')
   
   # draw lat/lon grid lines every 30 degrees.
   #map.drawmeridians(np.arange(0,360,30))
   #map.drawparallels(np.arange(-90,90,30))

   lons   = [0.0 for i in range(4)]
   lats   = [0.0 for i in range(4)]

   for i in range(nrows_in):
      for j in range(ncols_in):
         if iface_in[i][j] != -1:
#         if cidx[i][j][0]==3 and cidx[i][j][1]==4 and iface[i][j] == 6:
            for k in range(4):
               lons[k] = lonlat_in[i][j][k][0]/pi*180.0
               lats[k] = lonlat_in[i][j][k][1]/pi*180.0
            x, y = map(lons,lats)
            axis_in.scatter(x,y)

            npts  = 400/cube_in.ne
            lines = GetAlphaBetaSquare(cube_in, iface_in[i][j], npts, ab_in[i][j])
            for k in range(4):
               for ip in range(npts):
                  lines[k][0][ip] = lines[k][0][ip]/pi*180.0
                  lines[k][1][ip] = lines[k][1][ip]/pi*180.0
               x, y = map(lines[k][0], lines[k][1])
               axis_in.scatter(x,y,s=1)


def DrawGeneral1(axis_in, pos, upos, dpos, color_in, idx, tsize):
   ii, jj = pos[0], pos[1]
   iu, ju = upos[0], upos[1]
   id, jd = dpos[0], dpos[1]
   patch = CoordToBox(ii,jj,color_in)
   axis_in.add_patch(patch)
   axis_in.text(ii+id, jj+jd, str(idx), color='k', size=tsize, va='center', ha='center')


def DrawGeneral2(axis_in, pos, upos, dpos, color_in, idx1, idx2, tsize):
   ii, jj = pos[0], pos[1]
   iu, ju = upos[0], upos[1]
   id, jd = dpos[0], dpos[1]
   patch = CoordToBox(ii,jj,color_in)
   axis_in.add_patch(patch)
   if color_in == 'none':
      ucolor='k'
      dcolor='r'
   else:
      ucolor='w'
      dcolor='k'
   axis_in.text(ii+iu, jj+ju, str(idx1), color=ucolor, size=tsize*0.7, va='center', ha='center')#, label='%003d')
   axis_in.text(ii+id, jj+jd, str(idx2), color=dcolor, size=tsize, va='center', ha='center')#, label='%003d')


def ElementsPerProc(axis_in, pos, upos, dpos, color, iproc_in, vproc_in, nbr_in, idx, tsize, l_nbr=True, l_nbr_box=True):
   ii, jj = pos[0], pos[1]
   iu, ju = upos[0], upos[1]
   id, jd = dpos[0], dpos[1]
   if iproc_in == vproc_in:
      patch = CoordToBox(ii,jj,color)
      axis_in.add_patch(patch)
      axis_in.text(ii+id, jj+jd, str(idx), color='k', size=tsize, va='center', ha='center')
   else:
      if l_nbr:
         for idir in arange(8):
            cproc = nbr_in[idir]
            if cproc == vproc_in:
               patch = CoordToBox(ii,jj,color)
               axis_in.add_patch(patch)
               axis_in.text(ii+id, jj+jd, str(idx), color='k', size=tsize, va='center', ha='center')
            else:
               if l_nbr_box:
                  patch = CoordToBox(ii,jj,'none')
                  axis_in.add_patch(patch)
                  axis_in.text(ii+id, jj+jd, str(idx), color='k', size=tsize, va='center', ha='center')


def ElementsPerProc_Space(axis_in, pos, upos, dpos, color, iproc_in, vproc_in, nbr_in, idx, tsize, l_nbr=True, l_nbr_box=True):
   ii, jj = pos[0], pos[1]
   iu, ju = upos[0], upos[1]
   id, jd = dpos[0], dpos[1]
   if iproc_in == vproc_in:
      patch = CoordToBox(ii,jj,color)
      axis_in.add_patch(patch)
      axis_in.text(ii+id, jj+jd, str(idx), color='k', size=tsize, va='center', ha='center')
   else:
      if l_nbr:
         for idir in arange(8):
            cproc = nbr_in[idir]
            if cproc == vproc_in:
               patch = CoordToBox(ii,jj,color)
               axis_in.add_patch(patch)
               axis_in.text(ii+id, jj+jd, str(idx), color='k', size=tsize, va='center', ha='center')
            else:
               if l_nbr_box:
                  patch = CoordToBox(ii,jj,'none')
                  axis_in.add_patch(patch)
                  axis_in.text(ii+id, jj+jd, str(idx), color='k', size=tsize, va='center', ha='center')


def NewElementsPerProc(axis_in, pos, upos, dpos, color, ocolor, iproc_in, vproc_in, isOutin, nbr_in, idx, tsize):
   ii, jj = pos[0], pos[1]
   iu, ju = upos[0], upos[1]
   id, jd = dpos[0], dpos[1]
   if iproc_in == vproc_in:
      if isOutin == True:
         patch = CoordToBox(ii,jj,ocolor)
      else:
         patch = CoordToBox(ii,jj,color)
      axis_in.add_patch(patch)
      axis_in.text(ii+id, jj+jd, str(idx), color='k', size=tsize, va='center', ha='center')
   else:
      for idir in arange(8):
         cproc = nbr_in[idir]
         if cproc == vproc_in:
            patch = CoordToBox(ii,jj,color)
         else:
            patch = CoordToBox(ii,jj,'none')
         axis_in.add_patch(patch)
         axis.text(ii+id, jj+jd, str(idx), color='k', size=tsize, va='center', ha='center')


### Get Font's family
# [f.name for f in matplotlib.font_manager.fontManager.afmlist]
### Bold
# weight = ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']



####################################################
###                Main Program                  ###
####################################################
 
np      = 4
ne      = 2
nprocs  = 7

cube = CubedSphere(np, ne, nprocs)
cube.Ini()

doPlot = True

if doPlot:
   #---
   #### Quadrature Point Position 
   qp     = cube.qp
   qp[0]  = -0.94
   qp[-1] = +0.94
   x       = array([qp for i in arange(cube.np)])
   y       = array([[qp[i] for j in arange(cube.np)] for i in range(cube.np)])
   
   #---
   #### Text Colors
   # Method 1
#   icolors  = array(['#ff0000', '#00ff00', '#0000ff', '#ffff00', '#00ffff', '#ff00ff'])
#   iocolors = array(['#dd0000', '#00dd00', '#0000dd', '#dddd00', '#00dddd', '#dd00dd'])
#   colors  = []
#   ocolors = []
#   for i in range(30):
#      for j in range(6):
#         colors.append(icolors[j])
#         ocolors.append(iocolors[j])

   # Method 2
   colors  = []
   ocolors = []
   cstep  = ['ff', 'aa', '55', '00']
   #ocstep = ['bb', '77', '44', '11']
   ocstep = ['99', '77', '44', '11']
   for i in range(4):
      for j in range(3):
         colors.append('#'+cstep[i]+cstep[4-j-1]+cstep[4-j-1])
         colors.append('#'+cstep[4-j-1]+cstep[i]+cstep[4-j-1])
         colors.append('#'+cstep[4-j-1]+cstep[4-j-1]+cstep[i])
         colors.append('#'+cstep[i]+cstep[i]+cstep[4-j-1])
         colors.append('#'+cstep[4-j-1]+cstep[i]+cstep[i])
         colors.append('#'+cstep[i]+cstep[4-j-1]+cstep[i])
         ocolors.append('#'+ocstep[i]+ocstep[4-j-1]+ocstep[4-j-1])
         ocolors.append('#'+ocstep[4-j-1]+ocstep[i]+ocstep[4-j-1])
         ocolors.append('#'+ocstep[4-j-1]+ocstep[4-j-1]+ocstep[i])
         ocolors.append('#'+ocstep[i]+ocstep[i]+ocstep[4-j-1])
         ocolors.append('#'+ocstep[4-j-1]+ocstep[i]+ocstep[i])
         ocolors.append('#'+ocstep[i]+ocstep[4-j-1]+ocstep[i])
   
   #---
   #### Text Size
   txtsize = 80./float(ne)
    
   #---
   #### View Processor
   vproc = 1
   
   #---
   #### Text Postion
   #=[left, down]= old version
   #itU, jtU = +0.2, +0.63
   #itu, jtu = +0.3, +0.63
   #itd, jtd = +0.3, +0.3
   #[center, center]
   itu, jtu = +0.45, +0.75
   itd, jtd = +0.45, +0.45
    
   
   #---
   ### Plot Size
   nrows  = cube.nrows
   ncols  = cube.ncols
   
   #---
   ### Loading Data
   # old        0        1       2       3       4       5        6      7         8        9
   # ival = ['isOut', 'nlidx', 'iface', 'iproc', 'nbr', 'gidx', 'cidx', 'lidx', 'joiner', 'isfc']
   #           0         1       2         3      4      5        6       7      8        9        10       11
   # ival = ['isOut', 'nlidx', 'iface', 'iproc', 'nbr', 'gidx', 'cidx', 'ab', 'lonlat', 'lidx', 'joiner', 'isfc']
   iface  = cube.CubeToPlot(ival=2)
   gidx   = cube.CubeToPlot(ival=5)
   cidx   = cube.CubeToPlot(ival=6)
   joiner = cube.CubeToPlot(ival=10)
   isfc   = cube.CubeToPlot(ival=11)
   iproc  = cube.CubeToPlot(ival=3)
   nbridx = cube.CubeToPlot(ival=4)
   lidx   = cube.CubeToPlot(ival=9)
   nlidx  = cube.CubeToPlot(ival=1)
   isOut  = cube.CubeToPlot(ival=0)
   #mymesh.PrintMesh2D(nrows, ncols, joiner, 'joiner', True)
   ab     = cube.CubeToPlot(ival=7)
   lonlat = cube.CubeToPlot(ival=8)
   
   #---
   ### Domain Size
   margin = 0.2
   xmin = -0.1-margin
   xmax = nrows+margin
   ymin = -0.1-margin
   ymax = ncols+margin
   
   
   
   
   #------
   ### Plot
   
   nplots    = 11
   
   
   Titles    = ['Cube Element',                   'GridElem: Catesian, Global Index',  'GridElem: Global Index, Space Filling Curve', \
                'GridElem: Process Number',       'GridVertex: Global Index',          'MetaVertex(1): Local Index (all process)',    \
                'MetaVertex(1): Local Index',     'MetaVertex(1): Local Index(1)',     'MetaVertex(1): Local Index(2)', \
                'MetaVertex(1): New Local Index', '']
   
   OutDir    = 'FigCube/'
   Filenames = ['00_CubeElement.ps',                  '01_GridElem.ps',                 '02_GridElem_SFC.ps', \
                '03_GridElem_DomainDecomposition.ps', '04_GridVertex.ps',               '05_MetaVertex_allProcs.ps', \
                '06_MetaVertex.ps',                   '06_MetaVertex_1.ps',             '06_MetaVertex_2.ps', \
                '07_MetaVertex_NewLocalIndex.ps',     '08_CubedSphere.ps']
   
   
   print ' # START: PostProcessing          ...'
   
   pltfig = [plt.figure(figsize=(nrows*5.0/ne,ncols*5.0/ne)) for i in range(nplots)]
   
   
   for iplot in range(nplots):
      pltfig[iplot].clear()
      axis  = pltfig[iplot].add_subplot(1,1,1)
      IniAxis(axis,xmin,xmax,ymin,ymax)
      axis.set_title(Titles[iplot],fontsize =25)
      if iplot == 0:
         newx = 9./20.*(x+1.0)
         newy = 9./20.*(y+1.0)
      for i in range(nrows):
         for j in range(ncols):
            if iface[i][j] != -1:
               if iplot == 0:
                  DrawCubeStructure(axis, [i,j], 'none', newx, newy)
               elif iplot == 1:
                  DrawGeneral2(axis, [i,j], [itu,jtu], [itd,jtd], 'none', cidx[i][j], gidx[i][j], txtsize)
               elif iplot == 2:
                  DrawGeneral2(axis, [i,j], [itu,jtu], [itd,jtd], 'none', gidx[i][j], isfc[i][j], txtsize)
               elif iplot == 3:
                  DrawGeneral2(axis, [i,j], [itu,jtu], [itd,jtd], colors[iproc[i][j]-1], isfc[i][j], iproc[i][j], txtsize)
               elif iplot == 4:
                  DrawGeneral1(axis, [i,j], [itu,jtu], [itd,jtd], colors[iproc[i][j]-1], gidx[i][j], txtsize)
               elif iplot == 5:
                  DrawGeneral2(axis, [i,j], [itu,jtu], [itd,jtd], colors[iproc[i][j]-1], gidx[i][j], lidx[i][j], txtsize)
               elif iplot == 6:
                  ElementsPerProc(axis, [i,j], [itu,jtu], [itd,jtd], colors[iproc[i][j]-1], iproc[i][j], vproc, nbridx[i][j], lidx[i][j], txtsize)
               elif iplot == 7:
                  ElementsPerProc(axis, [i,j], [itu,jtu], [itd,jtd], colors[iproc[i][j]-1], iproc[i][j], vproc, nbridx[i][j], lidx[i][j], txtsize, True, False)
               elif iplot == 8:
                  ElementsPerProc(axis, [i,j], [itu,jtu], [itd,jtd], colors[iproc[i][j]-1], iproc[i][j], vproc, nbridx[i][j], lidx[i][j], txtsize, False, False)
               elif iplot == 9:
                  NewElementsPerProc(axis, [i,j], [itu,jtu], [itd,jtd], colors[iproc[i][j]-1], ocolors[iproc[i][j]-1], iproc[i][j], vproc, isOut[i][j], nbridx[i][j], nlidx[i][j], txtsize)
      if iplot == 10:
         DrawCubedSphere(axis, cube, nrows, ncols, iface, ab, lonlat)
      # Draw Face Number
      if iplot == 0 or iplot == 1: FaceNumber(axis, ne)
      # ne, np, nprocs options
      optp  = False
      optvp = False
      if iplot >= 3:
         optp  = True
      if iplot >= 6:
         optvp = True
      if iplot != 10:
         ViewSet(axis, ne, np, nprocs, vproc, optp, optvp)
      # Save Figure
      #pltfig[iplot].suptitle(Titles[iplot], fontsize =18)
      pltfig[iplot].savefig(OutDir+Filenames[iplot])
      # SFC Line
      if iplot == 2:
         for iiface in range(cube.nface):
            DrawLineSFC(axis, iiface, ne, joiner)
            if iiface == 0: pltfig[iplot].savefig(OutDir+'02_GridElem_SFC_1Face.ps')
         DrawLineSFC_FacePath(axis, ne)
         pltfig[iplot].savefig(OutDir+'02_GridElem_SFC_6Face.ps')
      print '   - Saving output file: '+str(Filenames[iplot])
   
   
   print ' # END  : PostProcessing          ...\n'
   plt.show(True)
