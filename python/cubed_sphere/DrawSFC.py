#!/usr/bin/env python
from numpy import *
from scipy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
import matplotlib.patches as patches
import time
from matplotlib import cm
import os
import sys
# load My Class
SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
sys.path.append(SHARE_DIR)
from LibSpaceFillingCurve import *
from LibControlMesh import *


def IniAxis(axis_in,xmin_in,xmax_in,ymin_in,ymax_in):
   axis_in.clear()
   axis_in.set_frame_on(False)
   axis_in.set_xlim(xmin_in,xmax_in)
   axis_in.set_ylim(ymin_in,ymax_in)
   axis_in.set_yticklabels([])
   axis_in.set_xticklabels([])
   axis_in.set_yticks([])
   axis_in.set_xticks([])

def CoordToBox(i_in, j_in, color_in, dx_in, ddxy_in):
   dl = dx_in-ddxy_in
   si = float(i_in)*dx_in
   ei = float(i_in)*dx_in + dl
   sj = float(j_in)*dx_in
   ej = float(j_in)*dx_in + dl

   #verts = [(si,sj), (si,ej), (ei,ej), (ei,sj), (si,sj)]
   verts = [[si,sj], [si,ej], [ei,ej], [ei,sj], [si,sj]]
   codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
   path = Path(verts, codes)

   return patches.PathPatch(path, facecolor=color_in, lw=2.5)


# Main Draw Algorithm
def DrawCubeStructure(axis_in, dx_in, ddxy_in, pos, color, cpos_in, isfc_in, tsize_in):
   ii, jj   = pos[0], pos[1]
   ic, jc   = cpos_in[0]-0.2, cpos_in[1]-0.2
   patch    = CoordToBox(ii, jj, color, dx_in, ddxy_in)
   axis_in.add_patch(patch)

   iit, jjt = dx_in*(ii+ic), dx_in*(jj+jc)
   axis.text(iit, jjt, str(isfc_in[ii][jj]), color='k', size=tsize_in, va='center', ha='center')

# Main Draw Algorithm
def DrawCubeStructure_old(axis_in, pos, color, ddxy_in):
   ii, jj = pos[0], pos[1]
   patch  = CoordToBox(ii, jj, color, ddxy_in)
   axis_in.add_patch(patch)

def DrawLineSFC(axis_in, ne_in, prim_in, joiner_in, dx_in, ddxy_in, larrow_in):
   moveto = [[-1,0], [1,0], [0,-1], [0,1]]
   dl     = (dx_in-ddxy_in)/dx/2.0
   verts  = []
   codes  = []

   start = [[0,0], [ne_in-1,ne_in-1], [0,0], [ne_in-1,ne_in-1]]
   i = start[prim_in][0]
   j = start[prim_in][1]

   verts.append([dx_in*(i+dl),dx_in*(j+dl)])
   codes.append(Path.MOVETO)
   dir = [0,0]
   size_arrow = dx_in*0.3
   for ie in range((ne_in*ne_in)-1):
      # arrow
      if larrow_in:
         joiner = joiner_in[i][j]
         if joiner == 0:
            lr = -size_arrow; ud = 0.0
         elif joiner == 1:
            lr =  size_arrow; ud = 0.0
         elif joiner == 2:
            lr = 0.0; ud = -size_arrow
         elif joiner == 3:
            lr = 0.0; ud =  size_arrow
         arrow = patches.Arrow(dx_in*(i+dl), dx_in*(j+dl), lr, ud, color='r', width=0.2)
         axis_in.add_patch(arrow)
      # for lines
      dir = moveto[joiner_in[i][j]]
      i   = i + dir[0]
      j   = j + dir[1]
      verts.append([dx_in*(i+dl),dx_in*(j+dl)])
      codes.append(Path.LINETO)
   path  = Path(verts, codes)
   patch = patches.PathPatch(path, lw=2.0, edgecolor='b', facecolor='none')
   axis_in.add_patch(patch)



######### Main Program ##########

ne     = 10
iprim  = 0
larrow = False

sfc    = SFC(ne)
isfci  = sfc.GenSFC()       # SFC index
jsfci  = sfc.GetJoiner()    # Joiner vectors

mymesh = Mesh()
xaxis, yaxis, xyaxis, mxyaxis, orgaxis = mymesh.GetAxisIndex()
if iprim == 0:
  isfc = isfci
  jsfc = jsfci
elif iprim == 1:
  isfc = mymesh.SymmetryMesh(ne, orgaxis, isfci)
  jsfc = mymesh.SymmetryDirMesh(ne, orgaxis, jsfci)
elif iprim == 2:
  isfc = mymesh.SymmetryMesh(ne, xyaxis, isfci)
  jsfc = mymesh.SymmetryDirMesh(ne, xyaxis, jsfci)
elif iprim == 3:
  isfc = mymesh.SymmetryMesh(ne, mxyaxis, isfci)
  jsfc = mymesh.SymmetryDirMesh(ne, mxyaxis, jsfci)


filename = 'FigSFC/SFC_n'+str(ne)+'_prim'+str(iprim)+'_'+str(larrow)+'.png'
print filename

### Plot Domain
dx    = 1.0   # width of boxes
ddxy  = 0.0   # distance between boxes
nrows = ne
ncols = ne

### Domain Size
margin = 0.2
xmin = -margin
xmax = dx*float(nrows)+2.0*margin
ymin = -margin
ymax = dx*float(ncols)+2.0*margin

### Text attributes
it, jt = +0.5*(dx-ddxy)/dx, +0.5*(dx-ddxy)/dx
tsize  = 30./sqrt(float(ne))

### Ploat
pltfig = plt.figure(figsize=(nrows*5.0/ne,ncols*5.0/ne))

pltfig.clear()
axis  = pltfig.add_subplot(1,1,1)
IniAxis(axis,xmin,xmax,ymin,ymax)
axis.set_title('Space-Filling Curve',fontsize =25)

for i in range(nrows):
   for j in range(ncols):
      DrawCubeStructure(axis, dx, ddxy, [i,j], 'none', [it, jt], isfc, tsize)
DrawLineSFC(axis, ne, iprim, jsfc, dx, ddxy, larrow)

pltfig.savefig(filename)

plt.show(True)
