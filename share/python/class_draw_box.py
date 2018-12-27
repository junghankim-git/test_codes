#!/usr/bin/env python
import os
import sys
import time
import numpy
import scipy
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDrawingArea




class class_draw_box:



   def __init__(self, np, pts, nx, ny):

      self.np = np
      self.nx = nx
      self.ny = ny

      # box dx, dl
      self.dx = 1.0  # do not modify
      self.dl = 0.9  # dl <= dx
  
      # points in the box: -1 <= pts <= 1
      self.pts_1d = pts
      y = numpy.array([self.pts_1d for i in range(np)])
      x = numpy.array([[self.pts_1d[i] for j in range(np)] for i in range(np)])
      self.pts_x = self.dl*(x+1.0)/2.0
      self.pts_y = self.dl*(y+1.0)/2.0
      del x, y
   
      # text size (box and point)
      self.tsize    = 192./float(ny)   # text size in center of box
      self.tsize_up = 168./float(ny)   # text size in up of box
      self.tsize_pt = self.tsize/3.2   # text size in point of box

      # text postion (center)
      self.tpos    = [+0.45, +0.45]  # position of center text
      self.tpos_up = [+0.45, +0.75]  # position of up text
      
      # text colors
      self.colors  = []
      self.ocolors = []
      cstep  = ['ff', 'aa', '55', '00']
      ocstep = ['99', '77', '44', '11']
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
       
      # domain size
      self.canv_ysize = 9.0
      margin = 0.2
      self.xmin = -0.1-margin
      self.xmax = nx+margin
      self.ymin = -0.1-margin
      self.ymax = ny+margin
      self.figs_dir = './fig'



   def get_box_color(self, kind):
      color = numpy.empty(shape=(self.nx,self.ny)).tolist()
      for i in range(self.nx):
         for j in range(self.ny):
            if kind[i][j] >= 0:
               color[i][j] = self.colors[kind[i][j]]
      return color



   # main draw algorithm
   def create_canvas(self, title, xmin=None, xmax=None, ymin=None, ymax=None):
      fig_xsize = self.canv_ysize*self.nx/self.ny
      fig_ysize = self.canv_ysize
      pltfig = plt.figure(figsize=(fig_xsize,fig_ysize))
      margin = 0.000
      pltfig.subplots_adjust(left=margin,bottom=margin,right=1.0-margin, \
                             top=1.0-margin-0.04,wspace=0.01,hspace=0.01)
  
      pltfig.clear()
      axis  = pltfig.add_subplot(1,1,1)
  
      if xmin==None: xmin=self.xmin
      if xmax==None: xmax=self.xmax
      if ymin==None: ymin=self.ymin
      if ymax==None: ymax=self.ymax
      self.axis_init(axis,xmin,xmax,ymin,ymax)
      axis.set_title(title,fontsize=25)

      return pltfig, axis
 
 
 
   def show(self):
      plt.show()
 
 
 
   def axis_init(self, axis, xmin_in, xmax_in, ymin_in, ymax_in):
      axis.clear()
      axis.set_frame_on(False)
      axis.set_xlim(xmin_in,xmax_in)
      axis.set_ylim(ymin_in,ymax_in)
      axis.set_yticklabels([])
      axis.set_xticklabels([])
      axis.set_yticks([])
      axis.set_xticks([])
   
   
 
   def coord_to_box(self, ii, jj, color):

      #    --------- ---------
      #    |       | |       |
      #  dl|       | |       |dl
      #    |  dl   | |  dl   |
      #    --------- ---------
      # (ie,je)  (ie+1,je)
      dx = self.dx
      dl = self.dl
      si = ii
      ei = ii+dl
      sj = jj
      ej = jj+dl
   
      verts = [[si,sj],[si,ej],[ei,ej],[ei,sj],[si,sj]]
      codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY]
      path = Path(verts,codes)
      if color=='none':
         lw    = 1.5
         alpha = 1.0
      else:
         lw    = 0.0
         alpha = 0.5

      return patches.PathPatch(path,facecolor=color,lw=lw,alpha=alpha)



   def point_position_in_box(self, x, y, np, xx, yy):

      ix = [[0.0 for j in range(np)] for i in range(np)]
      iy = [[0.0 for j in range(np)] for i in range(np)]
      for i in range(np):
         for j in range(np):
            ix[i][j] = x+xx[i][j]
            iy[i][j] = y+yy[i][j]
      return ix, iy



   def text_position_at_point_in_box(selfnp, np, xx, yy):

      ix = [[0.0 for j in range(np)] for i in range(np)]
      iy = [[0.0 for j in range(np)] for i in range(np)]
      for i in range(np):
         for j in range(np):
            dx = 7.0/100.0
            if i<=1:
               ix[i][j] = xx[i][j]+dx
            else:
               ix[i][j] = xx[i][j]-dx
            if j<=1:
               iy[i][j] = yy[i][j]+dx
            else:
               iy[i][j] = yy[i][j]-dx
      return ix, iy
 
 
   
 # draw box
 #  1) flag1 : on/off
 #  2) flag2 : specific box
 #  2) flag3 : vproc
 #  2) flag4 : neibors
   def do_draw_in_box(self, ie, je, flag1, flag2, flag3, flag4):
 
      draw_line  = False
      draw_color = False
      if flag1==False or flag1[ie][je]==True:
         if flag2==False or flag2[ie][je]==True:
            draw_line  = True
            if flag3==False:
               draw_color = True
            else:
               if flag3[ie][je]==True:
                  draw_color = True
               else:
                  if flag4==False or flag4[ie][je]==True:
                     draw_color = True
                  else:
                     draw_color = False
      return draw_line, draw_color
 
 
 
   def draw_box(self, axis, color=False,  \
                        flag1=False, flag2=False, flag3=False, flag4=False):
 
      for ie in range(self.nx):
         for je in range(self.ny):
            draw_line, draw_color = self.do_draw_in_box(ie,je,flag1,flag2,flag3,flag4)
            if draw_line:
               if color!=False and draw_color:
                  patch = self.coord_to_box(ie,je,color[ie][je])
                  axis.add_patch(patch)
               patch = self.coord_to_box(ie,je,'none')
               axis.add_patch(patch)
 
 
 
   def write_values_in_box(self, axis, value, upper=False, iscolor=False, alpha=1.0,\
                          flag1=False, flag2=False, flag3=False, flag4=False):

      if upper:
         it, jt = self.tpos_up[0], self.tpos_up[1]
         tsize  = self.tsize_up
         if iscolor:
            tcolor='w'
         else:
            tcolor='k'
      else:
         it, jt = self.tpos[0], self.tpos[1]
         tsize  = self.tsize
         if iscolor:
            tcolor='k'
         else:
            tcolor='r'
  
      for ie in range(self.nx):
         for je in range(self.ny):
            draw_line, draw_color = self.do_draw_in_box(ie,je,flag1,flag2,flag3,flag4)
            if draw_line:
               axis.text(ie+it,je+jt,str(value[ie][je]),color=tcolor,size=tsize,alpha=alpha,\
                         va='center',ha='center')
   
 
 
   def draw_points_in_box(self, axis, bigger=False, var=False, kind=False, size=1.0, \
                                 flag1=False, flag2=False, flag3=False, flag4=False):

      np = self.np
      for ie in range(self.nx):
         for je in range(self.ny):
            draw_line, draw_color = self.do_draw_in_box(ie,je,flag1,flag2,flag3,flag4)
            if draw_line:
               px, py = self.point_position_in_box(ie,je,np,self.pts_x,self.pts_y)
               if not bigger: axis.scatter(px,py,s=13,marker='o',c='k')
               draw_point = isinstance(var,list)
               if bigger or draw_point:
                  tx, ty = self.text_position_at_point_in_box(np,px,py)
                  for i in range(np):
                     for j in range(np):
                        if isinstance(kind,list):
                           if kind[ie][je][i][j]==1:
                              color = 'b'
                              alpha = 1.0
                           elif kind[ie][je][i][j]==2:
                              color = 'r'
                              alpha = 1.0
                           elif kind[ie][je][i][j]==3:
                              color = 'm'
                              alpha = 1.0
                           else:
                              color = 'k'
                              alpha = 0.2
                        else:
                           color = 'k'
                           alpha = 1.0

                        if bigger:
                           #circle = plt.Circle((px[i][j],py[i][j]),radius=0.035,color=color)
                           circle = plt.Circle((px[i][j],py[i][j]),radius=0.027, \
                                                                  color=color,alpha=alpha)
                           axis.add_patch(circle)

                        if draw_point:
                           if isinstance(var[ie][je][i][j],float):
                              string = '{0:3.1f}'.format(var[ie][je][i][j])
                           else:
                              string = var[ie][je][i][j]
                           axis.text(tx[i][j],ty[i][j],string,color=color, \
                              size=self.tsize_pt*size,alpha=alpha,va='center',ha='center')
 
   
 
   def draw_points_in_box_tmp(self, axis, var, isu, flag1=False, flag2=False, \
                flag3=False, flag4=False, flag5=False, log=[False,False,False,False]):

      np = self.np
      for ie in range(self.nx):
        for je in range(self.ny):
          if flag1==False or flag1[ie][je]==True:
            if flag2==False or flag2[ie][je]==True:
              do_plot = False
              if flag3==False:
                do_plot = True
              else:
                if flag3[ie][je]==True:
                  do_plot = True
                else:
                  if flag4==False or flag4[ie][je]==True:
                    do_plot = False
              if do_plot:
                for i in range(self.np):
                  for j in range(self.np):
                    axis.scatter(ie+self.pts_x[i][j],je+self.pts_y[i][j],s=13,marker ='o', \
                                 c='k',linewidths=3,edgecolors='face')
                tsize = self.tsize/2.5
                dx = 7.0/100.0
                isdraw=False
                if log[0] and ie==0 and je==0:
                    j=np-1
                    i=np-1
                    isdraw=True
                    str = log[0]
                elif log[1] and ie==1 and je==0:
                    j=0
                    i=np-1
                    isdraw=True
                    str = log[1]
                elif log[2] and ie==0 and je==1:
                    j=np-1
                    i=0
                    isdraw=True
                    str = log[2]
                elif log[3] and ie==1 and je==1:
                    j=0
                    i=0
                    isdraw=True
                    str = log[3]
                if isdraw:
                    ix = ie+self.pts_x[i][j]
                    iy = je+self.pts_y[i][j]
                    if i<=1:
                      iy=iy+dx
                    else:
                      iy=iy-dx
                    if j<=1:
                      ix=ix+dx
                    else:
                      ix=ix-dx
                    if not flag5[ie][je]:
                      axis.text(ix,iy,str,color='k',size=tsize,va='center',ha='center')
                    elif flag5[ie][je]==0:
                      axis.text(ix,iy,str,color='b',size=tsize,va='center',ha='center')
                    elif flag5[ie][je]==1:
                      axis.text(ix,iy,str,color='r',size=tsize,va='center',ha='center')
  
