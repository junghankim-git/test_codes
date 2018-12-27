#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy
# load my class
SHARE_DIR='/home/jhkim/work/share/python'
sys.path.append(SHARE_DIR)
from class_draw_box import *
from class_sfc      import *
from utilities      import *


fig_dir = 'fig'
method_name = ['ep','up','lup','eup']



def get_flag(nx,ny):
   flag = numpy.full(shape=(nx,ny),fill_value=1).tolist()
   return flag



def get_proc_dist(nx,ny,nproc,mm=None):
   proc = numpy.full(shape=(nx,ny),fill_value=0,dtype=int).tolist()
   if nproc<=0 or nx*ny<nproc:
      print('error 1: check nproc in get_proc_dist ({})'.format(nproc))
      quit()
   if nx!=ny:
      print('error 1: check nx,ny in get_proc_dist ({},{})'.format(nx,ny))
      quit()

   sfc = SFC(nx)
   mesh = sfc.generate_sfc()

   nn   = [0 for i in range(nproc)]
   ista = [0 for i in range(nproc)]
   iend = [0 for i in range(nproc)]
   if mm==None:
      for ip in range(nproc):
         nn[ip],ista[ip],iend[ip] = decompose_1d(1,nx*ny,nproc,ip)
   else:
      ii = 0
      for ip in range(nproc):
         if ip==0:
            ista[ip] = 0
         else:
            ista[ip] = iend[ip-1]+1
         iend[ip] = ista[ip]+mm[ip]-1

   # set proc number
   ii = -1
   for ie in range(nx):
      for je in range(ny):
         for ip in range(nproc):
            if mesh[ie][je]>=ista[ip] and mesh[ie][je]<=iend[ip]:
               proc[ie][je] = ip
   patch = sfc.get_sfc_line()
   
   '''
   if nproc<=0:
      return proc
   else:
      return proc, patch
   '''
   return proc, patch




def get_kind_with_proc(np,nx,ny,method,proc):
   '''
   method<=0: only EP, method=1: set UP, method=2: set UP&LUP, method=3: set UP&LUP&EUP
   '''
   cnt = []
   if method<0:
      kind = None
   elif method==0:
      kind = 0
   else:
      # halo
      # H0 0  1  1  1  H1
      # H0 0  1  1  1  H1
      # H0 0  0  1  1  H1
      # H0 0  0  0  0  H0
      # H0 H0 H0 H0 H0 H0
      kind = numpy.empty(shape=(nx,ny,np,np)).tolist()
      for ii in range(nx):
         for jj in range(ny):
            is_ww = ii-1>=0
            is_ee = ii+1<=nx-1
            is_ss = jj-1>=0
            is_nn = jj+1<=ny-1
            is_ws = is_ww and is_ss
            is_es = is_ee and is_ss
            is_wn = is_ww and is_nn
            is_en = is_ee and is_nn
            p_ij = proc[ii][jj]
            p_ww = proc[ii-1][jj]   if is_ww else -1
            p_ee = proc[ii+1][jj]   if is_ee else -1
            p_ss = proc[ii][jj-1]   if is_ss else -1
            p_nn = proc[ii][jj+1]   if is_nn else -1
            p_ws = proc[ii-1][jj-1] if is_ws else -1
            p_es = proc[ii+1][jj-1] if is_es else -1
            p_wn = proc[ii-1][jj+1] if is_wn else -1
            p_en = proc[ii+1][jj+1] if is_en else -1
            d_ww = p_ij != p_ww
            d_ee = p_ij != p_ee
            d_ss = p_ij != p_ss
            d_nn = p_ij != p_nn
            d_ws = p_ij != p_ws
            d_es = p_ij != p_es
            d_wn = p_ij != p_wn
            d_en = p_ij != p_en
            for i in range(np):
               for j in range(np):
                  value = 1   # initial value is UP
                  if i == 0 and j == 0:
                     if is_ww or is_ss or is_ws:
                        value = 0
                        if method >= 2:
                           if d_ww and d_ss and d_ws:
                              value = 2
                           elif method >= 3 and (d_ww or d_ss or d_ws) and ii!=0 and jj!=0:
                              value = 3
                  elif i == 0 and j != 0:
                     if is_ww:
                        value = 0
                        if method >= 2:
                           if d_ww:
                              value = 2
                           elif method >= 3 and j==np-1 and (d_ww or d_nn or d_wn) and jj!=ny-1:
                              value = 3
                  elif j == 0 and i != 0:
                     if is_ss:
                        value = 0
                        if method >= 2:
                           if d_ss and ((i != np-1) or d_es):
                                 value = 2
                           elif method >= 3 and i==np-1 and (d_ee or d_ss or d_es) and ii!=nx-1:
                              value = 3
                  kind[ii][jj][i][j] = value
                  cnt.append(value)

   nup  = cnt.count(1)
   nlup = nup+cnt.count(2)
   neup = nlup+cnt.count(3)
   nep  = neup+cnt.count(0)
   print 'n  EP = ',nep,100.0
   print 'n  UP = ',nup,float(nup)/float(nep)*100.0 if nup!=0 else 0.0
   print 'n LUP = ',nlup,float(nlup)/float(nep)*100.0 if nlup!=0 else 0.0
   print 'n EUP = ',neup,float(neup)/float(nep)*100.0 if neup!=0 else 0.0
   return kind



def get_value(set_value,np,nx,ny,kind=None):
   value = numpy.full(shape=(nx,ny,np,np),fill_value='',dtype=str).tolist()
   if set_value>=0:
      if kind!=None:
         ii = 0
         for jp in range(nx):
            for ip in range(ny):
               for j in range(np):
                  for i in range(np):
                     if set_value==0:
                        ii = ii+1
                        value[ip][jp][i][j] = ii
                     elif set_value==1 and kind[ip][jp][i][j]<=1:
                        ii = ii+1
                        value[ip][jp][i][j] = ii
                     elif set_value==2 and kind[ip][jp][i][j]<=2:
                        ii = ii+1
                        value[ip][jp][i][j] = ii
       
      else:
        if set_value==1:
           value[0][0][3][3] = 'A'
           value[1][0][0][3] = 'B'
           value[0][1][3][0] = 'C'
           value[1][1][0][0] = 'D'
   return value


def do_plot(np,nx,ny,flag,proc,ival,kind,fname,patch=False):
   if np==3:
      gll = [-1.0,0.0,1.0]
   elif np==4:
      gll = [-1.0,-0.5,0.5,1.0]
   plot         = class_draw_box(np,gll,nx,ny)
   color        = plot.get_box_color(proc)
   pltfig, axis = plot.create_canvas('')
   plot.draw_box(axis,color=color,flag1=flag)
   if nx==2 and ny==2:
      num = [[1,3],[2,4]]
   elif nx==3 and ny==3:
      num = [[1,4,7],[2,5,8],[3,6,9]]
   plot.write_values_in_box(axis,num,upper=False,iscolor=True,alpha=0.3,flag1=flag,flag4=proc)
   #plot.draw_points_in_box(axis,bigger=True,var=ival,kind=kind,flag1=flag,flag4=proc)
   plot.draw_points_in_box(axis,bigger=True,var=ival,kind=kind,size=0.6,flag1=flag,flag4=proc)
   if patch!=False:
      axis.add_patch(patch)
   pltfig.savefig(fname)



def generate_3x3(method=0):

   np    = 3
   nx    = 3
   ny    = 3
   nproc = 1
   set_v = -1
   flag  = get_flag(nx,ny)
   proc, patch = get_proc_dist(nx,ny,nproc)
   kind  = get_kind_with_proc(np,nx,ny,method,proc)
   ival  = get_value(set_v,np,nx,ny)
   fname = fig_dir+'/'+'3x3.png'
   return np, nx, ny, flag, proc, ival, kind, fname



def generate_big(method=3):

   np    = 4
   nx    = 10
   ny    = 10
   nproc = 10
   set_v = -1
   flag  = get_flag(nx,ny)
   proc, patch = get_proc_dist(nx,ny,nproc)
   kind  = get_kind_with_proc(np,nx,ny,method,proc)
   ival  = get_value(set_v,np,nx,ny)
   fname = fig_dir+'/'+'big.png'
   return np, nx, ny, flag, proc, ival, kind, fname, patch



def generate_exp(method=0):

   np    = 4
   nx    = 2
   ny    = 2
   nproc = 2
   set_v = 1
   flag  = get_flag(nx,ny)
   proc, patch = get_proc_dist(nx,ny,nproc)
   kind  = get_kind_with_proc(np,nx,ny,method,proc)
   ival  = get_value(set_v,np,nx,ny)
   fname = fig_dir+'/'+'exp.png'
   return np, nx, ny, flag, proc, ival, kind, fname



# method: 0(EP), 1(UP), 2(LUP), 3(EUP)
#  shape: processes distribution (1: all rank0, 2: 2xrank0+2xrank1, 3: 3xrank0+1xrank1)
def generate_fig(np,method,shape):

   if method!=0 and method!=1 and method!=2 and method!=3:
      print('check method in generate_fig ({})...'.format(method))
      quit()
   if shape!=0 and shape!=1 and shape!=2 and shape!=3:
      print('check shape in generate_fig ({})...'.format(shape))
      quit()

   if shape==0:
      nproc  = 2
      mm     = None
      set_v  = -1
      exname = ''
   if shape==1:
      nproc  = 1
      mm     = None
      set_v  = 1
      exname = '_label'+str(shape)
   elif shape==2:
      nproc  = 2
      mm     = None
      set_v  = 1
      exname = '_label'+str(shape)
   elif shape==3:
      nproc  = 2
      mm     = [3,1]
      set_v  = 1
      exname = '_label'+str(shape)

   #np    = 4
   if np==3:
      nx = 3
      ny = 3
   elif np==4:
      nx = 2
      ny = 2
   flag  = get_flag(nx,ny)
   proc, patch = get_proc_dist(nx,ny,nproc,mm)
   kind  = get_kind_with_proc(np,nx,ny,method,proc)
   ival  = get_value(set_v,np,nx,ny)
   fname = fig_dir+'/'+'fig_'+method_name[method]+exname+'.png'
   return np, nx, ny, flag, proc, ival, kind, fname



def generate_poster(method=0):

   np    = 4
   nx    = 2
   ny    = 2
   nproc = 2
   set_v = method
   flag  = get_flag(nx,ny)
   proc, patch = get_proc_dist(nx,ny,nproc)
   kind  = get_kind_with_proc(np,nx,ny,method,proc)
   ival  = get_value(set_v,np,nx,ny,kind)
   #ival  = get_value(0,np,nx,ny,kind)
   fname = fig_dir+'/'+'poster_'+method_name[method]+'.png'
   return np, nx, ny, flag, proc, ival, kind, fname




def main():
   
   '''
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_3x3()
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_exp()
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)

   np, nx, ny, flag, rank, value, is_uniq, fname = generate_fig(4,0,0)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_fig(4,1,0)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_fig(4,2,0)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_fig(4,0,1)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_fig(4,2,1)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_fig(4,2,2)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_fig(4,2,3)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_fig(4,3,1)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_fig(4,3,2)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_fig(4,3,3)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)

   np, nx, ny, flag, rank, value, is_uniq, fname, patch = generate_big()
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname,patch)

   np, nx, ny, flag, rank, value, is_uniq, fname = generate_poster(0)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_poster(1)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_poster(2)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   '''
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_fig(4,2,0)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   np, nx, ny, flag, rank, value, is_uniq, fname = generate_fig(3,2,0)
   do_plot(np,nx,ny,flag,rank,value,is_uniq,fname)
   
   plt.show()



if __name__ == '__main__':
   main()


