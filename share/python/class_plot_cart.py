import os
import sys
import time
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches



class plot_cart:

    def __init__(self,nx,ny,nprocs,vproc=1):
        self.nprocs  = nprocs
        self.nx      = nx
        self.ny      = ny

        print('* processes and domain')
        print(' - nprocs = '.format(nprocs))
        print(' - nx     = '.format(nx))
        print(' - ny     = '.format(ny))
        print(' ')

        # text size
        self.canv_size = 1.0
        if self.canv_size == 1.0:
            self.tsize  = 140./float(nx)
            self.utsize = 140./float(nx)*0.8
        else:
            self.tsize  = 160./float(nx)
            self.utsize = 160./float(nx)*0.8

        # view processor
        self.vproc = vproc

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
                            patch = self.coord2box(i,j,'none')
                            axis.add_patch(patch)
                    else:
                        patch = self.coord2box(i,j,color)
                        axis.add_patch(patch)




    def write_index_in_box(self,axis,mask,rank,index,isup=False,iscolor=False,isvproc=False):
        if isup:
            it,jt = self.itu,self.jtu
            tsize = self.utsize
            if iscolor:
                #tcolor='w'
                tcolor='k'
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


    def show(self):
        plt.show()




