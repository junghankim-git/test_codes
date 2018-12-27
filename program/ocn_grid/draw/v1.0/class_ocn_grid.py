import sys
sys.path.append('/home/jhkim/work/share/python')
from utilities       import decompose_1d
from class_plot_cart import *



class ocn_grid:

    def __init__(self):

        self.xproc    = 3
        self.nprocs   = 9
        self.nx       = 15
        self.ny       = 12
        self.gsize    = self.nx*self.ny

        nprocs        = self.nprocs
        nx            = self.nx
        ny            = self.ny

        self.figs_dir = './figs'
        self.cart     = plot_cart(nx,ny,nprocs,vproc=4)

        
        # cart index,1d index,mask
        self.idx_2d = [[[i+1,j+1] for j in range(ny)] for i in range(nx)]
        self.idx_1d = [[ j*nx+i+1 for j in range(ny)] for i in range(nx)]

        # cart index,1d index,mask
        self.mask, self.mask_l, self.mask_s = self.gen_lsmask(nx,ny)
        self.ranks                          = self.domain_decompose(nx,ny,nprocs,self.mask)



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





    def domain_decompose(self,nx,ny,nprocs,mask):

        xproc = self.xproc
        if nprocs%xproc!=0:
            print('error in nprocs ({})'.format(nprocs))
            quit()
        ijpr = nprocs
        ipr  = xproc
        jpr  = nprocs/xproc
        print('* number of processes')
        print('total nprocs = {}, x nprocs = {}, y nprocs = {}'.format(ijpr,ipr,jpr))
        print(' ')

        mnproc = [[-1,-1] for i in range(nprocs)]
        for ip in range(nprocs):
            mproc = (ip+1)%ipr
            if mproc==0: mproc = xproc
            nproc = ip/ipr+1 
            mnproc[ip] = [mproc-1,nproc-1]
            print(' mnproc[{}] = {}'.format(ip,mnproc[ip]))
        print(' ')

        ii_pe = [[-1,-1] for i in range(nprocs)]
        jj_pe = [[-1,-1] for i in range(nprocs)]
        print('* process grid ')
        for ip in range(nprocs):
            ii_pe[ip][1], ii_pe[ip][0], iend = decompose_1d(1,nx,ipr,mnproc[ip][0])
            jj_pe[ip][1], jj_pe[ip][0], iend = decompose_1d(1,ny,jpr,mnproc[ip][1])
            ii_pe[ip][0] = ii_pe[ip][0]-1
            jj_pe[ip][0] = jj_pe[ip][0]-1
            print('{:4d}, {:4d}, {:4d}, {:4d}'.format(ii_pe[ip][0],ii_pe[ip][1],jj_pe[ip][0],jj_pe[ip][1]))
        print(' ')

        print('* domain decomposition')
        ranks  = [[-1 for j in range(ny)] for i in range(nx)]
        tmp    = [0 for i in range(nx)]
        for ip in range(nprocs):
            for i in range(ii_pe[ip][0],ii_pe[ip][0]+ii_pe[ip][1]):
                for j in range(jj_pe[ip][0],jj_pe[ip][0]+jj_pe[ip][1]):
                    ranks[i][j] = ip
        for j in range(ny-1,-1,-1):
            for i in range(nx):
                tmp[i] = ranks[i][j]
            print(' ranks[:][{:3d}] = {}'.format(j,tmp[:]))
        print(' ')
   
        return ranks




    def draw_all(self):
        #pass
        self.draw_grid_structure('01_grid')
        self.draw_entire_process('02_entire_map')
        self.draw_specific_process('03_specific_map')



    def draw_grid_structure(self,title,isvproc=False):
        pltfig,axis = self.cart.create_canvas('')
        iscolor = False
        self.cart.grid_structure(axis,self.mask,self.ranks,iscolor=iscolor,isvproc=isvproc)
        self.cart.write_index_in_box(axis,self.mask,self.ranks,self.idx_2d, \
                                     isup=True,iscolor=iscolor,isvproc=isvproc)
        self.cart.write_index_in_box(axis,self.mask,self.ranks,self.idx_1d, \
                                     isup=False,iscolor=iscolor,isvproc=isvproc)
        pltfig.savefig(self.figs_dir+'/'+title)
        del axis



    def draw_entire_process(self,title,isvproc=False):
        pltfig,axis = self.cart.create_canvas('')
        iscolor = True
        self.cart.grid_structure(axis,self.mask,self.ranks,iscolor=iscolor,isvproc=isvproc)
        self.cart.write_index_in_box(axis,self.mask,self.ranks,self.idx_2d, \
                                     isup=False,iscolor=iscolor,isvproc=isvproc)
        pltfig.savefig(self.figs_dir+'/'+title)
        del axis



    def draw_specific_process(self,title,isvproc=True):
        pltfig,axis = self.cart.create_canvas('')
        iscolor = True
        self.cart.grid_structure(axis,self.mask,self.ranks,iscolor=True,isvproc=isvproc)
        self.cart.write_index_in_box(axis,self.mask,self.ranks,self.idx_1d, \
                                     isup=False,iscolor=iscolor,isvproc=isvproc)
        self.cart.write_index_in_box(axis,self.mask,self.ranks,self.idx_2d, \
                                     isup= True,iscolor=iscolor,isvproc=isvproc)
        pltfig.savefig(self.figs_dir+'/'+title)
        del axis



    def show(self):
        self.cart.show()


