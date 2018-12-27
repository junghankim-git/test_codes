import numpy as np
import matplotlib.pyplot as plt


base_size  = 4.0



def decompose_1d(n1,n2,nprocs,rank):
    domain = n2-n1+1
    res    = domain/nprocs
    extra  = domain%nprocs
    ista   = res*rank+n1+min(rank,extra)
    iend   = ista+res-1
    if rank<=extra-1:
        iend = iend+1
    res = iend-ista+1
    return res,ista,iend



def speedup_func(nprocs, serial):
    speed = 1.0/(serial+(1.0+serial)/nprocs)
    return speed



#def fit_func(func, xxx, yyy, n=10, xmin=min(xxx), xmax=max(xxx)):
#    popt, pcov = curve_fit(func, xxx, yyy)



def tolist_nd(var,nd):
    ndim = get_ndims(var)
    if ndim>nd:
        print('check dimension of var in tolist_nd..')
        quit()
    for i in range(nd-ndim):
       var = [var]
    return var



def get_ndims(var):
    ndim = 0
    if np.isscalar(var):
        return ndim
    if isinstance(var,list):
        ndim = 1
    else:
        return ndim
    if isinstance(var[0],list):
        ndim = 2
    else:
        return ndim
    if isinstance(var[0][0],list):
        ndim = 3
    else:
        return ndim
    if isinstance(var[0][0][0],list):
        ndim = 4
    else:
        return ndim
    if isinstance(var[0][0][0][0],list):
        ndim = 5
    else:
        return ndim
    print('dimension of var is greater than 5 or not list...')
    quit()




def get_dimsize(var):
    ndims = get_ndims(var)
    if ndims == 0:
        print 'error: variable is scalar..'
        quit()
    elif ndims == 1:
        ncol = len(var[:])
    elif ndims == 2:
        ncol = len(var[0][:])
    elif ndims == 3:
        ncol = len(var[0][0][:])
    elif ndims == 4:
        ncol = len(var[0][0][0][:])
    elif ndims == 5:
        ncol = len(var[0][0][0][0][:])
    return ncol



def len_file(filename):
    lines = 0
    for line in open(filename):
        lines = lines + 1
    return lines



def rotate_yaxis(axis):
    y_labels = axis.get_yticklabels()
    for label in y_labels:
        label.set_rotation(90)



def make_pltfig(ncols=1,nrows=1,base=base_size,mrg_wgt=1.0,yscale=0.75,title='',title_size=None):
    if title_size==None: title_size = 4.5*base
    if ncols==2 and nrows==3:
        mrg_wgt = 1.3
    lsize = base/7.0*mrg_wgt
    rsize = base/7.0*mrg_wgt
    dsize = base/7.0*mrg_wgt
    #tsize = base/8.0*mrg_wgt
    tsize = base/5.5*mrg_wgt
    wsize = base/2.5*mrg_wgt
    hsize = base/2.0*mrg_wgt
    unit_xsize = base
    unit_ysize = base*yscale 
    xsize = ncols*unit_xsize+lsize+rsize+(ncols-1)*wsize
    ysize = nrows*unit_ysize+dsize+tsize+(nrows-1)*hsize
    left,right,bottom,top = lsize/xsize,1.0-rsize/xsize,dsize/ysize,1.0-tsize/ysize
    wspace,hspace = wsize/xsize,hsize/ysize
    pltfig = plt.figure(figsize=(xsize,ysize))
    pltfig.subplots_adjust(left=left,right=right,bottom=bottom,top=top,wspace=wspace,hspace=hspace)
    pltfig.suptitle(title,fontsize=title_size)
    return pltfig



def create_plot(ncols=1,nrows=1,yscale=0.72,title='',base=base_size,title_size=None):

    if title_size==None: title_size = 4.5*base
    pltfig = make_pltfig(ncols,nrows,base,1.0,yscale,title,title_size)

    if nrows==1 and ncols==1:
        axis = pltfig.add_subplot(1,1,1)
        rotate_yaxis(axis)
    else:
        axis = [[pltfig.add_subplot(nrows,ncols,j*ncols+i+1) \
                 for j in range(nrows)] for i in range(ncols)]
        for i in range(ncols):
            for j in range(nrows):
                rotate_yaxis(axis[i][j])
    return pltfig,axis



def axis_init(axis,title,xlabel,ylabel,xmin,xmax,ymin,ymax,base=base_size,title_size=None,\
                                        xlabel_size=None,ylabel_size=None,tick_size=None):
    if  title_size==None:  title_size = 3.5*base
    if xlabel_size==None: xlabel_size = 3.0*base
    if ylabel_size==None: ylabel_size = 3.0*base
    if   tick_size==None:   tick_size = 2.5*base
    axis.set_title(title,fontsize=title_size)
    axis.set_xlabel(xlabel,fontsize=xlabel_size)
    axis.set_ylabel(ylabel,fontsize=ylabel_size)
    axis.tick_params(labelsize=tick_size)
    axis.set_xlim(xmin,xmax)
    axis.set_ylim(ymin,ymax)
    axis.grid(True)



def draw_plots(axis, title, carts, comps, xxx, yyy, fmts=None, \
               xlabel=None, ylabel=None, xmin=None, xmax=None, ymin=None, ymax=None, leg_loc=1, base=base_size):
    nmax_carts = 5
    nmax_comps = 4
    carts = tolist_nd(carts,1)   # ['cart1', 'cart2']
    comps = tolist_nd(comps,2)   # [['comp1-1','comp1-2'],['comp2-1','comp2-2','comp2-3']]
    xxx   = tolist_nd(xxx,2)     # [[1,2,3,4,5,....,9,10],[21,22,23,24,25,.....,29,30,31]]
    yyy   = tolist_nd(yyy,3)     # [[[y1,y2,...],[y1,y2,...]],[[y1,y2,...],[y1,y2,...],[y1,y2,...]]]
    ncarts = len(carts)
    if fmts!=None: fmts  = tolist_nd(fmts,2)
    if len(comps)!=ncarts:
        print('check dimensions of carts and comps')
        quit()
    if len(xxx)!=ncarts or len(yyy)!=ncarts:
        print('check dimensions of xxx({}) and yyy({})'.format(len(xxx),len(yyy)))
        quit()
    for i in range(ncarts):
        ncomps = len(comps[i])
        if len(yyy[i])!=ncomps:
            print('check {}th dimensions of yyy({})'.format(i,len(yyy[i])))
            quit()
        for j in range(ncomps):
            if len(xxx[i])!=len(yyy[i][j]):
                print('check ({},{})th dimensions of xxx({}) and yyy({})'.\
                                        format(i,j,len(xxx[i]),len(yyy[i][j])))
    vxmin = np.amin(xxx); vxmax = np.amax(xxx)
    vymin = np.amin(yyy); vymax = np.amax(yyy)
    xhalo = vxmax*0.05
    if xmin==None: xmin = vxmin-xhalo
    if xmax==None: xmax = vxmax+xhalo

    yhalo = vymax*0.05
    if ymin==None:
      if vymin>=0.0 and vymin<=vymax*0.1:
          ymin = 0.0
      else:
          ymin = vymin-yhalo
    if ymax==None: ymax = vymax+yhalo
    fmt_cart = ['k','r','b','g','c']
    fmt_comp = ['o-','^--','*:','.-','--']
    if xlabel==None: xlabel = 'x'
    if ylabel==None: ylabel = 'y'
    axis_init(axis,title,xlabel,ylabel,xmin,xmax,ymin,ymax)
    for i in range(ncarts):
       for j in range(len(comps[i])):
           if carts[i]=='':
               label = comps[i][j]
           else:
               label = carts[i]+': '+comps[i][j]
           if fmts==None:
               fmtt = fmt_cart[i]+fmt_comp[j]
           else:
               fmtt = fmts[i][j]
           axis.plot(xxx[i],yyy[i][j],fmtt,label=label)
    y_labels = axis.get_yticklabels()
    for label in y_labels:
        label.set_rotation(90)
    axis.grid(True)
    axis.legend(loc=leg_loc,fontsize=2.8*base,shadow=True)



def add_plots(axis, comps, xxx, yyy, fmts=None, leg_loc=1, base=base_size):
    nmax_carts = 5
    nmax_comps = 4
    comps  = tolist_nd(comps,1)   # ['comp1-1','comp1-2']
    xxx    = tolist_nd(xxx,1)     # [1,2,3,4,5,....,9,10]
    yyy    = tolist_nd(yyy,2)     # [[y1,y2,...],[y1,y2,...]]
    ncomps = len(comps)
    if len(yyy)!=ncomps:
        print('check dimensions of yyy'.format(len(yyy)))
        quit()
    for j in range(ncomps):
        if len(xxx)!=len(yyy[j]):
            print('check {}th dimensions of xxx({}) and yyy({})'.format(j,len(xxx),len(yyy[j])))
    if fmts==None:
        fmts = ['ko-','k^--','k*:']
    else:
        if len(fmts)!=ncomps:
            print('check dimensions of fmts'.format(ncomps))
            quit()

    for j in range(ncomps):
        label = comps[j]
        axis.plot(xxx,yyy[j],fmts[j],label=label)
    axis.legend(loc=leg_loc,fontsize=2.8*base,shadow=True)



def get_kim_profile(filename,tags):
    ndims = get_ndims(tags)
    if ndims>1:
        print('check dimension of tags in get_kim_profile...')
        quit()
    tags = tolist_nd(tags,1)
    ntags = len(tags)
    values = [0 for i in range(ntags)]

    file = open(filename,'r')
    do_find = False
    for line in file:
       slines = line.split()
       if len(slines)<4: continue
       if slines[0]=='name' and slines[1]=='ncalls': do_find = True
       if do_find:
           for i in range(ntags):
               if slines[0]==tags[i]:
                   values[i] = float(slines[3])
    file.close()

    return values





