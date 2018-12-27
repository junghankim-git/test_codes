#!/usr/bin/env python
import sys
sys.path.append('/home/jhkim/work/share/python')
from class_draw_remap import *


def read_namelist(fname):

    rfile      = ''
    sfile      = ''
    dfile      = ''
    nvars      = 0
    varnames   = []
    ndims_vars = []

    nfile = open(fname,'r')
    i = 0
    for line in nfile:
        i = i+1
        if  line=='': continue
        sline = line.split()
        if len(sline)==0: continue
        if sline[0]=='!': continue
        if sline[0]=='remapfile':
              rfile = sline[2].replace("'","")
        if sline[0]=='infile':
              sfile = sline[2].replace("'","")
        if sline[0]=='oufile':
              dfile = sline[2].replace("'","")
        if sline[0]=='nvars':
              nvars = int(sline[2])
        if sline[0]=='varnames':
              for i in range(nvars):
                  varnames.append((sline[2+i].replace(',','')).replace("'",""))
        if sline[0]=='ndims_vars':
              for i in range(nvars):
                  ndims_vars.append(int(sline[2+i].replace(',','')))
 
    nfile.close()
    return rfile, sfile, dfile, varnames, ndims_vars




def check_scrip(filename):
    infile = nc.Dataset(filename, mode='r')
    isscrip = True

    try:
        nsize = len(infile.dimensions['src_grid_size'])
    except:
        isscrip = False
    try:
        nsize = len(infile.dimensions['dst_grid_size'])
    except:
        isscrip = False
    try:
        nsize = len(infile.dimensions['num_links'])
    except:
        isscrip = False
    try:
        nsize = len(infile.dimensions['num_wgts'])
    except:
        isscrip = False

    infile.close()
    return isscrip




def get_grid_info(gsize):
    '''
    !    ne030:   48602
    !    ne060:  194402
    !    ne120:  777602
    !    ne240: 3110402
    !    ne320: 
    !  360x180:   64800
    !  720x330:  237600 (scrip)
    !  720x360:  259200
    ! 1024x768:  786432
    '''
    info = 'none'
    if   gsize==48602:
        info = 'CS(ne030)'
    elif gsize==194402:
        info = 'CS(ne060)'
    elif gsize==777602:
        info = 'CS(ne120)'
    elif gsize==3110402:
        info = 'CS(ne240)'
    elif gsize==64800:
        info = 'LL(360x180)'
    elif gsize==237600:
        info = 'LL(720x330)'
    elif gsize==259200:
        info = 'LL(720x360)'
    elif gsize==786432:
        info = 'LL(1024x768)'

    elif gsize==46080:
        info = 'HYCOM(240x192)'
    elif gsize==46320:
        info = 'HYCOM(240x193)'

    return info




def main():

    remapfile, srcfile, dstfile, varnames, ndims = read_namelist('inputs_remap.nl')
    isscrip = check_scrip(remapfile)
    print('isscrip   = ',isscrip)
    print('remapfile = ',remapfile)
    print('srcfile   = ',srcfile)
    print('dstfile   = ',dstfile)
    #print(nvars)
    print('varnames  = ',varnames)
    print('ndims     = ',ndims)
    
    src_info_file = './latlon.nc'
    src_latname   = 'src_lats'
    src_lonname   = 'src_lons'
    dst_info_file = './latlon.nc'
    dst_latname   = 'dst_lats'
    dst_lonname   = 'dst_lons'
    
    
    remap = class_draw_remap(src_info_file, src_latname, src_lonname, \
                             dst_info_file, dst_latname, dst_lonname)
    remap.set_variables(srcfile, dstfile, varnames, ndims)
    s_grid = get_grid_info(remap.src_nsize)
    d_grid = get_grid_info(remap.dst_nsize)
    title = 'remap: '+s_grid+' to '+d_grid
    remap.draw(title, 0)



if __name__ == '__main__':
    main()



