
#!/usr/bin/env python
import os
import sys
import netCDF4 as nc


'''
   call fail_message(err,'cannot open file : '//trim(remap_file))
   err = nf90_inq_dimid(fid,'src_grid_size',did)
   if (err.ne.nf90_noerr) isscrip = .false.
   err = nf90_inq_dimid(fid,'dst_grid_size',did)
   if (err.ne.nf90_noerr) isscrip = .false.
   err = nf90_inq_dimid(fid,'num_links',did)
   if (err.ne.nf90_noerr) isscrip = .false.
   err = nf90_inq_dimid(fid,'num_wgts',did)
'''


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



rfile1 = '/home/jhkim/work/program/scrip/3.remap_matrix/CS_WW3/ne030_rotated_to_WW3.nc'
rfile2 = '/data/KIM3.0/remap_matrix/ne030np4_rotated/bilinear/cs2ll_360x180_regular.nc'

print check_scrip(rfile1)
print check_scrip(rfile2)

