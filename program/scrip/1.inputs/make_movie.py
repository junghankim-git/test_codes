#!/bin/env /app/compilers/gcc/4.7.1/applib1/Python/2.7.3/bin/python
import os
import sys
import shutil


def cart_to_idx(nxx,nyy,ixx,iyy):
    return ixx+(iyy-1)*nxx
def idx_to_cart(nxx,nyy,idx):
    ixx = (idx-1)%nxx+1
    iyy = (idx-1)/nxx+1
    print('i = {}, x = {}, y = {}'.format(idx,ixx,iyy))
    return ixx, iyy

nx,ny = 240,193
#ys = [132,133,193]
ys = [191]


print(' ')
for iy in ys:
    print('* making gif: {}th row'.format(iy))
    filenames = ''
    for ix in range(nx):
        filenames = filenames+' figs/{0:06d}.png'.format(cart_to_idx(nx,ny,ix+1,iy))
    dirname = 'movies/{0:03d}'.format(iy)
    movie_file = './movie_{0:03d}.gif'.format(iy)

    print(' 1) making directory: {}'.format(dirname))
    if os.path.exists(dirname):
        #os.remove(dirname) # file
        #os.rmdir(dirname)  # direcotry
        shutil.rmtree(dirname)
    os.makedirs(dirname)

    print(' 2) copy files: {}th row'.format(iy))
    command = 'cp '+filenames+' '+dirname
    os.system(command)

    command = 'convert -delay 10 -loop 0 movies/{0:03d}/*.png {1:}'.format(iy,movie_file)
    print(' 3) makeing gif file: {}'.format(movie_file))
    os.system(command)

    print(' ')
