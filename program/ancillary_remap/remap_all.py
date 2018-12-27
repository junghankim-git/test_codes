#!/usr/bin/env python

import os
import sys


chkdir    = '/data/KIM2.4/inputdata/ne120np4_rotated'
basedir   = './namelists'
# 16 namelists
namelists = ['clim_aerosol.nl', 'landf.nl', 'maxsnowalb.nl', 'MODISems.nl', 'tg3.nl',          \
             'clim_omld.nl', 'gfrac.nl', 'lsmsk.nl', 'mingfrac.nl', 'mtnvar_srtm30.nl',        \
             'vtype.nl', 'clim_ozone.nl', 'hice.nl', 'maxgfrac.nl', 'MODISalb.nl', 'stype.nl']

for nl in namelists:
  print 'Processing namelist: ', nl
  #chkfile = chkdir+'/'+nl.replace('nl','nc')
  #print ' org file = ', chkfile
  #print os.path.isfile(chkfile)

  nlfile = basedir+'/'+nl
  print 'namelist file = ', nlfile
  if not os.path.isfile(nlfile):
    print 'check file : ', nlfile
    pass
  else:
    os.system('./remap.exe < '+nlfile)
  print ' '
