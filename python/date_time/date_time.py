#!/usr/bin/env python
import os
import sys
import numpy as np
import datetime
#from scanf import scanf



def set_datetime(yy,mm,dd,hh,tt,ss):
    '''
    # method 1
    date      = datetime.date(yy,mm,dd)
    time      = datetime.time(hh,tt,ss)
    date_time = datetime.datetime.combine(date,time)
    '''

    # method 2
    date_time = datetime.datetime.strptime( \
                '{0:04d}-{1:02d}-{2:02d} {3:02d}:{4:02d}:{5:02d}'.format(yy,mm,dd,hh,tt,ss), \
                '%Y-%m-%d %H:%M:%S')

    return date_time



def set_datetime_str(string):
    if len(string)>=8:
      yy = int(string[0:4])
      mm = int(string[4:6])
      dd = int(string[6:8])
      if len(string)==14:
          hh = int(string[8:10])
          tt = int(string[10:12])
          ss = int(string[12:14])
    else:
      print('check len(string) in set_datetime_str...')

    date_time = set_datetime(yy,mm,dd,hh,tt,ss)
    return date_time



def __main__():

    target_date = '20160626000000'
    in_path = '/s3/scratch/kiaps/kiaps-sys/jhkim/KIM3.0/h4dev3.0a_da/ENSfcst/'+target_date[0:10]
    ou_path = '/s4/home/kiaps/kiaps-sys/jhkim/Src/3.0/H4DEV/jhkim/ens2h4dev/exp'
    nsmpl = 50
    ntime = 7
    times = [-3,-2,-1,0,1,2,3]
    date_time = set_datetime_str(target_date)
    my_dt = []
    
    inlist = open('in_list.txt','w')
    oulist = open('ou_list.txt','w')
    for it in range(ntime):
        dt = datetime.timedelta(hours=times[it])
        my_dt.append(date_time+dt)
        #print(my_dt[it].strftime('%Y-%m-%d %H:%M:%S'))
        for ss in range(nsmpl):
            infile = '{0:}/{1:03d}/fcst-{2:}-{3:06d}.nc'.format(in_path,ss+1,my_dt[it].strftime('%Y%m%d%H%M%S'),it+3)
            oufile = '{0:}/ensSample_{1:}_{2:03d}.nc'.format(ou_path,my_dt[it].strftime('%Y%m%d%H'),ss+1)
            inlist.write(infile+'\n')
            oulist.write(oufile+'\n')
    inlist.close()
    oulist.close()
    


__main__()
