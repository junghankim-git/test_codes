&inputs

 ! scrip input with lat-lon infomations
 np1       = 4
 ne1       = 30
 nprocs1   = 1
 rotated1  = .true.
 filename1 = '/scratch/jhkim/TestBed/KIM/Output/3.0/ne30/gnu/10h/SW_G/101/UP-20110725220000-000010.nc'
 lonname1  = 'lons'
 latname1  = 'lats'

 np2       = 4
 ne2       = 30
 nprocs2   = 1
 rotated2  = .false.
 filename2 = './choi/UP-000001.nc'
 lonname2  = 'lon'
 latname2  = 'lat'


 ! generate remap matrix
 dogenremap = .true.


 ! scrip file
 scripfile1 = './scrip_1.nc'
 scripfile2 = './scrip_2.nc'

 ! remap matrix filename
 remapfile12 = './remap_12.nc'
 remapfile21 = './remap_21.nc'

/
