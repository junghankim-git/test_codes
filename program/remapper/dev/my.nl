&inputs


! remap files (fix)
!
!remapfile = '/home/jhkim/work/program/scrip/3.remap_matrix/CS_WW3/ne030_rotated_to_WW3.nc'
!remapfile = '/scratch/jhkim/TestBed/KIM/Input/remap_matrix/ne030np4_rotated/bilinear/cs2ll_1024x768_regular.nc'
!remapfile = '/scratch/jhkim/TestBed/KIM/Input/remap_matrix/ne030np4_rotated/bilinear/cs2ll_720x360_regular.nc'
remapfile = '/scratch/jhkim/TestBed/KIM/Input/remap_matrix/ne030np4_rotated/bilinear/cs2ll_360x180_regular.nc'
!remapfile = '/scratch/jhkim/TestBed/KIM/Input/remap_matrix/ne030np4_rotated/bilinear/ll2cs_1024x768_regular.nc'

! input / output
infile      = '/scratch/jhkim/TestBed/Data/3.0.07.02/10h/SW_G/101/UP-20110725140000-000002.nc'
!infile      = '/home/jhkim/TestBed/KIM/3.0.micros/3.0.02/11.dev/Exp_10h/SW_G/wav-file007200_charn.nc'
!infile      = '/home/jhkim/TestBed/KIM/3.0.micros/3.0.05/02.dev/Exp_10h/SW_G/wav_file_send_007200.nc'
!infile      = '/home/icna/KIM/3.0.05/Exp_10h/SW_G/wav_file_recv_007200.nc'

oufile      = './test.nc'


! dimensions
ndims          = 1
dimnames       = 'ncol'

! variables
nvars          = 2
varnames       = 'u10m', 'v10m'
ndims_vars     = 1, 1
dimnames_var01 = 'ncol'
dimnames_var02 = 'ncol'
dimnames_var03 = ''
dimnames_var04 = ''
dimnames_var05 = ''
dimnames_var06 = ''
dimnames_var07 = ''
dimnames_var08 = ''
dimnames_var09 = ''
dimnames_var10 = ''

/
