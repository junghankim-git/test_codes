&inputs

! remap files (fix)
!remapfile      = '/home/jhkim/work/program/scrip/3.remap_matrix/CS_CS/ne030_rotated_to_ne030_unrotated.nc'
!remapfile      = '../3.remap_matrix/CS_WW3/ne030_rotated_to_WW3.nc'
!remapfile      = '/home/jhkim/work/program/scrip/3.remap_matrix/CS_HYCOM/ne030_rotated_to_hycom_wo_pole.nc'
!remapfile      = '/home/jhkim/work/program/scrip/3.remap_matrix/CS_HYCOM/ne030_rotated_to_hycom.nc'
remapfile      = '/home/khkim/test_codes/bilinear/make_matrix_HYCOM/cs2hycom_ne030np4_rotated_0240x0193_bilinear.nc'

! input / output
infile         = '/scratch/jhkim/TestBed/KIM/Output/3.0/ne30/gnu/10h/SW_G/101/UP-20110725220000-000010.nc'
oufile         = './result.nc'


! dimensions
ndims          = 1
dimnames       = 'ncol'

! variables
nvars          = 2
varnames       = 'ps', 'u10m'
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
