&inputs

! remap files (fix)
isscrip        = .true.
!remapfile      = '/home/jhkim/work/program/gen_remap_mat_with_latlon/remap_matrix/ne030rotated_ne030variable/remap_12.nc'
remapfile      = '/home/jhkim/work/program/gen_remap_mat_with_latlon/remap_12.nc'

! input / output
infile         = '/data/KIM3.0/inputdata/ne030np4_rotated/clim_aerosol.nc'
oufile         = './outputs/clim_aerosol.nc'


! dimensions
ndims          = 3
dimnames       = 'ncol', 'nlev', 'month'

! variables
nvars          = 6
varnames       = 'BC', 'OC', 'SO4', 'SEASALT', 'DUST', 'pressure'
ndims_vars     = 3, 3, 3, 3, 3, 1
dimnames_var01 = 'ncol', 'nlev', 'month'
dimnames_var02 = 'ncol', 'nlev', 'month'
dimnames_var03 = 'ncol', 'nlev', 'month'
dimnames_var04 = 'ncol', 'nlev', 'month'
dimnames_var05 = 'ncol', 'nlev', 'month'
dimnames_var06 = 'nlev'
dimnames_var07 = ''
dimnames_var08 = ''
dimnames_var09 = ''
dimnames_var10 = ''

/
