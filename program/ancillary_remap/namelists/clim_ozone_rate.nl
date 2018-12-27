&inputs

! remap files (fix)
isscrip        = .true.
remapfile      = '/home/jhkim/work/program/gen_remap_mat_with_latlon/remap_matrix/ne120rotated_ne030variable/remap_12.nc'

! input / output
infile         = '/data/KIM2.4/inputdata/ne120np4_rotated/clim_ozone_rate.nc'
oufile         = './outputs/clim_ozone_rate.nc'


! dimensions
ndims          = 4
dimnames       = 'nlat', 'nlev', 'nday', 'nmonth'

! variables
nvars          = 4
varnames       = 'lats', 'pressure', 'production_rate', 'destruction_rate'
ndims_vars     = 1, 3, 4, 4
dimnames_var01 = 'nlat'
dimnames_var02 = nlev, nday, nmonth
dimnames_var03 = nlat, nlev, nday, nmonth
dimnames_var04 = nlat, nlev, nday, nmonth
dimnames_var05 = ''
dimnames_var06 = ''
dimnames_var07 = ''
dimnames_var08 = ''
dimnames_var09 = ''
dimnames_var10 = ''

/
