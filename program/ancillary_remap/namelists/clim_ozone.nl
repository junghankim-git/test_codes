&inputs

! remap files (fix)
isscrip        = .true.
remapfile      = '/home/jhkim/work/program/gen_remap_mat_with_latlon/remap_matrix/ne120rotated_ne030variable/remap_12.nc'

! input / output
infile         = '/data/KIM2.4/inputdata/ne120np4_rotated/clim_ozone.nc'
oufile         = './outputs/clim_ozone.nc'


! dimensions
ndims          = 2
dimnames       = 'nlev', 'ncol'

! variables
nvars          = 12
varnames       = 'ozone01', 'ozone02', 'ozone03', 'ozone04', 'ozone05', 'ozone06', 'ozone07', 'ozone08', 'ozone09', 'ozone10', 'ozone11', 'ozone12'
ndims_vars     = 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
dimnames_var01 = 'ncol', 'nlev'
dimnames_var02 = 'ncol', 'nlev'
dimnames_var03 = 'ncol', 'nlev'
dimnames_var04 = 'ncol', 'nlev'
dimnames_var05 = 'ncol', 'nlev'
dimnames_var06 = 'ncol', 'nlev'
dimnames_var07 = 'ncol', 'nlev'
dimnames_var08 = 'ncol', 'nlev'
dimnames_var09 = 'ncol', 'nlev'
dimnames_var10 = 'ncol', 'nlev'
dimnames_var11 = 'ncol', 'nlev'
dimnames_var12 = 'ncol', 'nlev'

/
