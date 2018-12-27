&inputs

! remap files (fix)
isscrip        = .true.
remapfile      = '/home/jhkim/work/program/gen_remap_mat_with_latlon/remap_matrix/ne120rotated_ne030variable/remap_12.nc'

! input / output
infile         = '/data/KIM2.4/inputdata/ne120np4_rotated/MODISalb.nc'
oufile         = './outputs/MODISalb.nc'


! dimensions
ndims          = 2
dimnames       = 'time', 'ncol'

! variables
nvars          = 6
varnames       = 'vis1', 'vis2', 'vis3', 'nir1', 'nir2', 'nir3'
ndims_vars     = 2, 2, 2, 2, 2, 2
dimnames_var01 = 'ncol', 'time'
dimnames_var02 = 'ncol', 'time'
dimnames_var03 = 'ncol', 'time'
dimnames_var04 = 'ncol', 'time'
dimnames_var05 = 'ncol', 'time'
dimnames_var06 = 'ncol', 'time'
dimnames_var07 = ''
dimnames_var08 = ''
dimnames_var09 = ''
dimnames_var10 = ''

/
