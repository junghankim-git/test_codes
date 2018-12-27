!-------------------------------------------------------------------------------
   program main_program
!-------------------------------------------------------------------------------
!
!  abstract :  (main program) NetCDF file compressor utility
!
!  history log :
!    2017-02-28  junghan kim    initial setup
!
!  structure :
!
!-------------------------------------------------------------------------------
   use kinds,      only: i4, r4, r8, l4
   use opt_parser, only: opt_args_t, max_len_str_, opt_bool_, opt_store_,      &
                         opt_logical_, opt_int_, opt_real_, opt_real8_, opt_string_, &
                         opt_initialize, opt_finalize, opt_add, opt_print,     &
                         opt_get, opt_get_arg, opt_show_help,                  &
                         max_len_arrs_, string2array
   use netcdf_utility, only: len_filenam, len_att_nam, len_var_nam,            &
                             ncu_t, ncu_initialize, ncu_finalize, ncu_run
!
   type(opt_args_t) :: opts
   type(ncu_t)      :: nc
   integer(i4)      :: compr_lev, iconvert
   character(len_filenam) :: infilename, oufilename
   character(max_len_str_) :: varnames_opt
   character(max_len_arrs_), dimension(:), allocatable :: varnames
!
   call opt_initialize(opts)
   call opt_add(opts,'-o',   '--output', opt_store_,'out.nc','output file name')
   call opt_add(opts,'-m',   '--convert', opt_store_,1,'0(nc4->nc4), 1(pnc->nc4), 2(nc4->pnc)')
   call opt_add(opts,'-c',    '--level', opt_store_,       0,  'compress level')
   call opt_add(opts,'-v','--variables', opt_store_,   'all','output variables')
   if (opts%nargs.le.0) then
     print *, 'check number of arguments of program....'
     call opt_show_help(opts)
   endif
   call opt_get_arg(opts,1,infilename)
   call opt_get(opts,'-o',oufilename)
   call opt_get(opts,'-m',iconvert)
   call opt_get(opts,'-c',compr_lev)
   call opt_get(opts,'-v',varnames_opt)
   call string2array(varnames_opt,varnames)
   call opt_finalize(opts)
!
   call ncu_initialize(nc,trim(infilename),trim(oufilename),compr_lev,iconvert,varnames)
   call ncu_run(nc)
   call ncu_finalize(nc)
!
   deallocate(varnames)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine dummy
   implicit none
!
!
   end subroutine dummy
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program main_program
!-------------------------------------------------------------------------------
