!-------------------------------------------------------------------------------
   program main
   use kinds,      only: i4, r4, r8, l4
   use opt_parser, only: opt_args_t, max_len_str_, opt_bool_, opt_store_,      &
                         opt_logical_, opt_int_, opt_real_, opt_real8_, opt_string_, &
                         opt_initialize, opt_finalize, opt_add, opt_print,     &
                         opt_get, opt_get_arg, opt_show_help
   use nc_control, only: netcdf_t,                                             &
                         nc_initialize, nc_finalize, nc_define, nc_compress
!
   type(opt_args_t)   :: opts
   type(netcdf_t)     :: nc
   integer(i4)        :: compr_lev
   character(len=max_len_str_) :: infilename, oufilename
!
   call opt_initialize(opts)
   call opt_add(opts,'-c', '--level',       1,  'compress level')
   call opt_add(opts,'-o','--output','out.nc','output file name')
   if (opts%nargs.le.0) then
     print *, 'check number of arguments of program....'
     call opt_show_help(opts)
   endif
   call opt_get_arg(opts,1,infilename)
   call opt_get(opts,'-c',compr_lev)
   call opt_get(opts,'-o',oufilename)
   call opt_finalize(opts)
!
   call nc_initialize(nc,trim(infilename),trim(oufilename),compr_lev) 
!
   call nc_define(nc)
   call nc_compress(nc)
!
   call nc_finalize(nc)
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
   end program main
!-------------------------------------------------------------------------------
