!-------------------------------------------------------------------------------
   program main_program
   use kinds,      only: i4, r4, r8, l4
   use opt_parser, only: opt_args_t, max_len_str_, opt_bool_, opt_store_,      &
                         opt_logical_, opt_int_, opt_real_, opt_real8_, opt_string_, &
                         opt_initialize, opt_finalize, opt_add, opt_print,     &
                         opt_get, opt_get_arg, opt_show_help
   use nfc_file_compressor, only: nfc_t, nfc_initialize, nfc_finalize, nfc_compress
!
   type(opt_args_t) :: opts
   type(nfc_t)      :: nc
   integer(i4)      :: compr_lev
   character(len=max_len_str_) :: infilename, oufilename
!
   call opt_initialize(opts)
   call opt_add(opts,'-c', '--level', opt_store_,       1,  'compress level')
   call opt_add(opts,'-o','--output', opt_store_,'out.nc','output file name')
   if (opts%nargs.le.0) then
     print *, 'check number of arguments of program....'
     call opt_show_help(opts)
   endif
   call opt_get_arg(opts,1,infilename)
   call opt_get(opts,'-c',compr_lev)
   call opt_get(opts,'-o',oufilename)
   call opt_finalize(opts)
!
   call nfc_initialize(nc,trim(infilename),trim(oufilename),compr_lev) 
!
   call nfc_compress(nc)
!
   call nfc_finalize(nc)
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
