!-------------------------------------------------------------------------------
   program main
   use kinds,      only: i4, r4, r8, l4
   use opt_parser, only: opt_args_t, max_len_str_, opt_bool_, opt_store_, opt_append_, &
                         opt_logical_, opt_int_, opt_real_, opt_real8_, opt_string_, &
                         opt_initialize, opt_finalize, opt_add, opt_print,     &
                         opt_get, opt_get_arg, opt_show_help
!
   type(opt_args_t) :: opts
   logical(l4)      :: val_l4
   integer(i4)      :: val_i4
   real(r4)         :: val_r4
   real(r8)         :: val_r8
   character(len=max_len_str_) :: val_str, arg
   integer(i4) :: i
!
   call opt_initialize(opts)
   call opt_add(opts,'-l', '--logical', opt_bool_ , .false., 'logical precision')
   call opt_add(opts,'-i', '--integer', opt_store_,       0, 'integer precision')
   call opt_add(opts,'-r', '--real'   , opt_store_,     0.0,  'single precision')
   call opt_add(opts,'-d', '--double' , opt_store_,  0.0_r8,  'double precision')
   call opt_add(opts,'-s', '--string' , opt_store_,  'none',            'string')
!
   call opt_get(opts,'-l',val_l4)
   call opt_get(opts,'-i',val_i4)
   call opt_get(opts,'-r',val_r4)
   call opt_get(opts,'-d',val_r8)
   call opt_get(opts,'-s',val_str)
!
   if (opts%nargs.ge.1) then
     call opt_get_arg(opts,1,arg)
   else
     arg = ''
   endif
!   call opt_print(opts)
   call opt_finalize(opts)
!
   print *, ' '
   print *, 'logical : ', val_l4
   print *, 'integer : ', val_i4
   print *, 'real    : ', val_r4
   print *, 'double  : ', val_r8
   print *, 'string  : ', trim(val_str)
   print *, 'args(1) : ', trim(arg)
   print *, ' '
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
