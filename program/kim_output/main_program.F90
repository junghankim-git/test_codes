!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
   use kinds, only : i4, l4, r4, r8
   use opt_parser, only : opt_args_t, max_len_str_, opt_bool_, opt_store_
   use opt_parser, only : opt_logical_, opt_int_, opt_real_, opt_real8_, opt_string_
   use opt_parser, only : opt_initialize, opt_finalize, opt_add, opt_print
   use opt_parser, only : opt_get, opt_get_arg, opt_show_help
   use kim_std_output, only : kim_out_t
   use kim_std_output, only : kimout_initialize, kimout_finalize, kimout_run, kimout_print
!-------------------------------------------------------------------------------
   implicit none
!
   type(opt_args_t) :: myopt
   type(kim_out_t) :: myout
! files
   integer(i4), parameter :: max_nfiles = 2
   integer(i4) :: nfiles
   character(len=512) :: filename1
   character(len=512) :: filename2
! variables
   integer(i4) :: nvars
   character(len=32), dimension(100) :: varnames
   logical(l4) :: usedb
!
   call opt_initialize(myopt)
   call opt_add(myopt,'-u','--usedb',opt_bool_,.false.,'using db')
!
   nfiles = myopt%nargs
   print*,'nfiles = ',nfiles
!
   if (nfiles==1) then
     call opt_get_arg(myopt,1,filename1)
   elseif (nfiles==2) then
     call opt_get_arg(myopt,1,filename1)
     call opt_get_arg(myopt,2,filename2)
     print*,trim(filename1)
     print*,trim(filename2)
   else
     print*,'# of arguments was not matched...'
     print*,'loading namelist'
     call readnl()
   endif
!
   call opt_get(myopt,'-u',usedb)
!
   call opt_finalize(myopt)
!
   if (nfiles==1) then
     call kimout_initialize(myout,filename1,nvars=nvars,varnames=varnames,usedb=usedb)
   elseif (nfiles==2) then
     call kimout_initialize(myout,filename1,filename2,nvars=nvars,varnames=varnames)
   endif
!
   call kimout_run(myout)
!
   call kimout_print(myout,0)
   !call kimout_print(myout,1)
!
   call kimout_finalize(myout)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine readnl()
!-------------------------------------------------------------------------------
   implicit none
!
   namelist/inputs/nfiles, filename1, filename2, nvars, varnames
! local variables
   integer(i4) :: err
!
   ! open
   open(21,file='./inputs.nl',status='old',iostat=err)
   if (err.ne.0) then
     print*,'some error in open namelist...',err
     stop
   endif
   ! read
   read(21,nml=inputs,iostat=err)
   if (err.ne.0) then
     print*,'some error in read namelist...',err
     stop
   endif
   print*,'nfiles    = ',nfiles
   print*,'filename1 = ',trim(filename1)
   print*,'filename2 = ',trim(filename2)
   print*,'nvars     = ',nvars
   print*,'varnames  = ',varnames(1:nvars)
   ! close
   close(21)
   print*,' '
!
   return
   end subroutine readnl
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end
!-------------------------------------------------------------------------------
