!-------------------------------------------------------------------------------
   program test
!
   use kinds,          only: i4, l4, r4, r8
   use check_var_in_netcdf, only: do_check
   implicit none
!
   character(len=256) :: fname1, fname2
   character(len=64)  :: vname1, vname2
!
fname1 = '/scratch/jhkim/TestBed/Data/3.2a_h/10h/SW_G/101/UP-20170601100000-000010.nc'
fname2 = '/scratch/jhkim/TestBed/Data/3.2a_h/10h/SW_G/101/UP-20170601100000-000010.nc'
vname1  = 'hgt'
vname2  = 'hgt_new'
call do_check(fname1,fname2,vname1,vname2)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine tmp(add_a, err, add_b)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(inout) :: add_a, err
   real(r8), intent(in   ) :: add_b
! local variables
   real(r8) :: tsum, new_err
!-------------------------------------------------------------------------------
!
   tsum = add_a+add_b
   new_err = (add_b-(tsum-add_a))+(add_a-(tsum-(tsum-add_a)))
   call scs(tsum,err,new_err)
   add_a = tsum
!
   return
   end subroutine tmp
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program test
