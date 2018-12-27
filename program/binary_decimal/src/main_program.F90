!-------------------------------------------------------------------------------
   program test
!
   use kinds,          only: i4, l4, r4, r8
   use binary_decimal, only: nbit32, nbit64, binary32_t, binary64_t,           &
                             bin_copy, bin_truncate, bin_add_bit, bin_remove_bit,&
                             dec_truncate, bin_roundoff, dec_roundoff,         &
                             dec_to_bin, bin_to_dec, bin_print
!                       dec_to_bin64, bin64_to_dec, bin64_copy, bin64_truncate, &
!                       dec_to_bin32, bin32_to_dec, bin32_copy, bin32_truncate, &
!                       dec_to_bin64_special,                                   &
!                       cut_64bit, cut_32bit, bin_print
   use netcdf
   implicit none
!
   type(binary32_t) :: bin32, bin32t
   type(binary64_t) :: bin64, bin64t
   real(r4) :: v4, w4
   real(r8) :: v8, w8
   integer(i4) :: i, n
   character(len=256) :: fname1, fname2
   character(len=64)  :: vname
!
!   call unit_test()
!   call test_roundoff_r8()
!   call test_5()
!   call test_64bit_32bit()
!
!fname1 = '/scratch/jhkim/TestBed/Data/3.2.01_up/10h/SW_G/1run/UP-20170601000300-000001.nc'
!fname2 = '/scratch/jhkim/TestBed/Data/3.2.01_up/10h/SW_G/2run/UP-20170601000300-000001.nc'
!fname1 = '/scratch/jhkim/TestBed/Data/3.2.01_up/10h/SW_G/1run/UP-20170601000600-000002.nc'
!fname2 = '/scratch/jhkim/TestBed/Data/3.2.01_up/10h/SW_G/2run/UP-20170601000600-000002.nc'
!fname1 = '/scratch/jhkim/TestBed/Data/3.2.01_up/10h/SW_G/1run/UP-20170601000900-000003.nc'
!fname2 = '/scratch/jhkim/TestBed/Data/3.2.01_up/10h/SW_G/2run/UP-20170601000900-000003.nc'
!fname1 = '/scratch/jhkim/TestBed/Data/3.2.01_up/10h/SW_G/1run/UP-20170601005700-000019.nc'
!fname2 = '/scratch/jhkim/TestBed/Data/3.2.01_up/10h/SW_G/2run/UP-20170601005700-000019.nc'
!fname1 = '/scratch/jhkim/TestBed/Data/3.2.01_up/10h/SW_G/1run/UP-20170601010000-000020.nc'
!fname2 = '/scratch/jhkim/TestBed/Data/3.2.01_up/10h/SW_G/2run/UP-20170601010000-000020.nc'
!fname1 = '/data/jhkim/TestBed/KIM/Output/3.2.01/ne30/gnu/10h/SW_G/101/UP-20170601000000-000000.nc'
!fname2 = '/scratch/jhkim/TestBed/Data/3.2.01.01/10h/SW_G/101/UP-20170601000000-000000.nc'
fname1 = '/scratch/jhkim/TestBed/Data/3.2.01.01/10h/SW_G/1run/UP-20170601000600-000002.nc'
fname2 = '/scratch/jhkim/TestBed/Data/3.2.01.01/10h/SW_G/2run/UP-20170601000600-000002.nc'
!fname2 = '/scratch/jhkim/TestBed/Data/3.2.01.01/10h/SW_G/restart/restart-20170601000300-000001.nc'
    !vname  = 'rw_tend'
    !vname  = 'thetap_m'
    !vname  = 'oz'
    !vname  = 't1_dyn'
    !vname  = 'frontgf'
    vname  = 'p_hyd_i'
    call check_kim_file(fname1,fname2,vname)
!
#if 0
   !v8 = -1.431243_r8
   !v8 = 3.33333333333333348_r8
   !v8 = 4.99999999999998579_r8
   v8 = 0.9999999999999999_r8
   !v8 = 4.99999999999999700_r8
   !v8 = 5.00000000000000300_r8
   v8 = 0.3892932132891738291_r8
   w8 = v8
   print *,v8
   call dec_to_bin(v8,bin64)
   call bin_print(bin64)
   call bin_copy(bin64,bin64t)
   call bin_roundoff(bin64t,3)
   call bin_print(bin64t)
   call bin_to_dec(bin64t,v8)
   print *,v8
   print *,' '
#endif
#if 0
   do i = 1,53-24
     call bin_truncate(bin64t,i)
     call bin_to_dec(bin64t,v8)
     call dec_truncate(w8,i)
     print *,v8,w8
   enddo
   !v4 = -0.999999_r8
   v4 = 0.999999_r8
   w4 = v4
   print *,v4
   call dec_to_bin(v4,bin32)
   call bin_copy(bin32,bin32t)
   call bin_to_dec(bin32t,v4)
   print *,v4
   do i = 1,12
     call bin_truncate(bin32t,i)
     call bin_to_dec(bin32t,v4)
     call dec_truncate(w4,i)
     print *,v4,w4
   enddo
#endif
!
#if 0
   !v8 = 0.99999_r8
   print *, 'fraction = ',fraction(v8)
   print *, 'radix    = ',radix(v8)
   print *, 'exponent = ',exponent(v8)
   call bin_print(bin64)
   print *,v8
   call bin64_to_dec(bin64,v8)
   print *,v8

   n = 5
   print *,n
   n = rshift(n,1)
   print *,n
   n = lshift(n,1)
   print *,n
   print *,v8
   v8 = rshift(v8,1)
   print *,v8
   v8 = lshift(v8,1)
   print *,v8
   print *,radix(v8)
#endif
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine unit_test()
   implicit none
! local variables
   type(binary32_t) :: bin32, bin32t
   type(binary64_t) :: bin64, bin64t
   real(r4) :: v4, w4
   real(r8) :: v8, w8, dv
   integer(i4) :: i, n, nroff
!
   v8 = 2.0
   call dec_to_bin(v8,bin64)
   print *,v8
   call bin_print(bin64)
   print *,' '

   v8 = 1.0
   call dec_to_bin(v8,bin64)
   print *,v8
   call bin_print(bin64)
   print *,' '

   v8 = 0.0
   call dec_to_bin(v8,bin64)
   print *,v8
   call bin_print(bin64)
   print *,' '

   v8 = 0.5
   call dec_to_bin(v8,bin64)
   print *,v8
   call bin_print(bin64)
   print *,' '

   v8 = 0.25
   call dec_to_bin(v8,bin64)
   print *,v8
   call bin_print(bin64)
   print *,' '

   v8 = 0.25
   call dec_to_bin(v8,bin64)
   print *,v8
   call bin_print(bin64)
   call bin_remove_bit(bin64,53)
   call bin_to_dec(bin64,v8)
   print *,v8
   call bin_print(bin64)
!
print *,'-------'
   nroff = 6
   dv = 1.0d-11

   v8 = +dv
   print *,v8
   call dec_roundoff(v8,nroff)
   print *,v8
   v8 = -dv
   print *,v8
   call dec_roundoff(v8,nroff)
   print *,v8
!
#if 0
   nroff = 3
   n = 100
   v8 = -0.3892932132891738291_r8
   do i = 1,n
     w8 = v8*real(100*(100-i),r8)+200.0_r8
     print *,w8
     call dec_to_bin(w8,bin64)
     call bin_print(bin64)
     call bin_copy(bin64,bin64t)
     call bin_roundoff(bin64t,nroff)
     call bin_print(bin64t)
     call bin_to_dec(bin64t,w8)
     print *,w8
     print *,''
   enddo
#endif
!
   end subroutine unit_test 
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine test_roundoff_r8()
   implicit none
! local variables
   real(r4)    :: v32
   real(r8)    :: v64
   real(r8)    :: frac, vmin, vmax
   integer(i4) :: nroff, n, ndiff
   logical(l4) :: prt, rnd
   type(binary32_t) :: bin32
   type(binary64_t) :: bin64
!
!
#if 1
   !nroff = 10
   nroff = 30
   n    = 1000000
   frac = 0.5d-15
   prt  = .false.
   rnd  = .true.
!   vmin = -10.0d0; vmax=10.0d0
   vmin = 0.0d0; vmax=1.0d0
!   vmin = 0.0d-10; vmax=1.0d-10
!   vmin = 0.0d100; vmax=1.0d100
!   vmin = 0.0d-100; vmax=1.0d-100
   call roundoff_r8(nroff,frac,prt,rnd,n,vmin,vmax,ndiff)
   write(*,'(a)')'========================================================='
   write(*,'(a,i5,a,f15.11,a)')'| ndiff = ', ndiff,  &
                         ', diff ratio = ',real(ndiff,r8)/real(n,r8)*100_r8,'%'
   write(*,'(a)')'========================================================='
#else
   nroff = 12
   n    = 100
   frac = 0.5d-15
   prt  = .false.
   rnd  = .true.
   vmin = 0.0d0; vmax=1.0d0
!   vmin = 0.0d0; vmax=1.0d0
!   vmin = 0.0d-10; vmax=1.0d-10
!   vmin = -1.0d100; vmax=1.0d100
!   vmin = -1.0d-100; vmax=1.0d-100
   call roundoff_r8(nroff,frac,prt,rnd,n,vmin,vmax,ndiff)
#endif
!
#if 0
   bin32%sig = .true.  
   bin32%exp = 0.0
   bin32%fra = 0.0
   bin32%bit(:) = 0
   bin32%bit(1) = 1
   bin32%bit(nbit32-5) = 1
   !bin32%bit(nbit32-1:nbit32) = 1
   call bin_to_dec(bin32,v32)
   print *, v32
   bin64%sig = .true.  
   bin64%exp = 0.0
   bin64%fra = 0.0
   bin64%bit(:) = 0
   bin64%bit(1) = 1
   bin64%bit(nbit64-5) = 1
   !bin64%bit(nbit64-1:nbit64) = 1
   call bin_to_dec(bin64,v64)
   print *, v64
#endif
!
   end subroutine test_roundoff_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine roundoff_r8(nroff, frac, prt, rnd, n, vmin, vmax, ndiff)
   implicit none
   integer(i4), intent(in   ) :: nroff
   real(r8)   , intent(in   ) :: frac
   logical(l4), intent(in   ) :: prt, rnd
   integer(i4), intent(in   ) :: n
   real(r8)   , intent(in   ) :: vmin, vmax
   integer(i4), intent(  out) :: ndiff
! local variables
   type(binary64_t) :: bin64, bin64v11, bin64v12, bin64v21, bin64v22
   real(r8)         :: var, v11, v12, v21, v22
   integer(i4)      :: i
   character(len=32) :: nstr, fstr
!
   write(nstr,'(a,i2.2,a)') '(a',59-nroff,',i2)'
   write(fstr,'(a,i2.2,a)') '(a',59-nroff,',a)'
!   print *,trim(fstr)
!
   ndiff = 0
   do i = 1,n
     call random_number(var)
     var = (vmax-vmin)*var+vmin
     if (prt) print *, var
     v11 = var+var*frac
     v21 = var-var*frac
     v12 = v11
     v22 = v21
!
     call dec_roundoff(v12,nroff)
     call dec_roundoff(v22,nroff)
!
     call dec_to_bin(var,bin64)
     call dec_to_bin(v11,bin64v11)
     call dec_to_bin(v12,bin64v12)
     call dec_to_bin(v21,bin64v21)
     call dec_to_bin(v22,bin64v22)
     !if (var.eq.v12.and.var.eq.v22) then
     if (v12.eq.v22) then
      if (prt) then
       print *,'-----------------------------------------------------'
       call bin_print(bin64)
       print *,'v1,v2 are same...'
       !print *,'round-off = ',nroff
       print trim(nstr),' ',nroff
       print trim(fstr),' ','|'
      ! print *,trim(fstr)
       call bin_print(bin64v11)
       call bin_print(bin64v21)
       print *,'   v11 = ',v11
       print *,'   v21 = ',v21
       call bin_print(bin64v12)
       call bin_print(bin64v22)
       print *,'   v12 = ',v12
       print *,'   v22 = ',v22
       print *,'-----------------------------------------------------'
       print *,' '
      endif
     else
       ndiff = ndiff+1
       print *,'-----------------------------------------------------'
       call bin_print(bin64)
       print *,'v1,v2 are diff...'
       !print *,'round-off = ',nroff
       print trim(nstr),' ',nroff
       print trim(fstr),' ','|'
      ! print *,trim(fstr)
       call bin_print(bin64v11)
       call bin_print(bin64v21)
       print *,bin64v11%fra
       print *,bin64v21%fra
       print *,'   v11 = ',v11
       print *,'   v21 = ',v21
       call bin_print(bin64v12)
       call bin_print(bin64v22)
       print *,bin64v12%fra
       print *,bin64v22%fra
       print *,'   v12 = ',v12
       print *,'   v22 = ',v22
       print *,'-----------------------------------------------------'
       print *,' '
     endif
   enddo
!
   end subroutine roundoff_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine test_5()
   implicit none
! local variables
   type(binary64_t) :: bin64, bin64t
   real(r8) :: v8, w8, dv
   integer(i4) :: i, n, nroff
!
   !    0.74987868499999932
   !+0.101111111111100000001100 10101100000010101110000001011
   v8 = 0.74987868499999999_r8
   call dec_to_bin(v8,bin64)
   print *,v8
   call bin_print(bin64)

   v8 = 0.749878685_r8
   call dec_to_bin(v8,bin64)
   print *,v8
   call bin_print(bin64)
!
   w8 = real(int(v8*10_r8**8,r8),r8)
   w8 = w8/10_r8**8
   print *,w8
   print *,v8-w8
!
   return
   end subroutine test_5
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine test_64bit_32bit()
   implicit none
! local variables
   type(binary32_t) :: bin32
   type(binary64_t) :: bin64
   real(r4) :: v4, w4
   real(r8) :: v8, w8
   integer(i4) :: i, n, nroff
!
   bin64%sig = .true.
   bin64%exp = 0
   bin64%bit(:) = 1
!
   call bin_to_dec(bin64,v8)
   print *,v8
   call bin_print(bin64)
   w8 = v8
   call dec_roundoff(w8,30)
   print *,w8
   print *,v8
   w8 = real(v8,r4)
   print *,w8
!
   return
   end subroutine test_64bit_32bit
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine check_kim_file(fname1, fname2, vname)
   implicit none
   character(len=*), intent(in   ) :: fname1, fname2, vname
! local variables
   integer(i4) :: i, k, ndiff
   integer(i4) :: err, fid, did, vid
   integer(i4) :: ncol, nlev
   real(r8) :: verr, vsum
   real(r8), dimension(:,:), allocatable :: v1, v2
   type(binary64_t) :: bin1, bin2
!
!  file 1
   err = nf90_open(trim(fname1),nf90_nowrite,fid)
   call nc_check(err)
   ! dim
   err = nf90_inq_dimid(fid,'ncol',did)
   call nc_check(err)
   err = nf90_inquire_dimension(fid,did,len=ncol)
   call nc_check(err)
   err = nf90_inq_dimid(fid,'nlev',did)
   call nc_check(err)
   err = nf90_inquire_dimension(fid,did,len=nlev)
   call nc_check(err)
   ! allocate
   allocate(v1(ncol,nlev+1))
   allocate(v2(ncol,nlev+1))
   ! var
   err = nf90_inq_varid(fid,trim(vname),vid)
   call nc_check(err)
   err = nf90_get_var(fid,vid,v1)
   call nc_check(err)
   err = nf90_close(fid)
   call nc_check(err)
!
!  file 2
   err = nf90_open(trim(fname2),nf90_nowrite,fid)
   call nc_check(err)
   err = nf90_inq_varid(fid,trim(vname),vid)
   call nc_check(err)
   err = nf90_get_var(fid,vid,v2)
   call nc_check(err)
   err = nf90_close(fid)
   call nc_check(err)
!
! algorithm
   ndiff = 0
   !do k = 1,nlev
   do k = 1,1
     do i = 1,ncol
#if 0
         call dec_to_bin(v1(i,k),bin1)
         call dec_to_bin(v2(i,k),bin2)
         print *,i,k,v1(i,k)
         print *,i,k,v2(i,k)
         call bin_print(bin1)
         call bin_print(bin2)
#endif
       if (v1(i,k).ne.v2(i,k)) then
         ndiff = ndiff+1
         call dec_to_bin(v1(i,k),bin1)
         call dec_to_bin(v2(i,k),bin2)
         print *,i,k,v1(i,k)
         print *,i,k,v2(i,k)
         call bin_print(bin1)
         call bin_print(bin2)
         print *,' '
       endif
     enddo
   enddo
   print*,ncol,ncol*nlev,ndiff
   vsum = 0.0_r8  
   verr = 0.0_r8
   do k=1,nlev+1
   do i=1,ncol
     call ddpdd(vsum,verr,v1(i,k))
   enddo
   enddo
   print*,minval(v1),maxval(v1),vsum
!
   deallocate(v1)
   deallocate(v2)
!
   return
   end subroutine check_kim_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_check(ncstat_in)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: ncstat_in
!
   if (ncstat_in.ne.nf90_noerr) then
     print *, trim(nf90_strerror(ncstat_in))
     stop
   endif
!
   return
   end subroutine nc_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine scs(add_a, err, add_b)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(inout) :: add_a, err
   real(r8), intent(in   ) :: add_b
! local variables
   real(r8) :: tsum, new_err, tmp
!-------------------------------------------------------------------------------
!
   tmp = err+add_b
   tsum = add_a+tmp
   new_err = tmp-(tsum-add_a)
   add_a = tsum
   err = new_err
!
   return
   end subroutine scs
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ddpdd(add_a, err, add_b)
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
   end subroutine ddpdd
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program test
