!-------------------------------------------------------------------------------
   program test
   use kinds,  only: i4, r8
   use io_ascii, only: open_file, close_file, write_ascii, write_ascii_lines
   use cubic_spline_interp, only: cs_interp_initialize, cs_interp_finalize, cs_interp_get
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), parameter :: n = 10
   real(r8), dimension(n) :: eta_i, psi_i
   real(r8), dimension(10*n) :: reta_i, rpsi_i, apsi_i
   real(r8) :: pi, dp1, dpn
   integer(i4), parameter :: f = 10, nwave = 2
!
   integer(i4) :: unitnum, testcase
   integer(i4) :: i
!
   print*,'######## start cubic spline interpolation ########'
   print*,' '
!
   pi = 3.14159265358979323846264_r8
   unitnum = 21
   testcase = 5
!
   if (testcase==1) then
   dp1 = 0.0_r8
   dpn = 0.0_r8
   elseif (testcase==2) then
   dp1 = 1.0_r8
   dpn = 1.0_r8
   elseif (testcase==3) then
   dp1 = 0.0_r8
   dpn = 2.0_r8
   elseif (testcase==4) then
   dp1 = 1.0_r8
   dpn = 1.0_r8
   elseif (testcase==5) then
   dp1 = 2.0_r8*pi*nwave
   dpn = 2.0_r8*pi*nwave
   !dp1 = 0.0_r8
   !dpn = 0.0_r8
   elseif (testcase==6) then
   dp1 = 0.0_r8
   dpn = 0.0_r8
   elseif (testcase==7) then
   dp1 = 0.0_r8
   dpn = 4.0_r8
   endif
! eta
   do i = 1,n
     eta_i(i) = dble(i-1)/dble(n-1)
   enddo
! psi: input
   do i = 1,n
     !psi = 1
     if (testcase==1) then
       psi_i(i) = 1.0_r8
     !psi = x
     elseif (testcase==2) then
       psi_i(i) = eta_i(i)
     !psi = x^2
     elseif (testcase==3) then
       psi_i(i) = eta_i(i)*eta_i(i)
     ! random
     elseif (testcase==4) then
       call random_number(psi_i(i))
     !psi = sin(2pix)
     elseif (testcase==5) then
       psi_i(i) = dsin(2.0_r8*pi*nwave*eta_i(i))
     !psi = cos(2pix)
     elseif (testcase==6) then
       psi_i(i) = dcos(2.0_r8*pi*nwave*eta_i(i))
     !psi = x^4
     elseif (testcase==7) then
       psi_i(i) = eta_i(i)*eta_i(i)*eta_i(i)*eta_i(i)
     endif
   enddo
!
! analytic solution
   do i = 1,10*n
     reta_i(i) = dble(i-1)/dble(10*n-1)
     if (testcase==1) then
       apsi_i(i) = 1.0_r8
     !psi = x
     elseif (testcase==2) then
       apsi_i(i) = reta_i(i)
     !psi = x^2
     elseif (testcase==3) then
       apsi_i(i) = reta_i(i)*reta_i(i)
     ! random
     elseif (testcase==4) then
       call random_number(apsi_i(i))
     !psi = sin(2pix)
     elseif (testcase==5) then
       apsi_i(i) = dsin(2.0_r8*pi*nwave*reta_i(i))
     !psi = cos(2pix)
     elseif (testcase==6) then
       apsi_i(i) = dcos(2.0_r8*pi*nwave*reta_i(i))
     !psi = x^4
     elseif (testcase==7) then
       apsi_i(i) = reta_i(i)*reta_i(i)*reta_i(i)*reta_i(i)
     endif
   enddo
!
   call cs_interp_initialize(n,eta_i,psi_i,dp1,dpn)
!
   rpsi_i = 0.0_r8
   do i = 1,10*n
     rpsi_i(i) = cs_interp_get(reta_i(i))
   enddo
!
   call open_file('Result.dat',unitnum)
   call write_ascii(unitnum,n)
   call write_ascii(unitnum,10*n)
   call write_ascii(unitnum,nwave)
   call write_ascii(unitnum,testcase)
   call write_ascii(unitnum,'IO : dimension info.')
   call write_ascii_lines(unitnum,n,eta_i,psi_i)
   call write_ascii(unitnum,'IO : input points')
   call write_ascii_lines(unitnum,10*n,reta_i,rpsi_i)
   call write_ascii(unitnum,'IO : output points')
   call write_ascii_lines(unitnum,10*n,reta_i,apsi_i)
   call write_ascii(unitnum,'IO : analytic solution')
   call close_file(unitnum)
!
   call cs_interp_finalize()
!
   end program test
!-------------------------------------------------------------------------------
