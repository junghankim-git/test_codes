!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    201x-0x-xx  your   name    code comment
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kinds
   use gmres_solver
   use ioascii
   use reprosum, only: scs, ddpdd
   use polynomial, only: inipoly, finpoly, argpoly, poly, diffpoly, intpoly
   implicit none
! local
!
   integer(i4), parameter :: n = 10
   real(r8), dimension(n+1) :: eta_i, psi_i, arg_i
   real(r8), dimension(n) :: eta_m, psi_m, arg_m
   real(r8), dimension(n+1, n+1) :: lu_i, luinv_i
   real(r8), dimension(n, n) :: lu_m, luinv_m
   real(r8), dimension(n+1, n+1) :: ainv_i, tmp_i
   real(r8), dimension(n, n) :: ainv_m, tmp_m
   real(r8), dimension(n+1) :: ipivot_i
   real(r8), dimension(n) :: ipivot_m
   real(r8), dimension(n+1) :: work_i
   real(r8), dimension(n) :: work_m
   integer(i4) :: ip, ie, ilev, info, i, j
   real(r8), dimension(n, n+1) :: inthf
   real(r8), dimension(n+1) :: intlag_i
   real(r8), dimension(n) :: intlag_m
!
   integer(i4) :: unitnum, testcase
   integer(i4), parameter :: f = 10, nwave = 2
   integer(i4), parameter :: nx_i = 10*(n+1)
   integer(i4), parameter :: nx_m = 10*n
   real(r8) :: errval, tmp, pi
   real(r8), dimension(nx_i) :: x_i, y_i
   real(r8), dimension(nx_m) :: x_m, y_m
! analytic
   real(r8), dimension(nx_i) :: diff_i, int_i
   real(r8), dimension(nx_m) :: diff_m, int_m
   real(r8), dimension(nx_i) :: intpsi_i
   real(r8), dimension(nx_m) :: intpsi_m
!
   print*,'######## start high order fdm setting ########'
   print*,' '
   pi = 3.14159265358979323846264_r8
   unitnum = 21
   testcase = 1
   call openfile('Result.dat',unitnum)
! eta
   do i = 1,n+1
     eta_i(i) = 1.0_r8/dble(n)*dble(i-1)
!eta_i(i)       = (1.0_r8/DBLE(n)*DBLE(i-1))**3.0_r8
   enddo
   do i = 1,n
     eta_m(i) = 0.5_r8*(eta_i(i)+eta_i(i+1))
   enddo
! psi
   do i = 1,n+1
!psi = 1
     if (testcase==1) then
       psi_i(i) = 1.0_r8
       if (i<n+1) psi_m(i) = 1.0_r8
!psi = x
     elseif (testcase==2) then
       psi_i(i) = eta_i(i)
       if (i<n+1) psi_m(i) = eta_m(i)
!psi = x^2
     elseif (testcase==3) then
       psi_i(i) = eta_i(i)*eta_i(i)
       if (i<n+1) psi_m(i) = eta_m(i)*eta_m(i)
! random
     elseif (testcase==4) then
       call random_number(psi_i(i))
       call random_number(psi_m(i))
!psi = sin(2pix)
     elseif (testcase==5) then
       psi_i(i) = dsin(2.0_r8*pi*nwave*eta_i(i))
       if (i<n+1) psi_m(i) = dsin(2.0_r8*pi*nwave*eta_m(i))
!psi = cos(2pix)
     elseif (testcase==6) then
       psi_i(i) = dcos(2.0_r8*pi*nwave*eta_i(i))
       if (i<n+1) psi_m(i) = dcos(2.0_r8*pi*nwave*eta_m(i))
     endif
   enddo
! check w
!           psi_i(1) = 0.0000000000000000_r8
!           psi_m(1) = 6.5534202333208391_r8
!           psi_i(2) = 5.7405187981472920_r8
!           psi_m(2) = 1.2784221009337351_r8
!           psi_i(3) =-3.4060469384031729_r8
!           psi_m(3) =-6.1726991875604824_r8
!           psi_i(4) =-5.5110040653391081_r8
!           psi_m(4) =-1.4520160693505386_r8
!           psi_i(5) = 3.8745432420133969_r8
!           psi_m(5) = 6.6341222138065961_r8
!           psi_i(6) = 0.0000000000000000_r8
! check F
!    psi_i(1) =  2.5338479573414562D0
!    psi_i(2) =    2.4086441477996119D0
!    psi_i(3) =    2.1688013038166813D0
!    psi_i(4) =    1.8325110114866204D0
!    psi_i(5) =    1.4281246489304665D0
!    psi_i(6) =    9.9168410351054694D-1
!    psi_i(7) =    5.6388554648404937D-1
!    psi_i(8) =    1.8685848932159918D-1
!    psi_i(9) =   -9.8596695336323008D-02
!    psi_i(10) =    -2.5441338009296833D-01
!    psi_i(11) =    -2.4429672745187760D-01
   print*,'Step 1'
   call inipoly(n,eta_i,eta_m)
   call argpoly(n+1,psi_i)
   call argpoly(n,psi_m)
   print*,'Step 2'
   call writedata_1int1(unitnum,n)
   call writedata_1int1(unitnum,f)
   call writedata_1int1(unitnum,testcase)
   call writedata_string(unitnum,'IO:dimension info.')
! Write eta_i,psi_i
   call writedata_2real(unitnum,n+1,eta_i,psi_i)
   call writedata_2real(unitnum,n,eta_m,psi_m)
   call writedata_string(unitnum,'IO:input data')
   print*,'Step 3'
! ### range Polynomial Test ( Interface )
! x
   do i = 1,nx_i
     x_i(i) = dble(i-1)/dble(nx_i-1)
   enddo
   call writedata_1real(unitnum,nx_i,x_i)
   call writedata_string(unitnum,'IO:x_i(half levels) ')
   print*,'Step 3-1'
! L_j(x)
   y_i(:) = 0.0_r8
   do j = 1,n+1
!      DO i = 1, nx_i
!        y_i(i) = Poly(n+1, x_i(i))
!      END DO
     call writedata_1real(unitnum,nx_i,y_i)
   enddo
   call writedata_string(unitnum,'IO:lag_i(half levels) ')
   print*,'Step 3-2'
! psi(\eta_i )
   y_i(:) = 0.0_r8
   do i = 1,nx_i! n
     y_i(i) = poly(n+1,x_i(i))
   enddo
   call writedata_1real(unitnum,nx_i,y_i)
   call writedata_string(unitnum,'IO:psi_i(half levels) ')
   print*,'Step 3-3'
! dL_j(x) / dx
   y_i(:) = 0.0_r8
   do j = 1,n+1
!      DO i = 1, nx_i
!        y_i(i) = DiffPoly(n+1, x_i(i))
!      END DO
     call writedata_1real(unitnum,nx_i,y_i)
   enddo
   call writedata_string(unitnum,'IO:diff_lag_i(half levels) ')
   print*,'Step 3-4'
! dpsi(x) / dx
   y_i(:) = 0.0_r8
   do i = 1,nx_i
!      DO j = 1, n+1
     y_i(i) = diffpoly(n+1,x_i(i))
!      END DO
   enddo
   call writedata_1real(unitnum,nx_i,y_i)
   call writedata_string(unitnum,'IO:diff_psi_i(half levels) ')
   print*,'Step 3-5'
! int{L_m(x)} dx
   do j = 1,n+1
     int_i(i) = 0.0_r8
!      DO i = 1,n+1
!        Int_i(i) = IntPoly(n+1,x_i(i))
!      END DO
!      CALL WriteData_1Real(unitnum,n+1,Int_i)
   enddo
!    CALL WriteData_String(unitnum,'IO: int_lag_i (half levels)')
   print*,'Step 3-6'
! int{psi(x)} dx
   do i = 1,n+1
     intpsi_i(i) = 0.0_r8
     do j = 1,nx_i
!        IntPsi_i(i) = IntPsi_i(i) + psi_i(j)*IntPoly(n+1, j, i)
     enddo
   enddo
!    CALL WriteData_1Real(unitnum,n+1,IntPsi_i)
!    CALL WriteData_String(unitnum,'IO: int_psi_i (half levels)')
   print*,'Step 3-7'
! int{psi(x)} dx
   intpsi_i(:) = 0.0_r8
   do i = 1,nx_i
     intpsi_i(i) = intpoly(n+1,x_i(i))
   enddo
   call writedata_1real(unitnum,nx_i,intpsi_i)
   call writedata_string(unitnum,'IO:int_psi_i(half levels) ')
!    DO i = 1,n
!      psi_m(i) = 0.0_r8
!      DO j = 1,n+1
!        psi_m(i) = psi_m(i) + psi_i(j)*Poly(n+1,j,eta_m(i))
!      END DO
!    END DO
!
!    PRINT *,' '
!    PRINT *,'### Check Interpolation(range Polynomial) : '
!    PRINT *,' '
!    PRINT *,'Half = '
!    DO j = 1,n
!      WRITE(*,'(F17.14,X)') psi_i(j)
!      WRITE(*,'(F17.14,X)') psi_m(j)
!    ENDDO
!      WRITE(*,'(F17.14,X)') psi_i(n+1)
   print*,'Step 4'
! ### range Polynomial Test ( Model )
! x
   do i = 1,nx_m
     x_m(i) =(eta_m(n)-eta_m(1))*dble(i-1)/dble(nx_m-1)+eta_m(1)
   enddo
   call writedata_1real(unitnum,nx_m,x_m)
   call writedata_string(unitnum,'IO:x_m(full levels) ')
! L_j(x)
   do j = 1,n
     do i = 1, nx_m
       y_m(i) = poly(n,x_m(i))
     enddo
     call writedata_1real(unitnum,nx_m,y_m)
   enddo
   call writedata_string(unitnum,'IO:lag_m(full levels) ')
! psi(\eta_i )
   y_m(:) = 0.0_r8
!    DO j = 1,n
   do i = 1,nx_m
     y_m(i) = poly(n,x_m(i))
   enddo
!    END DO
   call writedata_1real(unitnum,nx_m,y_m)
   call writedata_string(unitnum,'IO:psi_m(full levels) ')
! dL_j(x) / dx
   y_m(:) = 0.0_r8
   do j = 1,n
!      DO i = 1, nx_m
!        y_m(i) = DiffPoly(n, j, x_m(i))
!      END DO
     call writedata_1real(unitnum,nx_m,y_m)
   enddo
   call writedata_string(unitnum,'IO:diff_lag_m(full levels) ')
! dpsi(x) / dx
   y_m(:) = 0.0_r8
   do i = 1,nx_m
!      DO j = 1, n
     y_m(i) = diffpoly(n,x_m(i))
!      END DO
   enddo
   call writedata_1real(unitnum,nx_m,y_m)
   call writedata_string(unitnum,'IO:diff_psi_m(full levels) ')
! int{L_m(x)} dx
   do j = 1,n
     int_m(:) = 0.0_r8
!      DO i = 1,n
!        Int_m(i) = psi_m(j)*IntPoly(n,j,i)
!      END DO
!      CALL WriteData_1Real(unitnum,n,Int_m)
   enddo
!    CALL WriteData_String(unitnum,'IO: int_lag_m (full levels)')
! int{psi(x)} dx
   intpsi_m(:) = 0.0_r8
   do i = 1,nx_m
     intpsi_m(i) = intpoly(n,x_m(i))
   enddo
   call writedata_1real(unitnum,nx_m,intpsi_m)
   call writedata_string(unitnum,'IO:int_psi_m(full levels) ')
   print*,'Step 5'
! analytic differential,integration
   do i = 1,nx_i
!psi = 1
     if (testcase==1) then
       diff_i(i) = 0.0_r8
!psi = x
     elseif (testcase==2) then
       diff_i(i) = 1.0_r8
!psi = x^2
     elseif (testcase==3) then
       diff_i(i) = 2.0_r8*x_i(i)
!random
     elseif (testcase==4) then
       diff_i(i) = 0.0_r8
!psi = sin(2pix)
     elseif (testcase==5) then
       diff_i(i) = 2.0_r8*pi*nwave*dcos(2.0_r8*pi*nwave*x_i(i))
!psi = cos(2pix)
     elseif (testcase==6) then
       diff_i(i) =-2.0_r8*pi*nwave*dsin(2.0_r8*pi*nwave*x_i(i))
     endif
   enddo
   do i = 1,nx_m
!psi = 1
     if (testcase==1) then
       diff_m(i) = 0.0_r8
!psi = x
     elseif (testcase==2) then
       diff_m(i) = 1.0_r8
!psi = x^2
     elseif (testcase==3) then
       diff_m(i) = 2.0_r8*x_m(i)
!random
     elseif (testcase==4) then
       diff_m(i) = 0.0_r8
!psi = sin(2pix)
     elseif (testcase==5) then
       diff_m(i) = 2.0_r8*pi*nwave*dcos(2.0_r8*pi*nwave*x_m(i))
!psi = cos(2pix)
     elseif (testcase==6) then
       diff_m(i) =-2.0_r8*pi*nwave*dsin(2.0_r8*pi*nwave*x_m(i))
     endif
   enddo
   do i = 1,nx_i
!psi = 1
     if (testcase==1) then
       int_i(i) = x_i(i)
!psi = x
     elseif (testcase==2) then
       int_i(i) = 0.5_r8*x_i(i)**2.0_r8
!psi = x^2
     elseif (testcase==3) then
       int_i(i) =(x_i(i)**3.0_r8)/3.0_r8
!random
     elseif (testcase==4) then
       int_i(i) = 0.0_r8
!psi = sin(2pix)
     elseif (testcase==5) then
       int_i(i) =-dcos(2.0_r8*pi*nwave*x_i(i))/(2.0_r8*pi*nwave)+              &
                                                       1.0_r8/(2.0_r8*pi*nwave)
!psi = cos(2pix)
     elseif (testcase==6) then
       int_i(i) = dsin(2.0_r8*pi*nwave*x_i(i))/(2.0_r8*pi*nwave)
     endif
   enddo
   do i = 1,nx_m
!psi = 1
     if (testcase==1) then
       int_m(i) = x_m(i)
!psi = x
     elseif (testcase==2) then
       int_m(i) = 0.5_r8*x_m(i)**2.0_r8
!psi = x^2
     elseif (testcase==3) then
       int_m(i) =(x_m(i)**3.0_r8)/3.0_r8
!random
     elseif (testcase==4) then
       int_m(i) = 0.0_r8
!psi = sin(2pix)
     elseif (testcase==5) then
       int_m(i) =-dcos(2.0_r8*pi*nwave*x_m(i))/(2.0_r8*pi*nwave)+              &
                                                       1.0_r8/(2.0_r8*pi*nwave)
!psi = cos(2pix)
     elseif (testcase==6) then
       int_m(i) = dsin(2.0_r8*pi*nwave*x_m(i))/(2.0_r8*pi*nwave)
     endif
   enddo
   call writedata_1real(unitnum,nx_i,diff_i)
   call writedata_string(unitnum,'IO:analytic diff.(half levels) ')
   call writedata_1real(unitnum,nx_i,int_i)
   call writedata_string(unitnum,'IO:analytic int.(half levels) ')
   call writedata_1real(unitnum,nx_m,diff_m)
   call writedata_string(unitnum,'IO:analytic diff.(full levels) ')
   call writedata_1real(unitnum,nx_m,int_m)
   call writedata_string(unitnum,'IO:analytic int.(full levels) ')
   call finpoly()
   call closefile(unitnum)
   end program test
