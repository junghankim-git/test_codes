!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
   use kinds
   use read_hybrid_file, only : read_coefficient
   use std_atmosphere, only : std_atm_initialize,                              &
                              get_temperature, get_height, get_pressure, std_ps
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), parameter :: nlayers = 100
   integer(i4), parameter :: nlevs = 137
   real(r8), dimension(nlayers) :: h, t
   real(r8), dimension(nlevs+1) :: coefa, coefb, eta, p, hght, temp
   real(r8) :: p0
!
   integer(i4) :: unitnum, testcase
   integer(i4) :: i, k
!
   unitnum = 21
!
   print*,'######## start standard atmosphere ########'
   print*,' '
!
   unitnum = 21
   testcase = 7
!
   call std_atm_initialize(2)
   ! height - temperature
   do k = 1,nlayers
     h(k) = dble((k-1)*1000)
     t(k) = get_temperature(h(k))
   enddo
!
   open(unitnum,file = 'H-T.ascii',form = 'formatted',status = 'unknown')
   do k = 1,nlayers
     write(unitnum,*) h(k),t(k)
   enddo
   close(unitnum)
!
   call read_coefficient('/home/jhkim/work/program/hybrid_coord/dev/coefficient/ecmwf/ecmwf_half_0137nlevs.ascii',nlevs,coefa,coefb,eta,.true.)
!
   ! pressure - height
   p0 = 100000.0_r8
   do k = 1,nlevs+1
     p(k) = coefa(k)*p0+coefb(k)*std_ps
     if (k==1) p(k) = 0.5_r8*(coefa(k)+coefa(k+1))*p0+0.5_r8*(coefb(k)+coefb(k+1))*std_ps
     hght(k) = get_height(p(k))
     temp(k) = get_temperature(hght(k))
     print*,'P = ',p(k),get_pressure(hght(k))
     print*,'H = ',hght(k)
     print*,'T = ',temp(k)
     print*,' '
   enddo
!
   open(unitnum,file = 'P-H.ascii',form = 'formatted',status = 'unknown')
   do k = 1,nlevs+1
     write(unitnum,*) p(k),hght(k)
   enddo
   close(unitnum)
!
   end program test
