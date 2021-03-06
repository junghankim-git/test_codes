#include "KIM.h"
!-------------------------------------------------------------------------------
   module globalnorms
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
   use kiapsbase, only: iulog=>kim_iu_log, real_kind=>kim_real8_kind
   use edge, only: edgebuffer_t
   implicit none
!
   private
   public :: l1_snorm
   public :: l2_snorm
   public :: linf_snorm
   public :: l1_vnorm
   public :: l2_vnorm
   public :: linf_vnorm
!  public :: print_cfl
   public :: test_global_integral
   public :: global_integral
   public :: wrap_repro_sum
   private :: global_maximum
   type(edgebuffer_t), private :: edgebuf
!
   contains
! ================================
! global_integral:
!
! eq 81 in Williamson, et. al. p 218
! for spectral elements
!
! ================================
! --------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function global_integral(elem, h, hybrid, npts, nets, nete) result(i_sphere)
!use hybrid_mod,  only : hybrid_t
!-------------------------------------------------------------------------------
   use kiapsparallel, only: hybrid_t
   use element, only: element_t
   use dimensions, only: np, nelemd
   use physicalconstants, only: dd_pi
   use kiapsparallel, only: global_shared_buf, global_shared_sum
!
   type(element_t), intent(in   ) :: elem(:)
   integer, intent(in   ) :: npts, nets, nete
   real(real_kind), intent(in   ) :: h(npts, npts, nets:nete)
   type(hybrid_t), intent(in   ) :: hybrid
   real(real_kind) :: i_sphere
   real(real_kind) :: i_priv
   real(real_kind) :: i_shared
   common/gblintcom/i_shared
! Local variables
   integer :: ie, j, i
   real(real_kind) :: i_tmp(1)
   real(real_kind) :: da
   real(real_kind) :: j_tmp(nets:nete)
!
! This algorythm is independent of thread count and task count.
! This is a requirement of consistancy checking in cam.
!
!
   j_tmp = 0.0d0
   do ie = nets,nete
     do j = 1, np
       do i = 1, np
         da = elem(ie)%mp(i,j)*elem(ie)%metdet(i,j)
         j_tmp(ie) = j_tmp(ie)+da*h(i,j,ie)
       enddo
       enddo
   enddo
   do ie = nets,nete
     global_shared_buf(ie,1) = j_tmp(ie)
   enddo
   call wrap_repro_sum(nvars = 1,comm = hybrid%par%comm)
   i_tmp = global_shared_sum(1)
   i_sphere = i_tmp(1)/(4.0d0*dd_pi)
!
   end function global_integral
! ================================
! test_global_integral:
!
! test that the global integral of
! the area of the sphere is 1.
!
! ================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine test_global_integral(elem, hybrid, nets, nete, mindxout)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: hybrid_t, abortpar, global_shared_buf, global_shared_sum
   use element, only: element_t
   use dimensions, only: np, ne, nelem, nelemd
   use mesh, only: meshusemeshfile
   use reductions, only: parallelmin, parallelmax
   use physicalconstants, only: rrearth=>kim_rrearth, rearth=>kim_rearth, dd_pi
   use edge, only: edgebuffer_t, initedgebuffer, freeedgebuffer, edgevpack, edgevunpack
   use bndry, only: bndry_exchangev
!
   type(element_t), intent(inout) :: elem(:)
   integer, intent(in   ) :: nets, nete
   type(hybrid_t), intent(in   ) :: hybrid
   real(real_kind), intent(  out), optional :: mindxout
   real(real_kind) :: i_sphere
! Local variables
   real(real_kind), allocatable :: h(:,:,:)
! Element statisics
   real(real_kind) :: min_area, max_area, avg_area
   real(real_kind) :: min_max_dx, max_max_dx, avg_max_dx
   real(real_kind) :: min_min_dx, max_min_dx, avg_min_dx
   real(real_kind) :: min_unif_dx, max_unif_dx
   real(real_kind) :: min_max_eig, max_max_eig, avg_max_eig
   real(real_kind) :: min_min_eig, max_min_eig, avg_min_eig
   real(real_kind) :: min_len
   integer :: ie, corner, i, j
!   MNL: use these variables for calculating largest length scale of element with
!        smallest dx (in uniform meshes, these are in the corner, with
!        max_dx = sqrt(3)*min_dx
!    real (kind=real_kind) :: max_dx_at_min_dx
!    integer :: max_dx_at_min_dx_id
!
   allocate(h(np,np,nets:nete))
   h(:,:,nets:nete) = 1.0d0
! Calculate surface area by integrating 1.0d0 over sphere and dividing by 4*DD_PI
! (Should be 1)
   i_sphere = global_integral(elem,h(:,:,nets:nete),hybrid,np,nets,nete)
   min_area = 1d99
   max_area = 0
   avg_area = 0_real_kind
   min_max_eig = 1d99
   max_max_eig = 0
   avg_max_eig = 0_real_kind
   min_min_eig = 1d99
   max_min_eig = 0
   avg_min_eig = 0_real_kind
   min_max_dx = 1d99
   max_max_dx = 0
   avg_max_dx = 0_real_kind
   min_min_dx = 1d99
   max_min_dx = 0
   avg_min_dx = 0_real_kind
   do ie = nets,nete
     elem(ie)%area = sum(elem(ie)%spheremp(:,:))
     min_area = min(min_area,elem(ie)%area)
     max_area = max(max_area,elem(ie)%area)
     min_max_eig = min(min_max_eig,elem(ie)%max_eig)
     max_max_eig = max(max_max_eig,elem(ie)%max_eig)
     min_min_eig = min(min_min_eig,elem(ie)%min_eig)
     max_min_eig = max(max_min_eig,elem(ie)%min_eig)
     min_max_dx = min(min_max_dx,elem(ie)%dx_long)
     max_max_dx = max(max_max_dx,elem(ie)%dx_long)
     min_min_dx = min(min_min_dx,elem(ie)%dx_short)
     max_min_dx = max(max_min_dx,elem(ie)%dx_short)
! MNL: Want this output from all the uniform meshes
     global_shared_buf(ie,1) = elem(ie)%area
     global_shared_buf(ie,2) = elem(ie)%max_eig
     global_shared_buf(ie,3) = elem(ie)%min_eig
     global_shared_buf(ie,4) = elem(ie)%dx_long
     global_shared_buf(ie,5) = elem(ie)%dx_short
!      (Useful for setting level of variable hypervis)
!      (For now,only run on 1 proc if you want correct output!)
!       if (elem(ie)%dx_short.eq.min_min_dx) then
!           max_dx_at_min_dx = elem(ie)%dx_long
!           max_dx_at_min_dx_id = elem(ie)%GlobalId
!       end if
   enddo
! MNL See comment above re: output from uniform meshes
!    print*,"Min dx_short: ",min_min_dx
!    print*,"Occurs at elem ",max_dx_at_min_dx_id," where dx_long = ",max_dx_at_min_dx
   min_area = parallelmin(min_area,hybrid)
   max_area = parallelmax(max_area,hybrid)
   min_max_eig = parallelmin(min_max_eig,hybrid)
   max_max_eig = parallelmax(max_max_eig,hybrid)
   min_min_eig = parallelmin(min_min_eig,hybrid)
   max_min_eig = parallelmax(max_min_eig,hybrid)
   min_max_dx = parallelmin(min_max_dx,hybrid)
   max_max_dx = parallelmax(max_max_dx,hybrid)
   min_min_dx = parallelmin(min_min_dx,hybrid)
   max_min_dx = parallelmax(max_min_dx,hybrid)
   call wrap_repro_sum(nvars = 5,comm = hybrid%par%comm)
   avg_area = global_shared_sum(1)/dble(nelem)
   avg_max_eig = global_shared_sum(2)/dble(nelem)
   avg_min_eig = global_shared_sum(3)/dble(nelem)
   avg_max_dx = global_shared_sum(4)/dble(nelem)
   avg_min_dx = global_shared_sum(5)/dble(nelem)
! Physical units for area
   min_area = min_area*rearth*rearth/1000000_real_kind
   max_area = max_area*rearth*rearth/1000000_real_kind
   avg_area = avg_area*rearth*rearth/1000000_real_kind
! for an equation du/dt = i c u,leapfrog is stable for |c u dt| < 1
! Consider a gravity wave at the equator,c=340m/s
! u = exp(i kmax x/ a ) with x = longitude, and kmax =  pi a / dx,
! u = exp(i pi x / dx ),  so du/dt = c du/dx becomes du/dt = i c pi/dx u
! stable for dt < dx/(c*pi)
! CAM 26 level AMIP simulation: max gravity wave speed 341.75 m/s
   if (hybrid%ismasterthr) then
     write(iulog,*) ""
     write(iulog,*) "running global integral diagnostic..."
     write(iulog,*) "area of unit sphere is",i_sphere
     write(iulog,*) "should be 1.0 to round off..."
     write(iulog,'(a,f7.3) ') 'Element area:max/min',(max_area/min_area)
     if (.not.meshusemeshfile) then
       write(iulog,'(a,f6.3,f8.2) ') "average equatorial node spacing(deg,km) = ",&
dble(90)/dble(ne*(np-1)),dd_pi*rearth/(2000.0d0*dble(ne*(np-1)))
     endif
     write(iulog,'(a,2f8.4) ') 'Min eigenvalue of dinv(min,max):',min_min_eig,max_min_eig
     write(iulog,'(a,2f8.4) ') 'Max eigenvalue of dinv(min,max):',min_max_eig,max_max_eig
     write(iulog,'(a,3f8.2) ') 'dx for cfl(smallest scale per elem):ave,min,max = ',avg_min_dx,min_min_dx,max_min_dx
     write(iulog,'(a,3f8.2) ') 'dx for hypervis(largest scale per elem):ave,min,max = ',avg_max_dx,min_max_dx,max_max_dx
   endif
   deallocate(h)
   if (present(mindxout)) then
! min_len now based on norm(Dinv)
     mindxout = 1000_real_kind*min_len
     min_len=0.002d0*rearth/(dble(np-1)*max_max_eig)
   endif
!
   end subroutine test_global_integral
!  ! ================================
!  ! print_cfl:
!  !
!  ! Calculate / output CFL info
!  ! (both advective and based on
!  ! viscosity or hyperviscosity)
!  !
!  ! ================================
!
!  subroutine print_cfl(elem,hybrid,nets,nete,dtnu,dt_advection)
!
!    USE KiapsParallel, ONLY: hybrid_t, AbortPar, global_shared_buf, global_shared_sum
!    use Element, only : element_t
!    use Dimensions, only : np,ne,nelem,nelemd
!    use Quadrature, only : gausslobatto, quadrature_t
!
!    USE Reductions, only : ParallelMin,ParallelMax
!    USE PhysicalConstants, ONLY : rrearth=>KIM_RREARTH, rearth=>KIM_REARTH,dd_pi
!    use Control, only : nu, hypervis_order, nu_top, hypervis_power, &
!                            fine_ne, rk_stage_user, max_hypervis_courant
!    use Edge, only : EdgeBuffer_t, initedgebuffer, FreeEdgeBuffer, edgeVpack, edgeVunpack
!    use Bndry, only : bndry_exchangeV
!    use TimeStep, only : tstep
!
!    type(element_t)      , intent(inout) :: elem(:)
!    integer              , intent(in) :: nets,nete
!    type (hybrid_t)      , intent(in) :: hybrid
!    real (kind=real_kind), intent(in) :: dtnu
!    real (kind=real_kind), intent(in) :: dt_advection
!
!    ! Element statisics
!    real (kind=real_kind) :: min_max_dx, max_max_dx, avg_max_dx
!    real (kind=real_kind) :: min_min_dx, max_min_dx, avg_min_dx
!    real (kind=real_kind) :: max_unif_dx
!    real (kind=real_kind) :: min_max_eig, max_max_eig, avg_max_eig
!    real (kind=real_kind) :: min_min_eig, max_min_eig, avg_min_eig
!    real (kind=real_kind) :: min_hypervis, max_hypervis, avg_hypervis, stable_hv
!    real (kind=real_kind) :: max_hvcourant
!    real (kind=real_kind) :: x, y, noreast, nw, se, sw
!    real (kind=real_kind), dimension(np,np,nets:nete) :: zeta
!    real (kind=real_kind) :: lambda_max, lambda_vis, min_gw
!    real (kind=real_kind) :: dt,dt_eul
!    integer :: ie,corner, i, j
!    type (quadrature_t)    :: gp
!
!    ! Eigenvalues calculated by folks at UMich (Paul U & Jared W)
!    select case (np)
!        case (2)
!            lambda_max = 0.5d0
!        case (3)
!            lambda_max = 1.5d0
!            lambda_vis = 12.0d0
!        case (4)
!            lambda_max = 2.74d0
!            lambda_vis = 30.0d0
!        case (5)
!            lambda_max = 4.18d0
!            lambda_vis = 91.6742d0
!        case (6)
!            lambda_max = 5.86d0
!            lambda_vis = 190.1176d0
!        case (7)
!            lambda_max = 7.79d0
!            lambda_vis = 374.7788d0
!        case (8)
!            lambda_max = 10.0d0
!            lambda_vis = 652.3015d0
!        case DEFAULT
!            lambda_max = 0.0d0
!            lambda_vis = 0.0d0
!    end select
!
!    if ((lambda_max.eq.0d0).and.(hybrid%isMasterThr)) then
!        print*, "lambda_max not calculated for NP = ",np
!        print*, "Estimate of gravity wave timestep will be incorrect"
!    end if
!
!    do ie=nets,nete
!      elem(ie)%variable_hyperviscosity = 1;
!    end do
!
!    gp=gausslobatto(np)
!    min_gw = minval(gp%weights)
!
!    min_max_eig=1d99
!    max_max_eig=0
!    avg_max_eig=0_real_kind
!
!    min_min_eig=1d99
!    max_min_eig=0
!    avg_min_eig=0_real_kind
!
!    min_max_dx=1d99
!    max_max_dx=0
!    avg_max_dx=0_real_kind
!
!    min_min_dx=1d99
!    max_min_dx=0
!    avg_min_dx=0_real_kind
!
!    do ie=nets,nete
!        min_max_eig = min(min_max_eig,elem(ie)%max_eig)
!        max_max_eig = max(max_max_eig,elem(ie)%max_eig)
!
!        min_min_eig = min(min_min_eig,elem(ie)%min_eig)
!        max_min_eig = max(max_min_eig,elem(ie)%min_eig)
!
!        min_max_dx = min(min_max_dx,elem(ie)%dx_long)
!        max_max_dx = max(max_max_dx,elem(ie)%dx_long)
!
!        min_min_dx = min(min_min_dx,elem(ie)%dx_short)
!        max_min_dx = max(max_min_dx,elem(ie)%dx_short)
!
!       global_shared_buf(ie,1) = elem(ie)%max_eig
!       global_shared_buf(ie,2) = elem(ie)%min_eig
!       global_shared_buf(ie,3) = elem(ie)%dx_long
!       global_shared_buf(ie,4) = elem(ie)%dx_short
!    enddo
!
!    min_max_eig=ParallelMin(min_max_eig,hybrid)
!    max_max_eig=ParallelMax(max_max_eig,hybrid)
!
!    min_min_eig=ParallelMin(min_min_eig,hybrid)
!    max_min_eig=ParallelMax(max_min_eig,hybrid)
!
!    min_max_dx=ParallelMin(min_max_dx,hybrid)
!    max_max_dx=ParallelMax(max_max_dx,hybrid)
!
!    min_min_dx=ParallelMin(min_min_dx,hybrid)
!    max_min_dx=ParallelMax(max_min_dx,hybrid)
!
!    call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
!
!    avg_max_eig = global_shared_sum(1)/dble(nelem)
!    avg_min_eig = global_shared_sum(2)/dble(nelem)
!    avg_max_dx = global_shared_sum(3)/dble(nelem)
!    avg_min_dx = global_shared_sum(4)/dble(nelem)
!
!    if (hypervis_power /= 0) then
!
!        min_hypervis = 1d99
!        max_hypervis = 0
!        avg_hypervis = 0
!
!!       viscosity in namelist specified for smallest element:
!        max_unif_dx = min_max_dx
!        if (fine_ne>0) then
!           ! viscosity in namelist specified for regions with a resolution
!           ! equivilant to a uniform grid with ne=fine_ne
!           max_unif_dx = (111.28*30)/dble(fine_ne)
!        endif
!
!        dt = tstep
!        max_hvcourant = 0
!
!        do ie=nets,nete
!           ! variable viscosity based on map from ulatlon -> ucontra
!
!           ! dx_long
!           elem(ie)%variable_hyperviscosity = sqrt((elem(ie)%dx_long/max_unif_dx) ** hypervis_power)
!
!           ! dx_short
!           !elem(ie)%variable_hyperviscosity = sqrt((elem(ie)%dx_short/1.25/min_min_dx) ** hypervis_power)
!
!           ! geometric mean:
!           !elem(ie)%variable_hyperviscosity = sqrt(  (sqrt(elem(ie)%dx_short*elem(ie)%dx_long)/min_min_dx/1.9) ** hypervis_power)
!
!           ! average
!           !elem(ie)%variable_hyperviscosity = sqrt( ( (elem(ie)%dx_short+elem(ie)%dx_long)/3.8/min_min_dx) ** hypervis_power)
!
!           ! sqrt(area)
!           ! area = unit sphere.  other variables in km
!           !elem(ie)%variable_hyperviscosity = sqrt( ( sqrt(elem(ie)%area)*Rearth/2.4e3/max_unif_dx) ** hypervis_power)
!
!
!           elem(ie)%hv_courant = dtnu*(elem(ie)%variable_hyperviscosity(1,1)**2) * (lambda_vis**2) * ((rrearth*elem(ie)%max_eig)**4)
!
!
!            ! Check to see if this is stable
!            if (elem(ie)%hv_courant.gt.max_hypervis_courant) then
!                stable_hv = sqrt( max_hypervis_courant / &
!                     (  dtnu * (lambda_vis)**2 * (rrearth*elem(ie)%max_eig)**4 ) )
!
!#if 0
!        ! Useful print statements for debugging the adjustments to hypervis
!                print*, "Adjusting hypervis on elem ", elem(ie)%GlobalId
!                print*, "From ", nu*elem(ie)%variable_hyperviscosity(1,1)**2, " to ", nu*stable_hv
!                print*, "Difference = ", nu*(/elem(ie)%variable_hyperviscosity(1,1)**2-stable_hv/)
!                print*, "Factor of ", elem(ie)%variable_hyperviscosity(1,1)**2/stable_hv
!                print*, " "
!#endif
!
!!                make sure that: elem(ie)%hv_courant <=  max_hypervis_courant
!                elem(ie)%variable_hyperviscosity = stable_hv
!                elem(ie)%hv_courant = dtnu*(stable_hv**2) * (lambda_vis)**2 * (rrearth*elem(ie)%max_eig)**4
!            end if
!            max_hvcourant = max(max_hvcourant, elem(ie)%hv_courant)
!
!            min_hypervis = min(min_hypervis, elem(ie)%variable_hyperviscosity(1,1))
!            max_hypervis = max(max_hypervis, elem(ie)%variable_hyperviscosity(1,1))
!            global_shared_buf(ie,1) = elem(ie)%variable_hyperviscosity(1,1)
!        end do
!
!        min_hypervis = ParallelMin(min_hypervis, hybrid)
!        max_hypervis = ParallelMax(max_hypervis, hybrid)
!        call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
!        avg_hypervis = global_shared_sum(1)/dble(nelem)
!
!        max_hvcourant = ParallelMax(max_hvcourant, hybrid)
!
!        ! apply DSS (aka assembly procedure) to variable_hyperviscosity (makes continuous)
!        call initEdgeBuffer(edgebuf,1)
!        do ie=nets,nete
!            zeta(:,:,ie) = elem(ie)%variable_hyperviscosity(:,:)*elem(ie)%spheremp(:,:)
!            call edgeVpack(edgebuf,zeta(1,1,ie),1,0,elem(ie)%desc)
!        end do
!        call bndry_exchangeV(hybrid,edgebuf)
!        do ie=nets,nete
!            call edgeVunpack(edgebuf,zeta(1,1,ie),1,0,elem(ie)%desc)
!            elem(ie)%variable_hyperviscosity(:,:) = zeta(:,:,ie)*elem(ie)%rspheremp(:,:)
!        end do
!        call FreeEdgeBuffer(edgebuf)
!
!        ! replace hypervis w/ bilinear based on continuous corner values
!        do ie=nets,nete
!            noreast = elem(ie)%variable_hyperviscosity(np,np)
!            nw = elem(ie)%variable_hyperviscosity(1,np)
!            se = elem(ie)%variable_hyperviscosity(np,1)
!            sw = elem(ie)%variable_hyperviscosity(1,1)
!            do i=1,np
!                x = gp%points(i)
!                do j=1,np
!                    y = gp%points(j)
!                    elem(ie)%variable_hyperviscosity(i,j) = 0.25d0*( &
!                                            (1.0d0-x)*(1.0d0-y)*sw + &
!                                            (1.0d0-x)*(y+1.0d0)*nw + &
!                                            (x+1.0d0)*(1.0d0-y)*se + &
!                                            (x+1.0d0)*(y+1.0d0)*noreast)
!                end do
!            end do
!        end do
!    else
!        max_hvcourant = dtnu * (lambda_vis)**2 * (rrearth*max_max_eig)**4
!    end if
!
!    deallocate(gp%points)
!    deallocate(gp%weights)
!
!     if (hybrid%isMasterThr) then
!       if (rk_stage_user>0) then
!          dt_eul = (dt_advection/(rk_stage_user-1)) ! stepsize of internal euler step in RKSSP method
!          write(iulog,'(a,f10.2)') 'SSP Forward Euler (120m/s) Courant number=', (120.0d0*max_max_eig*rrearth*dt_eul)/min_gw
!       endif
!       write(iulog,'(a,f10.2)') 'advective (120m/s) Courant number=', 120.0d0*max_max_eig*lambda_max*rrearth*dt_advection
!       write(iulog,'(a,f10.2)') 'gravity wave (342m/s) Courant number=', 342.0d0*max_max_eig*lambda_max*rrearth*dt_advection
!       if (nu>0) then
!          if (hypervis_order==1) then
!              write(iulog,'(a,f10.2)') 'max viscosity courant number: ',dtnu*((rrearth*max_max_eig)**2)*lambda_vis
!          endif
!          if (hypervis_order==2) then
!             write(iulog,'(a,f12.6)') "max hyper viscosity courant number: ", max_hvcourant
!          endif
!       endif
!       if(nu_top>0) then
!          write(iulog,'(a,f10.2,a)') 'TOP3 viscosity CFL: dt < S*',1.0d0/(4*nu_top*((rrearth*max_max_eig)**2)*lambda_vis),'s'
!       end if
!      if (hypervis_power /= 0) then
!         write(iulog,'(a,3e11.4)') 'Hyperviscosity (dynamics): ave,min,max = ',nu*(/avg_hypervis**2,min_hypervis**2,max_hypervis**2/)
!!         print*, 'fine_ne = ', fine_ne
!!         print*, 'Using max_unif_dx = ', max_unif_dx
!      end if
!    end if
!
!  end subroutine print_cfl
! ================================
! global_maximum:
!
! Find global maximum on sphere
!
! ================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function global_maximum(h, hybrid, npts, nets, nete) result(max_sphere)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: hybrid_t
   use reductions, only: red_max, pmax_mt
!
   integer, intent(in   ) :: npts, nets, nete
   real(real_kind), intent(in   ) :: h(npts, npts, nets:nete)
   type(hybrid_t), intent(in   ) :: hybrid
   real(real_kind) :: max_sphere
! Local variables
   real(real_kind) :: redp(1)
!
   max_sphere = maxval(h(:,:,nets:nete))
   redp(1) = max_sphere
   call pmax_mt(red_max,redp,1,hybrid)
   max_sphere = red_max%buf(1)
!
   end function global_maximum
! ==========================================================
! l1_snorm:
!
! computes the l1 norm per Williamson et al, p. 218 eq(8)
! for a scalar quantity
! ===========================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function l1_snorm(elem, h, ht, hybrid, npts, nets, nete) result(l1)
!-------------------------------------------------------------------------------
   use element, only: element_t
   use kiapsparallel, only: hybrid_t
!
   type(element_t), intent(in   ) :: elem(:)
   integer, intent(in   ) :: npts, nets, nete
   real(real_kind), intent(in   ) :: h(npts, npts, nets:nete) ! computed soln
   real(real_kind), intent(in   ) :: ht(npts, npts, nets:nete) ! true soln
   type(hybrid_t), intent(in   ) :: hybrid
   real(real_kind) :: l1
! Local variables
   real(real_kind) :: dhabs(npts, npts, nets:nete)
   real(real_kind) :: htabs(npts, npts, nets:nete)
   real(real_kind) :: dhabs_int
   real(real_kind) :: htabs_int
   integer i, j, ie
!
   do ie = nets, nete
     do j = 1, npts
       do i = 1, npts
         dhabs(i,j,ie) = abs(h(i,j,ie)-ht(i,j,ie))
         htabs(i,j,ie) = abs(ht(i,j,ie))
       enddo
       enddo
   enddo
!
   dhabs_int = global_integral(elem,dhabs(:,:,nets:nete),hybrid,npts,nets,nete)
   htabs_int = global_integral(elem,htabs(:,:,nets:nete),hybrid,npts,nets,nete)
   l1 = dhabs_int/htabs_int
!
   end function l1_snorm
! ===========================================================
! l1_vnorm:
!
! computes the l1 norm per Williamson et al, p. 218 eq(97),
! for a contravariant vector quantity on the velocity grid.
!
! ===========================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function l1_vnorm(elem, v, vt, hybrid, npts, nets, nete) result(l1)
!-------------------------------------------------------------------------------
   use element, only: element_t
!    use hybrid_mod, only : hybrid_t
   use kiapsparallel, only: hybrid_t
!
   type(element_t), intent(in   ), target :: elem(:)
   integer, intent(in   ) :: npts, nets, nete
   real(real_kind), intent(in   ) :: v(npts, npts, 2, nets:nete) ! computed soln
   real(real_kind), intent(in   ) :: vt(npts, npts, 2, nets:nete) ! true soln
   type(hybrid_t), intent(in   ) :: hybrid
   real(real_kind) :: l1
! Local variables
   real(real_kind), dimension(:,:,:,:), pointer :: met
   real(real_kind) :: dvsq(npts, npts, nets:nete)
   real(real_kind) :: vtsq(npts, npts, nets:nete)
   real(real_kind) :: dvco(npts, npts, 2) ! covariant velocity
   real(real_kind) :: vtco(npts, npts, 2) ! covariant velocity
   real(real_kind) :: dv1, dv2
   real(real_kind) :: vt1, vt2
   real(real_kind) :: dvsq_int
   real(real_kind) :: vtsq_int
   integer i, j, ie
!
   do ie = nets, nete
     met=>elem(ie)%met
     do j = 1,npts
       do i = 1, npts
         dv1 = v(i,j,1,ie)-vt(i,j,1,ie)
         dv2 = v(i,j,2,ie)-vt(i,j,2,ie)
         vt1 = vt(i,j,1,ie)
         vt2 = vt(i,j,2,ie)
         dvco(i,j,1) = met(1,1,i,j)*dv1+met(1,2,i,j)*dv2
         dvco(i,j,2) = met(2,1,i,j)*dv1+met(2,2,i,j)*dv2
         vtco(i,j,1) = met(1,1,i,j)*vt1+met(1,2,i,j)*vt2
         vtco(i,j,2) = met(2,1,i,j)*vt1+met(2,2,i,j)*vt2
         dvsq(i,j,ie) = sqrt(dvco(i,j,1)*dv1+dvco(i,j,2)*dv2)
         vtsq(i,j,ie) = sqrt(vtco(i,j,1)*vt1+vtco(i,j,2)*vt2)
       enddo
     enddo
   enddo
!
   dvsq_int = global_integral(elem,dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
   vtsq_int = global_integral(elem,vtsq(:,:,nets:nete),hybrid,npts,nets,nete)
   l1 = dvsq_int/vtsq_int
!
   end function l1_vnorm
! ==========================================================
! l2_snorm:
!
! computes the l2 norm per Williamson et al, p. 218 eq(83)
! for a scalar quantity on the pressure grid.
!
! ===========================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function l2_snorm(elem, h, ht, hybrid, npts, nets, nete) result(l2)
!-------------------------------------------------------------------------------
   use element, only: element_t
!    use hybrid_mod, only : hybrid_t
   use kiapsparallel, only: hybrid_t
!
   type(element_t), intent(in   ) :: elem(:)
   integer, intent(in   ) :: npts, nets, nete
   real(real_kind), intent(in   ) :: h(npts, npts, nets:nete) ! computed soln
   real(real_kind), intent(in   ) :: ht(npts, npts, nets:nete) ! true soln
   type(hybrid_t), intent(in   ) :: hybrid
   real(real_kind) :: l2
! Local variables
   real(real_kind) :: dh2(npts, npts, nets:nete)
   real(real_kind) :: ht2(npts, npts, nets:nete)
   real(real_kind) :: dh2_int
   real(real_kind) :: ht2_int
   integer i, j, ie
!
   do ie = nets, nete
     do j = 1, npts
       do i = 1, npts
         dh2(i,j,ie) =(h(i,j,ie)-ht(i,j,ie))**2
         ht2(i,j,ie) = ht(i,j,ie)**2
       enddo
       enddo
   enddo
!
   dh2_int = global_integral(elem,dh2(:,:,nets:nete),hybrid,npts,nets,nete)
   ht2_int = global_integral(elem,ht2(:,:,nets:nete),hybrid,npts,nets,nete)
   l2 = sqrt(dh2_int)/sqrt(ht2_int)
!
   end function l2_snorm
! ==========================================================
! l2_vnorm:
!
! computes the l2 norm per Williamson et al, p. 219 eq(98)
! for a contravariant vector quantity on the velocity grid.
!
! ===========================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function l2_vnorm(elem, v, vt, hybrid, npts, nets, nete) result(l2)
!-------------------------------------------------------------------------------
   use element, only: element_t
!    use hybrid_mod, only : hybrid_t
   use kiapsparallel, only: hybrid_t
!
   type(element_t), intent(in   ), target :: elem(:)
   integer, intent(in   ) :: npts, nets, nete
   real(real_kind), intent(in   ) :: v(npts, npts, 2, nets:nete) ! computed soln
   real(real_kind), intent(in   ) :: vt(npts, npts, 2, nets:nete) ! true soln
   type(hybrid_t), intent(in   ) :: hybrid
   real(real_kind) :: l2
! Local variables
   real(real_kind), dimension(:,:,:,:), pointer :: met
   real(real_kind) :: dvsq(npts, npts, nets:nete)
   real(real_kind) :: vtsq(npts, npts, nets:nete)
   real(real_kind) :: dvco(npts, npts, 2) ! covariant velocity
   real(real_kind) :: vtco(npts, npts, 2) ! covariant velocity
   real(real_kind) :: dv1, dv2
   real(real_kind) :: vt1, vt2
   real(real_kind) :: dvsq_int
   real(real_kind) :: vtsq_int
   integer i, j, ie
!
   do ie = nets, nete
     met=>elem(ie)%met
     do j = 1,npts
       do i = 1, npts
         dv1 = v(i,j,1,ie)-vt(i,j,1,ie)
         dv2 = v(i,j,2,ie)-vt(i,j,2,ie)
         vt1 = vt(i,j,1,ie)
         vt2 = vt(i,j,2,ie)
         dvco(i,j,1) = met(1,1,i,j)*dv1+met(1,2,i,j)*dv2
         dvco(i,j,2) = met(2,1,i,j)*dv1+met(2,2,i,j)*dv2
         vtco(i,j,1) = met(1,1,i,j)*vt1+met(1,2,i,j)*vt2
         vtco(i,j,2) = met(2,1,i,j)*vt1+met(2,2,i,j)*vt2
         dvsq(i,j,ie) = dvco(i,j,1)*dv1+dvco(i,j,2)*dv2
         vtsq(i,j,ie) = vtco(i,j,1)*vt1+vtco(i,j,2)*vt2
       enddo
     enddo
   enddo
!
   dvsq_int = global_integral(elem,dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
   vtsq_int = global_integral(elem,vtsq(:,:,nets:nete),hybrid,npts,nets,nete)
   l2 = sqrt(dvsq_int)/sqrt(vtsq_int)
!
   end function l2_vnorm
! ==========================================================
! linf_snorm:
!
! computes the l infinity norm per Williamson et al, p. 218 eq(84)
! for a scalar quantity on the pressure grid...
!
! ===========================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function linf_snorm(h, ht, hybrid, npts, nets, nete) result(linf)
!    use hybrid_mod, only : hybrid_t
!-------------------------------------------------------------------------------
   use kiapsparallel, only: hybrid_t
!
   integer, intent(in   ) :: npts, nets, nete
   real(real_kind), intent(in   ) :: h(npts, npts, nets:nete) ! computed soln
   real(real_kind), intent(in   ) :: ht(npts, npts, nets:nete) ! true soln
   type(hybrid_t), intent(in   ) :: hybrid
   real(real_kind) :: linf
! Local variables
   real(real_kind) :: dhabs(npts, npts, nets:nete)
   real(real_kind) :: htabs(npts, npts, nets:nete)
   real(real_kind) :: dhabs_max
   real(real_kind) :: htabs_max
   integer i, j, ie
!
   do ie = nets, nete
     do j = 1, npts
       do i = 1, npts
         dhabs(i,j,ie) = abs(h(i,j,ie)-ht(i,j,ie))
         htabs(i,j,ie) = abs(ht(i,j,ie))
       enddo
       enddo
   enddo
!
   dhabs_max = global_maximum(dhabs(:,:,nets:nete),hybrid,npts,nets,nete)
   htabs_max = global_maximum(htabs(:,:,nets:nete),hybrid,npts,nets,nete)
   linf = dhabs_max/htabs_max
!
   end function linf_snorm
! ==========================================================
! linf_vnorm:
!
! computes the linf norm per Williamson et al, p. 218 eq(99),
! for a contravariant vector quantity on the velocity grid.
!
! ===========================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function linf_vnorm(elem, v, vt, hybrid, npts, nets, nete) result(linf)
!    use hybrid_mod, only : hybrid_t
!-------------------------------------------------------------------------------
   use kiapsparallel, only: hybrid_t
   use element, only: element_t
!
   type(element_t), intent(in   ), target :: elem(:)
   integer, intent(in   ) :: npts, nets, nete
   real(real_kind), intent(in   ) :: v(npts, npts, 2, nets:nete) ! computed soln
   real(real_kind), intent(in   ) :: vt(npts, npts, 2, nets:nete) ! true soln
   type(hybrid_t), intent(in   ) :: hybrid
   real(real_kind) :: linf
! Local variables
   real(real_kind), dimension(:,:,:,:), pointer :: met
   real(real_kind) :: dvsq(npts, npts, nets:nete)
   real(real_kind) :: vtsq(npts, npts, nets:nete)
   real(real_kind) :: dvco(npts, npts, 2) ! covariant velocity
   real(real_kind) :: vtco(npts, npts, 2) ! covariant velocity
   real(real_kind) :: dv1, dv2
   real(real_kind) :: vt1, vt2
   real(real_kind) :: dvsq_max
   real(real_kind) :: vtsq_max
   integer i, j, ie
!
   do ie = nets, nete
     met=>elem(ie)%met
     do j = 1,npts
       do i = 1, npts
         dv1 = v(i,j,1,ie)-vt(i,j,1,ie)
         dv2 = v(i,j,2,ie)-vt(i,j,2,ie)
         vt1 = vt(i,j,1,ie)
         vt2 = vt(i,j,2,ie)
         dvco(i,j,1) = met(1,1,i,j)*dv1+met(1,2,i,j)*dv2
         dvco(i,j,2) = met(2,1,i,j)*dv1+met(2,2,i,j)*dv2
         vtco(i,j,1) = met(1,1,i,j)*vt1+met(1,2,i,j)*vt2
         vtco(i,j,2) = met(2,1,i,j)*vt1+met(2,2,i,j)*vt2
         dvsq(i,j,ie) = sqrt(dvco(i,j,1)*dv1+dvco(i,j,2)*dv2)
         vtsq(i,j,ie) = sqrt(vtco(i,j,1)*vt1+vtco(i,j,2)*vt2)
       enddo
     enddo
   enddo
!
   dvsq_max = global_maximum(dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
   vtsq_max = global_maximum(vtsq(:,:,nets:nete),hybrid,npts,nets,nete)
   linf = dvsq_max/vtsq_max
!
   end function linf_vnorm
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine wrap_repro_sum(nvars, comm, nsize)
!-------------------------------------------------------------------------------
   use dimensions, only: nelemd
   use reprosummod, only: reprosum
   use kiapsparallel, only: global_shared_buf, global_shared_sum, nrepro_vars, abortpar
   implicit none
!
   integer :: nvars ! number of variables to be summed(cannot exceed nrepro_vars)
   integer :: comm ! mpi communicator
   integer, optional :: nsize ! local buffer size(defaults to nelemd-number of elements in mpi task)
   integer nsize_use
!
   if (present(nsize)) then
     nsize_use = nsize
   else
     nsize_use = nelemd
   endif
!
   if (nvars.gt.nrepro_vars) call abortpar(message = 'ERROR:repro_sum_buffer_size exceeded')
! Repro_sum contains its own OpenMP,so only one thread should call it (AAM)
#if (! defined 	ELEMENT_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
   call reprosum(global_shared_buf,global_shared_sum,nsize_use,nelemd,nvars,commid = comm)
#if (! defined 	ELEMENT_OPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif
!
   end subroutine wrap_repro_sum
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module globalnorms
!-------------------------------------------------------------------------------
