!-------------------------------------------------------------------------------
   module cubed_sphere
!-------------------------------------------------------------------------------
!
!  abstract : cubed-sphere grid generate module
!
!  history log :
!    2015-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds,        only: i4, l4, r4, r8
   use mesh_control, only: ww, ee, ss, nn, ws, es, wn, en,                     &
                           axis_x, axis_y, axis_xy, axis_mxy, axis_0,          &
                           mesh_convert, mesh_dir_convert, mesh_print
   use coordinates,  only: deg2rad, rad2deg, zero_thr,                         &
                           coordinate_t, cartesian2lonlat, lonlat2cartesian,   &
                           rotation_lonlat_ll0_to_axis, & !inv_rotation_lonlat_ll0_to_axis, &
                           rotation_lonslats_ll0_to_axis, inv_rotation_lonslats_ll0_to_axis, &
                           rotation_lonlat_xaxis_to_ll0,                       &
                           lonlat_to_great_circle_angle, great_circle_angle_to_lonlat
   use quadrature,   only: qp_t=>quadrature_t, gausslobatto
   use bilinear,     only: bilinear_t, bilinear_set, bilinear_interpolation
   use scrip_input,  only: scrip_t, none_, scrip_initialize, scrip_finalize,   &
                           scrip_check, scrip_write
   use space_filling_curve, only: sfc_t, sfc_initialize, sfc_finalize, sfc_gen_curve, sfc_get_sfc
   use netcdf
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   real(r8)   , parameter :: pi = acos(-1.0_r8)
   integer(i4), parameter :: nnbrs = 8
   integer(i4), parameter :: nvrts = 4 ! # of element's vertices
   integer(i4), parameter :: n_nbr = 5 ! # of info in neighbor(iface, ie, je, global, rank)
   integer(i4), parameter :: ncell = 4 ! # of maximum cell vetices
!
   character(len=2), dimension(nnbrs), parameter :: nbr_str =(/'ww','ee','ss','nn','ws','es','wn','en'/)
!
   type nbr_quadrature_t
     real(r8) :: lon, lat
     real(r8) :: x, y, z
   end type nbr_quadrature_t
!
   type quadrature_t
     real(r8) :: alpha, beta
     real(r8) :: lon, lat
     real(r8) :: x, y, z
     logical(l4) :: isunique = .false.
     logical(l4) :: hasunique = .false.
     ! for unique point
     integer(i4) :: iup, u_face, u_ie, u_je, u_ip, u_jp
     type(nbr_quadrature_t), dimension(nnbrs) :: nbrs
     type(nbr_quadrature_t), dimension(ncell) :: cell
   end type quadrature_t
!
   type neighbor_t
     integer(i4) :: iface  ! face number
     integer(i4) :: ie, je ! element index(in face)
     integer(i4) :: global ! global index
     integer(i4) :: rank   ! process number
   end type neighbor_t
!
   type element_t
     type(neighbor_t), dimension(nnbrs) :: nbrs ! neighbors
     integer(i4) :: iface ! face number
     integer(i4) :: global ! global index
     integer(i4) :: local !
     ! index
     integer(i4) :: rank ! process number
     integer(i4) :: ie, je ! element index(in face)
     integer(i4) :: isfc ! sfc index
     integer(i4) :: joiner ! joiner's direction
     logical(l4) :: isout ! outer element
     real(r8), dimension(nvrts) :: alpha, beta ! ld, rd, ru, lu
     real(r8), dimension(nvrts) :: lon, lat ! ld, rd, ru, lu
     type(quadrature_t), dimension(:,:), allocatable :: qp
   end type element_t
!
   type cubed_sphere_t
     integer(i4) :: nface = 6
     integer(i4) :: np, ne, nprocs
     integer(i4) :: nelem
     integer(i4) :: nnbrs = nnbrs
     integer(i4) :: ncell = ncell
     integer(i4) :: pnx, pny
     integer(i4) :: nvrts = nvrts ! # of vertices
     integer(i4), dimension(:), allocatable :: nelemd, nelemo ! # of procs in each elements
     integer(i4) :: nups, neps
     logical(l4) :: isrotated
     logical(l4) :: uselonlatinfo
     character(len=512) :: lonlatfilename
     character(len= 16) :: lonvarname
     character(len= 16) :: latvarname
     real(r8)       , dimension(:)    , allocatable :: qp
     type(element_t), dimension(:,:,:), allocatable :: elem
   end type cubed_sphere_t
!
   interface get_value_for_plot
     module procedure get_value_for_plot_i4_2d
     module procedure get_value_for_plot_i4_4d
     module procedure get_value_for_plot_qp
     module procedure get_value_for_plot_r8_3d
   end interface get_value_for_plot
!
   public :: cubed_sphere_t, cs_initialize, cs_finalize, cs_run
   public :: cs_print_coordinates, cs_print_directions
   public :: cs_write_cube_infos, cs_write_lonslats, cs_write_scrip
   public :: cs_do_check, cs_get_distance, cs_find_index_point
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_initialize(cs, np, ne, nprocs, isrotated, uselonlatinfo, lonlatfilename, lonvarname, latvarname)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t),       intent(inout) :: cs
   integer(i4),                intent(in   ) :: np, ne, nprocs
   logical(l4),      optional, intent(in   ) :: isrotated
   logical(l4),      optional, intent(in   ) :: uselonlatinfo
   character(len=*), optional, intent(in   ) :: lonlatfilename
   character(len=*), optional, intent(in   ) :: lonvarname, latvarname
! local variables
   integer(i4) :: iface, ie, je
   type(qp_t)  :: qp
   logical(l4) :: islonlatinfo
!
   cs%nface = 6
   cs%np = np
   cs%ne = ne
   cs%nelem = cs%nface*ne*ne
   cs%nprocs = nprocs
   cs%pnx = 4*ne
   cs%pny = 3*ne
   if (present(isrotated)) then
     cs%isrotated = isrotated
   else
     cs%isrotated = .false.
   endif
   allocate(cs%nelemd(nprocs))
   allocate(cs%nelemo(nprocs))
   allocate(cs%qp(np))
   allocate(cs%elem(ne,ne,cs%nface))
   do iface = 1,6
   do je = 1,ne
   do ie = 1,ne
     allocate(cs%elem(ie,je,iface)%qp(np,np))
   enddo
   enddo
   enddo
!
   qp = gausslobatto(np)
   cs%qp(:) = qp%points(:)
!
   deallocate(qp%points,qp%weights)
!
   if (present(uselonlatinfo)) then
     islonlatinfo = uselonlatinfo
   else
     islonlatinfo = .false.
   endif
!
   if (islonlatinfo.and.present(lonlatfilename).and.present(lonvarname).and.present(latvarname)) then
     cs%uselonlatinfo  = .true.
     cs%lonlatfilename = lonlatfilename
     cs%lonvarname     = lonvarname
     cs%latvarname     = latvarname
   else
     cs%uselonlatinfo  = .false.
     cs%lonlatfilename = 'None'
     cs%lonvarname     = 'None'
     cs%latvarname     = 'None'
   endif
!
   return
   end subroutine cs_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_finalize(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: iface, ie, je
!
   deallocate(cs%nelemd)
   deallocate(cs%nelemo)
   deallocate(cs%qp)
   do iface = 1,6
   do je = 1,cs%ne
   do ie = 1,cs%ne
     deallocate(cs%elem(ie,je,iface)%qp)
   enddo
   enddo
   enddo
   deallocate(cs%elem)
!
   return
   end subroutine cs_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_run(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
!
   print*,'Set coordinates...'
   if (.not.cs%uselonlatinfo) then
   call set_coordinates(cs) ! iface,global index(global),cartesian index(ie,je),alpha,beta,lon,lat
   endif
   print*,'Set sfc index...'
   call set_sfc_index(cs) ! isfc,joiner
   print*,'Set domain decomposition...'
   call set_rank(cs) ! rank
   print*,'Set local index...'
   call set_local_index(cs) !
   print*,'Set neighbor index...'
   call set_neighbors(cs)
   print*,'Set unique and entire points...'
   call set_quadrature_ups_eps(cs) !
!
   print*,'Set quadrature points...'
   if (cs%uselonlatinfo) then
     call set_quadrature_coordinate_with_llfile(cs) ! not operated
   else
     call set_quadrature_coordinate(cs) !
   endif
   print*,'Set quadrature points neghbors...'
   call set_quadrature_neighbors(cs)
   print*,'Set quadrature points cells...'
   call set_quadrature_cells(cs)
!
   return
   end subroutine cs_run
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_coordinates(cs, usestr)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t),  intent(inout) :: cs
   logical(l4), optional, intent(in   ) :: usestr
! local variables
   integer(i4) :: ne, nface, iface, ie, je, i
   real(r8) :: sa, sb, da, db
   real(r8) :: a1, a2, b1, b2
   real(r8), dimension(:), allocatable :: da_s, db_s
   logical(l4) :: luse
!
   ne = cs%ne
   nface = cs%nface
!
   if (present(usestr)) then
     luse = usestr
   else
     luse = .false.
   endif
!
   ! face,global index,catesian index
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
   cs%elem(ie,je,iface)%iface = iface
   cs%elem(ie,je,iface)%global = ne*ne*(iface-1)+ne*(je-1)+(ie-1)+1
   cs%elem(ie,je,iface)%ie = ie
   cs%elem(ie,je,iface)%je = je
   enddo
   enddo
   enddo
!
   ! alpha,beta,lon,lat
   if (.not.luse) then
!
     sa =-pi/4.0_r8
     sb =-pi/4.0_r8
     da = pi/2.0_r8/real(ne,8)
     db = da
     do iface = 1,nface
       do je = 1,ne
       do ie = 1,ne
       a1 = sa+da*real(ie-1,8)
       a2 = sa+da*real(ie,8)
       b1 = sb+db*real(je-1,8)
       b2 = sb+db*real(je,8)
       cs%elem(ie,je,iface)%alpha =(/a1,a2,a2,a1/)
       cs%elem(ie,je,iface)%beta =(/b1,b1,b2,b2/)
       do i = 1,cs%nvrts
         call alphabeta2lonlat(iface,cs%elem(ie,je,iface)%alpha(i),cs%elem(ie,je,iface)%beta(i), &
                                     cs%elem(ie,je,iface)%lon(i),  cs%elem(ie,je,iface)%lat(i),  &
                               cs%isrotated)
       enddo
     enddo
     enddo
     enddo
!
   else
!
     allocate(da_s(ne),db_s(ne))
     ! face: 1-4
     do iface = 1,4
       sa =-pi/4.0_r8
       sb =-pi/4.0_r8
       da = pi/2.0_r8/real(ne,8)
       db = 5.0_r8*pi/8.0_r8/real(ne,8)
       do je = 1,ne
       do ie = 1,ne
         a1 = sa+da*real(ie-1,8)
         a2 = sa+da*real(ie,8)
         b1 = sb+db*real(je-1,8)
         b2 = sb+db*real(je,8)
         cs%elem(ie,je,iface)%alpha =(/a1,a2,a2,a1/)
         cs%elem(ie,je,iface)%beta =(/b1,b1,b2,b2/)
         do i = 1,cs%nvrts
         call alphabeta2lonlat(iface,cs%elem(ie,je,iface)%alpha(i),cs%elem(ie,je,iface)%beta(i), &
                                     cs%elem(ie,je,iface)%lon(i),  cs%elem(ie,je,iface)%lat(i),  &
                               cs%isrotated)
         enddo
       enddo
       enddo
     enddo
     ! face: 5,6
     do iface = 5,6
!
       if (iface==6) then
         sa =-pi/8.0_r8
         sb =-pi/8.0_r8
         da_s(:) = pi/4.0_r8/real(ne,8)
         db_s(:) = da_s(:)
       elseif (iface==5) then
         sa =-pi/4.0_r8
         sb =-pi/4.0_r8
         da_s(:) = pi/2.0_r8/real(ne,8)
         db_s(:) = da_s(:)
       endif
!
       do je = 1,ne
       do ie = 1,ne
         a1 = sa+da_s(ie)*real(ie-1,8)
         a2 = sa+da_s(ie)*real(ie,8)
         b1 = sb+db_s(je)*real(je-1,8)
         b2 = sb+db_s(je)*real(je,8)
         cs%elem(ie,je,iface)%alpha =(/a1,a2,a2,a1/)
         cs%elem(ie,je,iface)%beta =(/b1,b1,b2,b2/)
         do i = 1,cs%nvrts
         call alphabeta2lonlat(iface,cs%elem(ie,je,iface)%alpha(i),cs%elem(ie,je,iface)%beta(i), &
                                     cs%elem(ie,je,iface)%lon(i),  cs%elem(ie,je,iface)%lat(i),  &
                               cs%isrotated)
         enddo
       enddo
       enddo
     enddo
     deallocate(da_s,db_s)
!
   endif
!
   return
   end subroutine set_coordinates
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine alphabeta2lonlat(iface, alpha, beta, lon, lat, isrotated)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: iface
   real(r8),    intent(in   ) :: alpha, beta
   real(r8),    intent(  out) :: lon, lat
   logical(l4), optional, intent(in   ) :: isrotated
! local variables
   real(r8) :: x, y, rs, lon0, lat0, rlon, rlat
!
   x = tan(alpha)
   y = tan(beta)
   call xy2lonlat(iface,x,y,lon,lat)
   if (present(isrotated)) then
     if (isrotated) then
       lon0 = 127.0_r8*deg2rad
       lat0 =  38.0_r8*deg2rad
!       print *, 'not yet supprt rotation...'
!       stop
       call rotation_lonlat_xaxis_to_ll0(lon0,lat0,1.0_r8,lon,lat,rs,rlon,rlat)
       lon = rlon; lat = rlat
     endif
   endif
!
   return
   end subroutine alphabeta2lonlat
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine xy2lonlat(iface, x, y, lon, lat)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: iface
   real(r8),    intent(in   ) :: x, y
   real(r8),    intent(  out) :: lon, lat
!
   if (iface.ge.1.and.iface.le.4) then
     lon = atan(x)+pi/2.0_r8*real(iface-1,8)
     lat = atan(y*cos(lon-pi/2.0_r8*real(iface-1,8)))
   elseif (iface.eq.5) then
     if (x==0.0_r8.and.y>0.0_r8) then
       lon = 0.0_r8
       lat =-pi/2.0_r8+atan(y)
     elseif (x==0.0_r8.and.y<0.0_r8) then
       lon = pi
       lat =-pi/2.0_r8-atan(y)
     elseif (x>0.0_r8.and.y==0.0_r8) then
       lon = pi/2.0_r8
       lat =-pi/2.0_r8+atan(x)
     elseif (x<0.0_r8.and.y==0.0_r8) then
       lon = 3.0*pi/2.0_r8
       lat =-pi/2.0_r8-atan(x)
     elseif (x==0.0_r8.and.y==0.0_r8) then
       lon = 0.0_r8
       lat =-pi/2.0_r8
     elseif (x/= 0.0_r8.and.y>0.0_r8) then
       lon = atan(x/y)
       lat = atan(-sin(lon)/x)
     elseif (x/= 0.0_r8.and.y<0.0_r8) then
       lon = pi+atan(x/y)
       lat = atan(-sin(lon)/x)
     endif
   elseif (iface.eq.6) then
     if (x==0.0_r8.and.y>0.0_r8) then
       lon = pi
       lat = pi/2.0_r8-atan(y)
     elseif (x==0.0_r8.and.y<0.0_r8) then
       lon = 0.0_r8
       lat = pi/2.0_r8+atan(y)
     elseif (x>0.0_r8.and.y==0.0_r8) then
       lon = pi/2.0_r8
       lat = pi/2.0_r8-atan(x)
     elseif (x<0.0_r8.and.y==0.0_r8) then
       lon = 3.0*pi/2.0_r8
       lat = pi/2.0_r8+atan(x)
     elseif (x==0.0_r8.and.y==0.0_r8) then
       lon = 0.0_r8
       lat = pi/2.0_r8
     elseif (x/= 0.0_r8.and.y>0.0_r8) then
       lon = pi-atan(x/y)
       lat = atan(sin(lon)/x)
     elseif (x/= 0.0_r8.and.y<0.0_r8) then
       lon = atan(-x/y)
       lat = atan(sin(lon)/x)
     endif
   else
     print*,'check iface in xy2lonlat....'
     stop
   endif
!
   if (lon<0.0_r8)    lon = lon+2.0_r8*pi
   if (lon>2.0_r8*pi) lon = lon-2.0_r8*pi
!
   return
   end subroutine xy2lonlat
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_print_coordinates(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: ne, iface, ie, je, i
!
   ne = cs%ne
!
   ! face,global index,catesian index
   do iface = 1,cs%nface
     print*,'iface = ',iface
     do je = 1,ne
     do ie = 1,ne
       write(*,'(a,x,i2,x,i2)') '(ie,je) = ',ie,je
       write(*,'(x,a,x,i3)')    '(rank)  = ',cs%elem(ie,je,iface)%rank
       write(*,'(x,a,x,i3)')    '(local) = ',cs%elem(ie,je,iface)%local
     enddo
     enddo
   enddo
!
   return
   end subroutine cs_print_coordinates
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_print_directions(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: ne, iface, ie, je, i
!
   ne = cs%ne
!
   ! face,global index,catesian index
   do iface = 1,cs%nface
     print*,'iface = ',iface
     call mesh_print(cs%elem(:,:,iface)%joiner,.true.)
   enddo
!
   return
   end subroutine cs_print_directions
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_sfc_index(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: ne, iface, ie, je, i
   integer(i4) :: base, iaxis, naxis
   integer(i4), dimension(3) :: axises
   type(sfc_t) :: sfc
   integer(i4), dimension(:,:), allocatable :: b_sfc, b_joiner, isfc, joiner, sfc_m, joiner_m
!
   ne = cs%ne
!
   call sfc_initialize(sfc,ne)
   call sfc_gen_curve(sfc)
   call sfc_get_sfc(sfc,b_sfc)
   allocate(b_joiner(ne,ne),isfc(ne,ne),joiner(ne,ne),sfc_m(ne,ne),joiner_m(ne,ne))
   b_joiner = sfc%curve%joiner
!
   do iface = 1,cs%nface
     !
     if (iface.eq.1) then
       base = 0
       naxis = 1
       axises(1) = axis_x
     elseif (iface.eq.2) then
       base = ne*ne
       naxis = 1
       axises(1) = axis_x
     elseif (iface.eq.3) then
       base = 5*ne*ne
       naxis = 1
       axises(1) = 0
     elseif (iface.eq.4) then
       base = 3*ne*ne
       naxis = 2
       axises(1) = axis_xy
       axises(2) = axis_x
     elseif (iface.eq.5) then
       base = 4*ne*ne
       naxis = 1
       axises(1) = 0
     elseif (iface.eq.6) then
       base = 2*ne*ne
       naxis = 1
       axises(1) = axis_0
     endif
     !
     sfc_m = b_sfc
     joiner_m = b_joiner
     do i = 1,naxis
       call mesh_convert(axises(i),sfc_m,isfc)
       call mesh_dir_convert(axises(i),joiner_m,joiner)
       sfc_m = isfc
       joiner_m = joiner
     enddo
     !
     do je = 1,ne
     do ie = 1,ne
       cs%elem(ie,je,iface)%isfc = base+isfc(ie,je)
       cs%elem(ie,je,iface)%joiner = joiner(ie,je)
     enddo
     enddo
!
   enddo
!
   call sfc_finalize(sfc)
   deallocate(b_sfc,b_joiner,isfc,joiner)
!
   return
   end subroutine set_sfc_index
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
! need optimize
   subroutine set_rank(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: ne, iface, ie, je, isfc
   integer(i4) :: nprocs, rank, ista, iend, dprocs
!
   nprocs = cs%nprocs
   ne = cs%ne
!
   do rank = 0,nprocs-1
     cs%nelemd(rank+1) = decompose1d(1,cs%nelem,nprocs,rank,ista,iend)
     do iface = 1,cs%nface
     do je = 1,ne
     do ie = 1,ne
       isfc = cs%elem(ie,je,iface)%isfc
       if (isfc.ge.ista.and.isfc.le.iend) then
         cs%elem(ie,je,iface)%rank = rank
       endif
     enddo
     enddo
     enddo
   enddo
!
   return
   end subroutine set_rank
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function decompose1d(n1, n2, nprocs, rank, ista, iend) result(nd)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: n1, n2
   integer(i4),           intent(in   ) :: nprocs, rank
   integer(i4), optional, intent(  out) :: ista, iend
   integer(i4) :: nd
! local variables
   integer(i4) :: l_sta,l_end,domain,extra
!
   domain = n2-n1+1
   nd = domain/nprocs
   extra = mod(domain,nprocs)
   l_sta = nd*rank+n1+min(rank,extra)
   l_end = l_sta+nd-1
   if (rank<=extra-1) then
     l_end = l_end+1
   endif
   nd = l_end-l_sta+1
!
   if (present(ista)) then
     ista = l_sta
   endif
   if (present(iend)) then
     iend = l_end
   endif
!
   return
!
   end function decompose1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_local_index(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: ne, iface, ie, je, nprocs, iproc
   integer(i4), dimension(:), allocatable :: local
!
   nprocs = cs%nprocs
   ne = cs%ne
   allocate(local(nprocs))
   local(:) = 0
!
   do iface = 1,cs%nface
   do je = 1,ne
   do ie = 1,ne
     iproc = cs%elem(ie,je,iface)%rank+1
     local(iproc) = local(iproc)+1
     cs%elem(ie,je,iface)%local = local(iproc)
   enddo
   enddo
   enddo
!
   deallocate(local)
!
   return
   end subroutine set_local_index
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_neighbors(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: ne, nface, nprocs
   integer(i4) :: iface, proc, i, j
   integer(i4) :: ww_f, ww_i, ww_j, ee_f, ee_i, ee_j, nn_f, nn_i, nn_j, ss_f, ss_i, ss_j
   integer(i4) :: ws_f, ws_i, ws_j, wn_f, wn_i, wn_j, es_f, es_i, es_j, en_f, en_i, en_j
!
   nface = cs%nface
   nprocs = cs%nprocs
   ne = cs%ne
!
   ! west, east, south, north
   do iface = 1,nface
   do i = 1,ne
   do j = 1,ne
     call get_neighbor_coord(iface,ne,i,j,ww,ww_f,ww_i,ww_j)
     call get_neighbor_coord(iface,ne,i,j,ee,ee_f,ee_i,ee_j)
     call get_neighbor_coord(iface,ne,i,j,ss,ss_f,ss_i,ss_j)
     call get_neighbor_coord(iface,ne,i,j,nn,nn_f,nn_i,nn_j)
     !
     ! ww
     cs%elem(i,j,iface)%nbrs(ww)%iface = ww_f
     cs%elem(i,j,iface)%nbrs(ww)%ie = ww_i
     cs%elem(i,j,iface)%nbrs(ww)%je = ww_j
     cs%elem(i,j,iface)%nbrs(ww)%global = cs%elem(ww_i,ww_j,ww_f)%global
     cs%elem(i,j,iface)%nbrs(ww)%rank = cs%elem(ww_i,ww_j,ww_f)%rank
     ! ee
     cs%elem(i,j,iface)%nbrs(ee)%iface = ee_f
     cs%elem(i,j,iface)%nbrs(ee)%ie = ee_i
     cs%elem(i,j,iface)%nbrs(ee)%je = ee_j
     cs%elem(i,j,iface)%nbrs(ee)%global = cs%elem(ee_i,ee_j,ee_f)%global
     cs%elem(i,j,iface)%nbrs(ee)%rank = cs%elem(ee_i,ee_j,ee_f)%rank
     ! ss
     cs%elem(i,j,iface)%nbrs(ss)%iface = ss_f
     cs%elem(i,j,iface)%nbrs(ss)%ie = ss_i
     cs%elem(i,j,iface)%nbrs(ss)%je = ss_j
     cs%elem(i,j,iface)%nbrs(ss)%global = cs%elem(ss_i,ss_j,ss_f)%global
     cs%elem(i,j,iface)%nbrs(ss)%rank = cs%elem(ss_i,ss_j,ss_f)%rank
     ! nn
     cs%elem(i,j,iface)%nbrs(nn)%iface = nn_f
     cs%elem(i,j,iface)%nbrs(nn)%ie = nn_i
     cs%elem(i,j,iface)%nbrs(nn)%je = nn_j
     cs%elem(i,j,iface)%nbrs(nn)%global = cs%elem(nn_i,nn_j,nn_f)%global
     cs%elem(i,j,iface)%nbrs(nn)%rank = cs%elem(nn_i,nn_j,nn_f)%rank
   enddo
   enddo
   enddo
!
   ! west-south,east-south,west-north,east-north
   do iface = 1,nface
   do i = 1,ne
   do j = 1,ne
!
     call get_neighbor_coord(iface,ne,i,j,ws,ws_f,ws_i,ws_j)
     call get_neighbor_coord(iface,ne,i,j,wn,wn_f,wn_i,wn_j)
     call get_neighbor_coord(iface,ne,i,j,es,es_f,es_i,es_j)
     call get_neighbor_coord(iface,ne,i,j,en,en_f,en_i,en_j)
     ! set ws
     cs%elem(i,j,iface)%nbrs(ws)%iface = ws_f
     cs%elem(i,j,iface)%nbrs(ws)%ie = ws_i
     cs%elem(i,j,iface)%nbrs(ws)%je = ws_j
     if (ws_f.ne.-1) then
       cs%elem(i,j,iface)%nbrs(ws)%global = cs%elem(ws_i,ws_j,ws_f)%global
       cs%elem(i,j,iface)%nbrs(ws)%rank = cs%elem(ws_i,ws_j,ws_f)%rank
     else
       cs%elem(i,j,iface)%nbrs(ws)%global =-1
       cs%elem(i,j,iface)%nbrs(ws)%rank =-1
     endif
     ! set es
     cs%elem(i,j,iface)%nbrs(es)%iface = es_f
     cs%elem(i,j,iface)%nbrs(es)%ie = es_i
     cs%elem(i,j,iface)%nbrs(es)%je = es_j
     if (es_f.ne.-1) then
       cs%elem(i,j,iface)%nbrs(es)%global = cs%elem(es_i,es_j,es_f)%global
       cs%elem(i,j,iface)%nbrs(es)%rank = cs%elem(es_i,es_j,es_f)%rank
     else
       cs%elem(i,j,iface)%nbrs(es)%global =-1
       cs%elem(i,j,iface)%nbrs(es)%rank =-1
     endif
     ! set wn
     cs%elem(i,j,iface)%nbrs(wn)%iface = wn_f
     cs%elem(i,j,iface)%nbrs(wn)%ie = wn_i
     cs%elem(i,j,iface)%nbrs(wn)%je = wn_j
     if (wn_f.ne.-1) then
       cs%elem(i,j,iface)%nbrs(wn)%global = cs%elem(wn_i,wn_j,wn_f)%global
       cs%elem(i,j,iface)%nbrs(wn)%rank = cs%elem(wn_i,wn_j,wn_f)%rank
     else
       cs%elem(i,j,iface)%nbrs(wn)%global =-1
       cs%elem(i,j,iface)%nbrs(wn)%rank =-1
     endif
     ! set en
     cs%elem(i,j,iface)%nbrs(en)%iface = en_f
     cs%elem(i,j,iface)%nbrs(en)%ie = en_i
     cs%elem(i,j,iface)%nbrs(en)%je = en_j
     if (en_f.ne.-1) then
       cs%elem(i,j,iface)%nbrs(en)%global = cs%elem(en_i,en_j,en_f)%global
       cs%elem(i,j,iface)%nbrs(en)%rank = cs%elem(en_i,en_j,en_f)%rank
     else
       cs%elem(i,j,iface)%nbrs(en)%global =-1
       cs%elem(i,j,iface)%nbrs(en)%rank =-1
     endif
!
   enddo
   enddo
   enddo
!
   return
   end subroutine set_neighbors
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_quadrature_ups_eps(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: np, ne, nface, nprocs
   integer(i4) :: iface, proc, ie, je, ip, jp
   logical(l4), dimension(4) :: isuniques, hasuniques
   integer(i4) :: wif, eif, sif, nif
   integer(i4) :: wie, eie, sie, nie
   integer(i4) :: wje, eje, sje, nje
   integer(i4) :: wip, eip, sip, nip
   integer(i4) :: wjp, ejp, sjp, njp
!
   nface = cs%nface
   ne = cs%ne
   np = cs%np
!
! set unique points
   cs%nups = 0
   cs%neps = 0
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
   do jp = 1,np
   do ip = 1,np
!
     call get_neighbor_quadrature(iface,ne,np,ie,je,ip,jp,ww,wif,wie,wje,wip,wjp,.true.) ! west
     call get_neighbor_quadrature(iface,ne,np,ie,je,ip,jp,ee,eif,eie,eje,eip,ejp,.true.) ! east
     call get_neighbor_quadrature(iface,ne,np,ie,je,ip,jp,ss,sif,sie,sje,sip,sjp,.true.) ! south
     call get_neighbor_quadrature(iface,ne,np,ie,je,ip,jp,nn,nif,nie,nje,nip,njp,.true.) ! north
     !
     isuniques(1) = cs%elem(wie,wje,wif)%qp(wip,wjp)%isunique
     isuniques(2) = cs%elem(eie,eje,eif)%qp(eip,ejp)%isunique
     isuniques(3) = cs%elem(sie,sje,sif)%qp(sip,sjp)%isunique
     isuniques(4) = cs%elem(nie,nje,nif)%qp(nip,njp)%isunique
     hasuniques(1) = cs%elem(wie,wje,wif)%qp(wip,wjp)%hasunique
     hasuniques(2) = cs%elem(eie,eje,eif)%qp(eip,ejp)%hasunique
     hasuniques(3) = cs%elem(sie,sje,sif)%qp(sip,sjp)%hasunique
     hasuniques(4) = cs%elem(nie,nje,nif)%qp(nip,njp)%hasunique
     !
     if (.not.hasuniques(1).and..not.hasuniques(2).and..not.hasuniques(3).and..not.hasuniques(4)) then
       cs%elem(ie,je,iface)%qp(ip,jp)%isunique  = .true.
       cs%elem(ie,je,iface)%qp(ip,jp)%hasunique = .true.
       cs%nups = cs%nups+1
       cs%elem(ie,je,iface)%qp(ip,jp)%iup    = cs%nups
       cs%elem(ie,je,iface)%qp(ip,jp)%u_face = iface
       cs%elem(ie,je,iface)%qp(ip,jp)%u_ie   = ie
       cs%elem(ie,je,iface)%qp(ip,jp)%u_je   = je
       cs%elem(ie,je,iface)%qp(ip,jp)%u_ip   = ip
       cs%elem(ie,je,iface)%qp(ip,jp)%u_jp   = jp
     elseif (hasuniques(1).or.hasuniques(2).or.hasuniques(3).or.hasuniques(4)) then
       cs%elem(ie,je,iface)%qp(ip,jp)%isunique  = .false.
       cs%elem(ie,je,iface)%qp(ip,jp)%hasunique = .true.
       if (hasuniques(1)) then
         cs%elem(ie,je,iface)%qp(ip,jp)%iup    = cs%elem(wie,wje,wif)%qp(wip,wjp)%iup
         cs%elem(ie,je,iface)%qp(ip,jp)%u_face = cs%elem(wie,wje,wif)%qp(wip,wjp)%u_face
         cs%elem(ie,je,iface)%qp(ip,jp)%u_ie   = cs%elem(wie,wje,wif)%qp(wip,wjp)%u_ie
         cs%elem(ie,je,iface)%qp(ip,jp)%u_je   = cs%elem(wie,wje,wif)%qp(wip,wjp)%u_je
         cs%elem(ie,je,iface)%qp(ip,jp)%u_ip   = cs%elem(wie,wje,wif)%qp(wip,wjp)%u_ip
         cs%elem(ie,je,iface)%qp(ip,jp)%u_jp   = cs%elem(wie,wje,wif)%qp(wip,wjp)%u_jp
       elseif (hasuniques(2)) then
         cs%elem(ie,je,iface)%qp(ip,jp)%iup    = cs%elem(eie,eje,eif)%qp(eip,ejp)%iup
         cs%elem(ie,je,iface)%qp(ip,jp)%u_face = cs%elem(eie,eje,eif)%qp(eip,ejp)%u_face
         cs%elem(ie,je,iface)%qp(ip,jp)%u_ie   = cs%elem(eie,eje,eif)%qp(eip,ejp)%u_ie
         cs%elem(ie,je,iface)%qp(ip,jp)%u_je   = cs%elem(eie,eje,eif)%qp(eip,ejp)%u_je
         cs%elem(ie,je,iface)%qp(ip,jp)%u_ip   = cs%elem(eie,eje,eif)%qp(eip,ejp)%u_ip
         cs%elem(ie,je,iface)%qp(ip,jp)%u_jp   = cs%elem(eie,eje,eif)%qp(eip,ejp)%u_jp
       elseif (hasuniques(3)) then
         cs%elem(ie,je,iface)%qp(ip,jp)%iup    = cs%elem(sie,sje,sif)%qp(sip,sjp)%iup
         cs%elem(ie,je,iface)%qp(ip,jp)%u_face = cs%elem(sie,sje,sif)%qp(sip,sjp)%u_face
         cs%elem(ie,je,iface)%qp(ip,jp)%u_ie   = cs%elem(sie,sje,sif)%qp(sip,sjp)%u_ie
         cs%elem(ie,je,iface)%qp(ip,jp)%u_je   = cs%elem(sie,sje,sif)%qp(sip,sjp)%u_je
         cs%elem(ie,je,iface)%qp(ip,jp)%u_ip   = cs%elem(sie,sje,sif)%qp(sip,sjp)%u_ip
         cs%elem(ie,je,iface)%qp(ip,jp)%u_jp   = cs%elem(sie,sje,sif)%qp(sip,sjp)%u_jp
       elseif (hasuniques(4)) then
         cs%elem(ie,je,iface)%qp(ip,jp)%iup    = cs%elem(nie,nje,nif)%qp(nip,njp)%iup
         cs%elem(ie,je,iface)%qp(ip,jp)%u_face = cs%elem(nie,nje,nif)%qp(nip,njp)%u_face
         cs%elem(ie,je,iface)%qp(ip,jp)%u_ie   = cs%elem(nie,nje,nif)%qp(nip,njp)%u_ie
         cs%elem(ie,je,iface)%qp(ip,jp)%u_je   = cs%elem(nie,nje,nif)%qp(nip,njp)%u_je
         cs%elem(ie,je,iface)%qp(ip,jp)%u_ip   = cs%elem(nie,nje,nif)%qp(nip,njp)%u_ip
         cs%elem(ie,je,iface)%qp(ip,jp)%u_jp   = cs%elem(nie,nje,nif)%qp(nip,njp)%u_jp
       endif
     else
       cs%elem(ie,je,iface)%qp(ip,jp)%isunique  = .false.
       cs%elem(ie,je,iface)%qp(ip,jp)%hasunique = .false.
       cs%elem(ie,je,iface)%qp(ip,jp)%iup       = -1
       cs%elem(ie,je,iface)%qp(ip,jp)%u_face    = -1
       cs%elem(ie,je,iface)%qp(ip,jp)%u_ie      = -1
       cs%elem(ie,je,iface)%qp(ip,jp)%u_je      = -1
       cs%elem(ie,je,iface)%qp(ip,jp)%u_ip      = -1
       cs%elem(ie,je,iface)%qp(ip,jp)%u_jp      = -1
       print*,'WHAT!!!'
       stop
     endif
!
     cs%neps = cs%neps+1
     if (cs%elem(ie,je,iface)%qp(ip,jp)%iup==33423) then
       print*, cs%elem(ie,je,iface)%qp(ip,jp)%iup,iface,ie,je,ip,jp
     endif
!
   enddo
   enddo
   enddo
   enddo
   enddo
!
   print*,'nEPs = ',cs%neps
   print*,'nUPs = ',cs%nups
!
   return
   end subroutine set_quadrature_ups_eps
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_quadrature_coordinate(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: np, ne, nface, nprocs
   integer(i4) :: iface, proc, ie, je, ip, jp
   real(r8), dimension(4) :: a, b ! alpha, beta ! for element vetices
!   real(r8), dimension(cs%np) :: a, b ! alpha, beta ! for element vetices
   real(r8), dimension(cs%np) :: u, v ! reference domain ! for quadrature points
   type(bilinear_t) :: a_t, b_t ! bilinear interpolation for alpha, beta
   integer(i4) :: iup, u_face, u_ie, u_je, u_ip, u_jp
   real(r8) :: error
!
   nface = cs%nface
   ne = cs%ne
   np = cs%np
! set coordinate for quadrature points
   u = cs%qp
   v = cs%qp
!
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
     a(:) = cs%elem(ie,je,iface)%alpha(:)
     b(:) = cs%elem(ie,je,iface)%beta(:)
     call bilinear_set(a_t,a(:))
     call bilinear_set(b_t,b(:))
     do jp = 1,np
     do ip = 1,np
       call bilinear_interpolation(a_t,cs%elem(ie,je,iface)%qp(ip,jp)%alpha,u(ip),v(jp))
       call bilinear_interpolation(b_t,cs%elem(ie,je,iface)%qp(ip,jp)%beta,u(ip),v(jp))
       call alphabeta2lonlat(iface,cs%elem(ie,je,iface)%qp(ip,jp)%alpha,cs%elem(ie,je,iface)%qp(ip,jp)%beta,&
                                   cs%elem(ie,je,iface)%qp(ip,jp)%lon,  cs%elem(ie,je,iface)%qp(ip,jp)%lat, cs%isrotated)
   enddo
   enddo
   enddo
   enddo
   enddo
!
! check unique points
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
   do jp = 1,np
   do ip = 1,np
     if (.not.cs%elem(ie,je,iface)%qp(ip,jp)%isunique) then
       u_face = cs%elem(ie,je,iface)%qp(ip,jp)%u_face
       u_ie = cs%elem(ie,je,iface)%qp(ip,jp)%u_ie
       u_je = cs%elem(ie,je,iface)%qp(ip,jp)%u_je
       u_ip = cs%elem(ie,je,iface)%qp(ip,jp)%u_ip
       u_jp = cs%elem(ie,je,iface)%qp(ip,jp)%u_jp
       error =(cs%elem(ie,je,iface)%qp(ip,jp)%lon-cs%elem(u_ie,u_je,u_face)%qp(u_ip,u_jp)%lon)/cs%elem(ie,je,iface)%qp(ip,jp)%lon
       if (cs%elem(ie,je,iface)%qp(ip,jp)%lon.ne.cs%elem(u_ie,u_je,u_face)%qp(u_ip,u_jp)%lon.and. &
           abs(error).gt.zero_thr) then
         print*,'warning in set_quadrature_coordinate...'
         print*,iface,ie,je,ip,jp
         print*,u_face,u_ie,u_je,u_ip,u_jp
         print*,cs%elem(ie,je,iface)%qp(ip,jp)%lon,cs%elem(u_ie,u_je,u_face)%qp(u_ip,u_jp)%lon,error
         print*,' '
       endif
     endif
   enddo
   enddo
   enddo
   enddo
   enddo
!
   return
   end subroutine set_quadrature_coordinate
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_quadrature_coordinate_with_llfile(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: np, ne, nface, nprocs
   integer(i4) :: iface, proc, ie, je, ip, jp, iup
   ! netcdf
   integer(i4) :: ierr, file_id, dim_id, dimsize, lon_id, lat_id
   real(r8), dimension(:), allocatable :: lons, lats
   logical(l4) :: isdegree
!
   allocate(lons(cs%nups))
   allocate(lats(cs%nups))
!
   ierr = nf90_open(trim(cs%lonlatfilename),nf90_nowrite,file_id)
   if (ierr.ne.0) then
     print*,'Check file...',trim(cs%lonlatfilename)
     stop
   endif
! dimension
   ierr = nf90_inq_dimid(file_id,'ncol',dim_id)
   if (ierr.ne.0) then
     ierr = nf90_inq_dimid(file_id,'ncols',dim_id)
     if (ierr.ne.0) then
       ierr = nf90_inq_dimid(file_id,'grid_size',dim_id)
       if (ierr.ne.0) then
         ierr = nf90_inq_dimid(file_id,'nUPs',dim_id)
         if (ierr.ne.0) then
           print*,'Check dimension name... in ',trim(cs%lonlatfilename)
           stop
         endif
       endif
     endif
   endif
!
   ierr = nf90_inquire_dimension(file_id,dim_id,len=dimsize)
   if (ierr.ne.0.or.dimsize.ne.cs%nups) then
     print*,'Check dimsize.. ',dimsize,' in ',trim(cs%lonlatfilename)
     stop
   endif
! variable
   ierr = nf90_inq_varid(file_id,trim(cs%lonvarname),lon_id)
   if (ierr.ne.0) then
     print*,'Check lonname...',trim(cs%lonvarname)
     stop
   endif
   ierr = nf90_inq_varid(file_id,trim(cs%latvarname),lat_id)
   if (ierr.ne.0) then
     print*,'Check latname...',trim(cs%latvarname)
     stop
   endif
!
   ierr = nf90_get_var(file_id,lon_id,lons)
   ierr = nf90_get_var(file_id,lat_id,lats)
!
   if (maxval(lons).gt.2.0_r8*pi) then
     isdegree = .true.
   else
     isdegree = .false.
   endif
!
   nface = cs%nface
   ne = cs%ne
   np = cs%np
!
! set coordinate for quadrature points
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
   do jp = 1,np
   do ip = 1,np
     cs%elem(ie,je,iface)%qp(ip,jp)%alpha =-1.0_r8
     cs%elem(ie,je,iface)%qp(ip,jp)%beta =-1.0_r8
  
     iup = cs%elem(ie,je,iface)%qp(ip,jp)%iup
     if (isdegree) then
       cs%elem(ie,je,iface)%qp(ip,jp)%lon = lons(iup)*pi/180.0_r8
       cs%elem(ie,je,iface)%qp(ip,jp)%lat = lats(iup)*pi/180.0_r8
     else
       cs%elem(ie,je,iface)%qp(ip,jp)%lon = lons(iup)
       cs%elem(ie,je,iface)%qp(ip,jp)%lat = lats(iup)
     endif
   enddo
   enddo
   enddo
   enddo
   enddo
!
   deallocate(lons)
   deallocate(lats)
!
   ierr = nf90_close(file_id)
!
   return
   end subroutine set_quadrature_coordinate_with_llfile
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_quadrature_coordinate_wo_alphabeta(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: np, ne, nface, nprocs
   integer(i4) :: iface, proc, ie, je, ip, jp, ii
   real(r8), dimension(4) :: u, v ! reference domain for quadrature points
   real(r8), dimension(4) :: a, b ! reference domain for quadrature points
   real(r8) :: r
   real(r8), dimension(4) :: rs, lons, lats, rrs, rlons, rlats, x, y, z ! lons, lats for element vetices
   real(r8) :: cr, clon, clat, cx, cy, cz ! center vector
   real(r8), dimension(cs%np, cs%np) :: alphaq, betaq, rq, lonq, latq, rrq, rlonq, rlatq
   type(bilinear_t) :: a_t, b_t
!
   nface = cs%nface
   ne = cs%ne
   np = cs%np
! set coordinate for quadrature points
   u = cs%qp
   v = cs%qp
!
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
     r = 1.0_r8
     rs(:) = 1.0_r8
     rrs(:) = 1.0_r8
     !
     lons(:) = cs%elem(ie,je,iface)%lon(:)
     lats(:) = cs%elem(ie,je,iface)%lat(:)
     do ii = 1,4
       call lonlat2cartesian(r,lons(ii),lats(ii),x(ii),y(ii),z(ii))
     enddo
     cx = 0.0_r8
     cy = 0.0_r8
     cz = 0.0_r8
     do ii = 1,4
       cx = cx+x(ii)
       cy = cy+y(ii)
       cz = cz+z(ii)
     enddo
     call cartesian2lonlat(cx,cy,cz,cr,clon,clat)
     ! lon,lat -> arbitrary alpha,beta
     call rotation_lonslats_ll0_to_axis(clon,clat,'x',4,rs,lons,lats,rrs,rlons,rlats)
     do ii = 1,4
       call lonlat_to_great_circle_angle(rlons(ii),rlats(ii),a(ii),b(ii))
     enddo
     call bilinear_set(a_t,a(:))
     call bilinear_set(b_t,b(:))
     !
     do jp = 1,np
     do ip = 1,np
       call bilinear_interpolation(a_t,alphaq(ip,jp),u(ip),v(jp))
       call bilinear_interpolation(b_t,betaq(ip,jp),u(ip),v(jp))
       call great_circle_angle_to_lonlat(alphaq(ip,jp),betaq(ip,jp),lonq(ip,jp),latq(ip,jp))
     enddo
     enddo
     rq(:,:) = 1.0_r8
     call inv_rotation_lonslats_ll0_to_axis(clon,clat,'x',np*np,rq,lonq,latq,rrq,rlonq,rlatq)
     !
     do jp = 1,np
     do ip = 1,np
       cs%elem(ie,je,iface)%qp(ip,jp)%alpha = alphaq(ip,jp)
       cs%elem(ie,je,iface)%qp(ip,jp)%beta = betaq(ip,jp)
       cs%elem(ie,je,iface)%qp(ip,jp)%lon = rlonq(ip,jp)
       cs%elem(ie,je,iface)%qp(ip,jp)%lat = rlatq(ip,jp)
     enddo
     enddo
!
   enddo
   enddo
   enddo
!
   return
   end subroutine set_quadrature_coordinate_wo_alphabeta
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine convert_up2ep(cs, ups, eps)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t),                                  intent(inout) :: cs
   real(r8), dimension(:),                                intent(in   ) :: ups
   real(r8), dimension(cs%np,cs%np,cs%ne,cs%ne,cs%nface), intent(  out) :: eps
! local variables
   integer(i4) :: nsize, iup
   integer(i4) :: np, ne, nface, nprocs
   integer(i4) :: iface, ie, je, ip, jp
!
   nface = cs%nface
   ne = cs%ne
   np = cs%np
!
   if (size(ups).ne.cs%nups) then
     print*,'Check up size in convert_up2ep...'
     stop
   endif
   if (size(eps).ne.cs%neps) then
     print*,'Check ep size in convert_up2ep...'
     stop
   endif
!
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
   do jp = 1,np
   do ip = 1,np
     iup = cs%elem(ie,je,iface)%qp(ip,jp)%iup
     eps(ip,jp,ie,je,iface) = ups(iup)
   enddo
   enddo
   enddo
   enddo
   enddo
!
   return
   end subroutine convert_up2ep
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_quadrature_neighbors(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: np, ne, nface!, nnbrs
   integer(i4) :: iface, ie, je, ip, jp, id, nbr
   logical(l4), dimension(nnbrs) :: isuniques
   integer(i4), dimension(nnbrs) :: dirs
   integer(i4) :: face, iie, jje, iip, jjp
!
   nface = cs%nface
   ne = cs%ne
   np = cs%np
!
   dirs =(/ww,ws,ss,es,ee,en,nn,wn/)
! set neighbors
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
   do jp = 1,np
   do ip = 1,np
!
     do nbr = 1,nnbrs
       call get_neighbor_quadrature(iface,ne,np,ie,je,ip,jp,nbr,face,iie,jje,iip,jjp)
       if (face.ge.1) then
         cs%elem(ie,je,iface)%qp(ip,jp)%nbrs(nbr)%lon = cs%elem(iie,jje,face)%qp(iip,jjp)%lon
         cs%elem(ie,je,iface)%qp(ip,jp)%nbrs(nbr)%lat = cs%elem(iie,jje,face)%qp(iip,jjp)%lat
       else
         cs%elem(ie,je,iface)%qp(ip,jp)%nbrs(nbr)%lon =-9999.0_r8
         cs%elem(ie,je,iface)%qp(ip,jp)%nbrs(nbr)%lat =-9999.0_r8
       endif
     enddo
!
   enddo
   enddo
   enddo
   enddo
   enddo
!
   return
   end subroutine set_quadrature_neighbors
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_quadrature_cells(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
! local variables
   integer(i4) :: np, ne, nface!, nnbrs, ncell
   integer(i4) :: iface, ie, je, ip, jp, d, ii, ic
   integer(i4), dimension(nvrts, ncell) :: dirs
   real(r8) :: r
   real(r8), dimension(nvrts) :: lons, lats
   real(r8), dimension(nvrts) :: x, y, z
   real(r8) :: sx, sy, sz, sr, slon, slat
   logical(l4) :: iscell
   integer(i4) :: ncs
!
   nface = cs%nface
   ne = cs%ne
   np = cs%np
!
   dirs(:,1) =(/-1,ww,ws,ss/) ! ws
   dirs(:,2) =(/-1,ss,es,ee/) ! es
   dirs(:,3) =(/-1,ee,en,nn/) ! en
   dirs(:,4) =(/-1,nn,wn,ww/) ! wn
!
   r = 1.0_r8
!
! set cells
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
   do jp = 1,np
   do ip = 1,np
!
     ncs = 0
     lons(1) = cs%elem(ie,je,iface)%qp(ip,jp)%lon
     lats(1) = cs%elem(ie,je,iface)%qp(ip,jp)%lat
  
     do ic = 1,ncell
       ii = 1
       iscell = .true.
       !
       do d = 2,4
         if (cs%elem(ie,je,iface)%qp(ip,jp)%nbrs(dirs(d,ic))%lon.ne.-9999.0_r8) then
           ii = ii+1
           lons(ii) = cs%elem(ie,je,iface)%qp(ip,jp)%nbrs(dirs(d,ic))%lon
           lats(ii) = cs%elem(ie,je,iface)%qp(ip,jp)%nbrs(dirs(d,ic))%lat
         else
           iscell = .false.
         endif
       enddo
       !
       if (iscell) then
       !
         do d = 1,4
           call lonlat2cartesian(r,lons(d),lats(d),x(d),y(d),z(d))
         enddo
         sx = sum(x) ; sy = sum(y) ; sz = sum(z)
         sr = sqrt(sx*sx+sy*sy+sz*sz)
         sx = r/sr*sx; sy = r/sr*sy; sz = r/sr*sz
         call cartesian2lonlat(sx,sy,sz,sr,slon,slat)
         if (abs(r-sr)>1.0d-14) then
           print*,'check sum lon,lat...',r,sr
           stop
         endif
         !
         ncs = ncs+1
         cs%elem(ie,je,iface)%qp(ip,jp)%cell(ncs)%lon = slon
         cs%elem(ie,je,iface)%qp(ip,jp)%cell(ncs)%lat = slat
       endif
     enddo
     ! checking ncell
     if (ncs.ne.ncell) then
       if (ncs.ne.ncell-1) then
         print*,'check ncells',ncs,iface,ie,je,ip,jp
         stop
       endif
       cs%elem(ie,je,iface)%qp(ip,jp)%cell(ncell)%lon = cs%elem(ie,je,iface)%qp(ip,jp)%cell(ncell-1)%lon
       cs%elem(ie,je,iface)%qp(ip,jp)%cell(ncell)%lat = cs%elem(ie,je,iface)%qp(ip,jp)%cell(ncell-1)%lat
     endif
!
   enddo
   enddo
   enddo
   enddo
   enddo
!
   return
   end subroutine set_quadrature_cells
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   recursive subroutine get_neighbor_coord(iface, n, i, j, direction, dst_f, dst_i, dst_j)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: iface, n, i, j, direction
   integer(i4), intent(  out) :: dst_f, dst_i, dst_j
! local variables
   integer(i4) :: tmp_f, tmp_i, tmp_j, dir1, dir2
!
   if (direction.eq.ww) then
!
     dst_f = iface
     dst_i = i-1
     dst_j = j
     !
     if (i.eq.1) then
       if (iface.eq.1) then
         dst_f = 4
         dst_i = n
       elseif (iface.eq.2) then
         dst_f = 1
         dst_i = n
       elseif (iface.eq.3) then
         dst_f = 2
         dst_i = n
       elseif (iface.eq.4) then
         dst_f = 3
         dst_i = n
       elseif (iface.eq.5) then
         dst_f = 4
         dst_i = j
         dst_j = 1
       elseif (iface.eq.6) then
         dst_f = 4
         dst_i = n-j+1
         dst_j = n
       endif
     endif
!
   elseif (direction.eq.ee) then
!
     dst_f = iface
     dst_i = i+1
     dst_j = j
     if (i.eq.n) then
       if (iface.eq.1) then
         dst_f = 2
         dst_i = 1
       elseif (iface.eq.2) then
         dst_f = 3
         dst_i = 1
       elseif (iface.eq.3) then
         dst_f = 4
         dst_i = 1
       elseif (iface.eq.4) then
         dst_f = 1
         dst_i = 1
       elseif (iface.eq.5) then
         dst_f = 2
         dst_i = n-j+1
         dst_j = 1
       elseif (iface.eq.6) then
         dst_f = 2
         dst_i = j
         dst_j = n
       endif
     endif
!
   elseif (direction.eq.ss) then
!
     dst_f = iface
     dst_i = i
     dst_j = j-1
     if (j.eq.1) then
       if (iface.eq.1) then
         dst_f = 5
         dst_j = n
       elseif (iface.eq.2) then
         dst_f = 5
         dst_i = n
         dst_j = n-i+1
       elseif (iface.eq.3) then
         dst_f = 5
         dst_i = n-i+1
         dst_j = 1
       elseif (iface.eq.4) then
         dst_f = 5
         dst_i = 1
         dst_j = i
       elseif (iface.eq.5) then
         dst_f = 3
         !dst_j = n! ?
         !dst_j = 1
         dst_i = n-i+1
         dst_j = 1
       elseif (iface.eq.6) then
         dst_f = 1
         dst_j = n
       endif
     endif
!
   elseif (direction.eq.nn) then
!
     dst_f = iface
     dst_i = i
     dst_j = j+1
     if (j.eq.n) then
       if (iface.eq.1) then
         dst_f = 6
         dst_j = 1
       elseif (iface.eq.2) then
         dst_f = 6
         dst_i = n
         dst_j = i
       elseif (iface.eq.3) then
         dst_f = 6
         dst_i = n-i+1
         dst_j = n
       elseif (iface.eq.4) then
         dst_f = 6
         dst_i = 1
         dst_j = n-i+1
       elseif (iface.eq.5) then
         dst_f = 1
         dst_j = 1
       elseif (iface.eq.6) then
         !dst_f = 5   !??????????
         !dst_j = 1   !??????????
         dst_f = 3
         dst_i = n-i+1
         dst_j = n
       endif
     endif
!
   elseif (direction.eq.ws) then
!
     if (i==1.and.j==1) then
       dst_f =-1
       dst_i =-1
       dst_j =-1
       return
     endif
     !
     if ((iface==5).or.(iface==6).and.(i==1)) then
       dir1 = ss
       dir2 = ww
     else
       dir1 = ww
       dir2 = ss
     endif
     call get_neighbor_coord(iface,n,i,j,dir1,tmp_f,tmp_i,tmp_j)
     call get_neighbor_coord(tmp_f,n,tmp_i,tmp_j,dir2,dst_f,dst_i,dst_j)
!
   elseif (direction.eq.wn) then
!
     if (i==1.and.j==n) then
       dst_f =-1
       dst_i =-1
       dst_j =-1
       return
     endif
     !
     if ((iface==5).or.(iface==6).and.(i==1)) then
       dir1 = nn
       dir2 = ww
     else
       dir1 = ww
       dir2 = nn
     endif
     call get_neighbor_coord(iface,n,i,j,dir1,tmp_f,tmp_i,tmp_j)
     call get_neighbor_coord(tmp_f,n,tmp_i,tmp_j,dir2,dst_f,dst_i,dst_j)
!
   elseif (direction.eq.es) then
!
     if (i==n.and.j==1) then
       dst_f =-1
       dst_i =-1
       dst_j =-1
       return
     endif
     !
     if ((iface==5).or.(iface==6).and.(i==n)) then
       dir1 = ss
       dir2 = ee
     else
       dir1 = ee
       dir2 = ss
     endif
     call get_neighbor_coord(iface,n,i,j,dir1,tmp_f,tmp_i,tmp_j)
     call get_neighbor_coord(tmp_f,n,tmp_i,tmp_j,dir2,dst_f,dst_i,dst_j)
!
   elseif (direction.eq.en) then
!
     if (i==n.and.j==n) then
       dst_f =-1
       dst_i =-1
       dst_j =-1
       return
     endif
     !
     if ((iface==5).or.(iface==6).and.(i==n)) then
       dir1 = nn
       dir2 = ee
     else
       dir1 = ee
       dir2 = nn
     endif
     call get_neighbor_coord(iface,n,i,j,dir1,tmp_f,tmp_i,tmp_j)
     call get_neighbor_coord(tmp_f,n,tmp_i,tmp_j,dir2,dst_f,dst_i,dst_j)
!
   else
     print*,'not supported direction...(get_neighbor_coord) '
     stop
   endif
!
   return
   end subroutine get_neighbor_coord
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   recursive subroutine get_neighbor_quadrature(iface, ne, np, ie, je, ip, jp, direction, dst_f, dst_ie, dst_je, dst_ip, dst_jp, chkup)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: iface, ne, np, ie, je, ip, jp, direction
   integer(i4),           intent(  out) :: dst_f, dst_ie, dst_je, dst_ip, dst_jp
   logical(l4), optional, intent(in   ) :: chkup
! local variables
   integer(i4) :: tmp_f, tmp_ie, tmp_je, tmp_ip, tmp_jp
   integer(i4) :: dir1, dir2
   logical(l4) :: l_chkup
!
   if (present(chkup)) then
     l_chkup = chkup
   else
     l_chkup = .false.
   endif
!--------------
! inner
   dst_f = iface
   dst_ie = ie
   dst_je = je
   if (l_chkup) then ! unique point
     dst_ip = ip
     dst_jp = jp
   else ! neighbor
     if (direction==ww) then
       dst_ip = ip-1
       dst_jp = jp
     elseif (direction==ee) then
       dst_ip = ip+1
       dst_jp = jp
     elseif (direction==ss) then
       dst_ip = ip
       dst_jp = jp-1
     elseif (direction==nn) then
       dst_ip = ip
       dst_jp = jp+1
     elseif (direction==ws) then
       dst_ip = ip-1
       dst_jp = jp-1
     elseif (direction==wn) then
       dst_ip = ip-1
       dst_jp = jp+1
     elseif (direction==es) then
       dst_ip = ip+1
       dst_jp = jp-1
     elseif (direction==en) then
       dst_ip = ip+1
       dst_jp = jp+1
     else
       print*,'check direction...'
       stop
     endif
   endif
!
!--------------
! outer
!  - west
   if (direction==ww.and.ip.eq.1) then
!
     if (ie.eq.1) then
       call get_neighbor_coord(iface,ne,ie,je,ww,dst_f,dst_ie,dst_je)
       call get_neighbor_coord(iface,np,ip,jp,ww,dst_f,dst_ip,dst_jp)
       if (.not.l_chkup) then
         if (iface.eq.5) then
           dst_jp = dst_jp+1
         elseif (iface.eq.6) then
           dst_jp = dst_jp-1
         else
           dst_ip = dst_ip-1
         endif
       endif
     else
       dst_f = iface
       dst_ie = ie-1
       dst_je = je
       dst_ip = np
       dst_jp = jp
       if (.not.l_chkup) then
         dst_ip = dst_ip-1
       endif
     endif
!
   endif
!
!  - east
   if (direction==ee.and.ip.eq.np) then
!
     if (ie.eq.ne) then
       call get_neighbor_coord(iface,ne,ie,je,ee,dst_f,dst_ie,dst_je)
       call get_neighbor_coord(iface,np,ip,jp,ee,dst_f,dst_ip,dst_jp)
       if (.not.l_chkup) then
         if (iface.eq.5) then
           dst_jp = dst_jp+1
         elseif (iface.eq.6) then
           dst_jp = dst_jp-1
         else
           dst_ip = dst_ip+1
         endif
       endif
     else
       dst_f = iface
       dst_ie = ie+1
       dst_je = je
       dst_ip = 1
       dst_jp = jp
       if (.not.l_chkup) then
         dst_ip = dst_ip+1
       endif
     endif
!
   endif
!
!  - south
   if (direction==ss.and.jp.eq.1) then
!
     if (je.eq.1) then
       call get_neighbor_coord(iface,ne,ie,je,ss,dst_f,dst_ie,dst_je)
       call get_neighbor_coord(iface,np,ip,jp,ss,dst_f,dst_ip,dst_jp)
       if (.not.l_chkup) then
         if (iface.eq.2) then
           dst_ip = dst_ip-1
         elseif (iface.eq.3) then
           dst_jp = dst_jp+1
         elseif (iface.eq.4) then
           dst_ip = dst_ip+1
         elseif (iface.eq.5) then
           dst_jp = dst_jp+1
         else
           dst_jp = dst_jp-1
         endif
       endif
     else
       dst_f = iface
       dst_ie = ie
       dst_je = je-1
       dst_ip = ip
       dst_jp = np
       if (.not.l_chkup) then
         dst_jp = dst_jp-1
       endif
     endif
!
   endif
!
!  - north
   if (direction==nn.and.jp.eq.np) then
!
     if (je.eq.ne) then
       call get_neighbor_coord(iface,ne,ie,je,nn,dst_f,dst_ie,dst_je)
       call get_neighbor_coord(iface,np,ip,jp,nn,dst_f,dst_ip,dst_jp)
       if (.not.l_chkup) then
         if (iface.eq.2) then
           dst_ip = dst_ip-1
         elseif (iface.eq.3) then
           dst_jp = dst_jp-1
         elseif (iface.eq.4) then
           dst_ip = dst_ip+1
         elseif (iface.eq.6) then
           dst_jp = dst_jp-1
         else
           dst_jp = dst_jp+1
         endif
       endif
     else
       dst_f = iface
       dst_ie = ie
       dst_je = je+1
       dst_ip = ip
       dst_jp = 1
       if (.not.l_chkup) then
         dst_jp = dst_jp+1
       endif
     endif
!
   endif
!
!  - west south
   if (direction==ws) then
!
     if ((ie.eq.1.and.je.eq.1).and.(ip.eq.1.and.jp.eq.1)) then
       dst_f =-1
       dst_ie =-1
       dst_je =-1
       dst_ip =-1
       dst_jp =-1
       return
     endif
     !
     dir1 = ww
     dir2 = ss
     call get_neighbor_quadrature(iface,ne,np,ie,je,ip,jp,dir1,tmp_f,tmp_ie,tmp_je,tmp_ip,tmp_jp,l_chkup)
     call get_neighbor_quadrature(tmp_f,ne,np,tmp_ie,tmp_je,tmp_ip,tmp_jp,dir2,dst_f,dst_ie,dst_je,dst_ip,dst_jp,l_chkup)
!
   endif
!
!  - east south
   if (direction==es.and.(ip.eq.np.or.jp.eq.1)) then
!
     if ((ie.eq.ne.and.je.eq.1).and.(ip.eq.np.and.jp.eq.1)) then
       dst_f =-1
       dst_ie =-1
       dst_je =-1
       dst_ip =-1
       dst_jp =-1
       return
     endif
  
     dir1 = ee
     dir2 = ss
     call get_neighbor_quadrature(iface,ne,np,ie,je,ip,jp,dir1,tmp_f,tmp_ie,tmp_je,tmp_ip,tmp_jp,l_chkup)
     call get_neighbor_quadrature(tmp_f,ne,np,tmp_ie,tmp_je,tmp_ip,tmp_jp,dir2,dst_f,dst_ie,dst_je,dst_ip,dst_jp,l_chkup)

!
   endif
!
!  - west north
   if (direction==wn.and.(ip.eq.1.or.jp.eq.np)) then
!
     if ((ie.eq.1.and.je.eq.ne).and.(ip.eq.1.and.jp.eq.np)) then
       dst_f =-1
       dst_ie =-1
       dst_je =-1
       dst_ip =-1
       dst_jp =-1
       return
     endif
  
     dir1 = ww
     dir2 = nn
     call get_neighbor_quadrature(iface,ne,np,ie,je,ip,jp,dir1,tmp_f,tmp_ie,tmp_je,tmp_ip,tmp_jp,l_chkup)
     call get_neighbor_quadrature(tmp_f,ne,np,tmp_ie,tmp_je,tmp_ip,tmp_jp,dir2,dst_f,dst_ie,dst_je,dst_ip,dst_jp,l_chkup)
!
   endif
!
!  - east north
   if (direction==en.and.(ip.eq.np.or.jp.eq.np)) then
!
     if ((ie.eq.ne.and.je.eq.ne).and.(ip.eq.np.and.jp.eq.np)) then
       dst_f =-1
       dst_ie =-1
       dst_je =-1
       dst_ip =-1
       dst_jp =-1
       return
     endif
  
     dir1 = ee
     dir2 = nn
     call get_neighbor_quadrature(iface,ne,np,ie,je,ip,jp,dir1,tmp_f,tmp_ie,tmp_je,tmp_ip,tmp_jp,l_chkup)
     call get_neighbor_quadrature(tmp_f,ne,np,tmp_ie,tmp_je,tmp_ip,tmp_jp,dir2,dst_f,dst_ie,dst_je,dst_ip,dst_jp,l_chkup)
!
   endif
!
   return
   end subroutine get_neighbor_quadrature
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_write_cube_infos(cs, filename)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(in   ) :: cs
   character(len=*),     intent(in   ) :: filename
! local variables
   integer(i4) :: ierr
   integer(i4) :: fileid
   integer(i4) :: did_nface, did_np, did_ne, did_nprocs, did_pnx, did_pny, did_nvrts, did_nnbrs, did_n_nbr
   integer(i4) :: vid_iface, vid_ie, vid_je, vid_global, vid_isfc, vid_rank, vid_local, vid_joiner, vid_nbrs
   integer(i4) :: vid_alpha, vid_beta, vid_lon, vid_lat, vid_iup, vid_isu
   integer(i4) :: np, ne, nface, nprocs, pnx, pny, iface, ie, je
   integer(i4), dimension(:,:), allocatable :: output_i4_2d
   integer(i4), dimension(:,:,:,:), allocatable :: output_i4_4d
   integer(i4), dimension(:,:,:,:), allocatable :: output_l4_4d
   real(r8), dimension(:,:), allocatable :: output_r8_2d
   real(r8), dimension(:,:,:), allocatable :: output_r8_3d
!
   np = cs%np
   ne = cs%ne
   nface = cs%nface
   nprocs = cs%nprocs
   pnx = cs%pnx
   pny = cs%pny
!
   ! file create
   ierr = nf90_create(trim(filename),ior(nf90_clobber,nf90_64bit_offset),fileid)
   ! attribute
   if (cs%isrotated) then
     ierr = nf90_put_att(fileid,nf90_global,'rotation',1)
   else
     ierr = nf90_put_att(fileid,nf90_global,'rotation',0)
   endif
   ! demensions
   ierr = nf90_def_dim(fileid,'nface',nface,did_nface)
   ierr = nf90_def_dim(fileid,'np',np,did_np)
   ierr = nf90_def_dim(fileid,'ne',ne,did_ne)
   ierr = nf90_def_dim(fileid,'nprocs',nprocs,did_nprocs)
   ierr = nf90_def_dim(fileid,'pnx',pnx,did_pnx)
   ierr = nf90_def_dim(fileid,'pny',pny,did_pny)
   ierr = nf90_def_dim(fileid,'nvrts',nvrts,did_nvrts)
   ierr = nf90_def_dim(fileid,'nnbrs',nnbrs,did_nnbrs)
   ierr = nf90_def_dim(fileid,'n_nbr',n_nbr,did_n_nbr)
   ! define variables (i4)
   ierr = nf90_def_var(fileid,'iface',nf90_real8,(/did_pnx,did_pny/),vid_iface)
   ierr = nf90_def_var(fileid,'ie',nf90_int,(/did_pnx,did_pny/),vid_ie)
   ierr = nf90_def_var(fileid,'je',nf90_int,(/did_pnx,did_pny/),vid_je)
   ierr = nf90_def_var(fileid,'global',nf90_int,(/did_pnx,did_pny/),vid_global)
   ierr = nf90_def_var(fileid,'isfc',nf90_int,(/did_pnx,did_pny/),vid_isfc)
   ierr = nf90_def_var(fileid,'rank',nf90_int,(/did_pnx,did_pny/),vid_rank)
   ierr = nf90_def_var(fileid,'local',nf90_int,(/did_pnx,did_pny/),vid_local)
   ierr = nf90_def_var(fileid,'joiner',nf90_int,(/did_pnx,did_pny/),vid_joiner)
   ierr = nf90_def_var(fileid,'nbrs',nf90_int,(/did_n_nbr,did_nnbrs,did_pnx,did_pny/),vid_nbrs)
   ierr = nf90_def_var(fileid,'iup',nf90_int,(/did_np,did_np,did_pnx,did_pny/),vid_iup)
   ierr = nf90_def_var(fileid,'isunique',nf90_int,(/did_np,did_np,did_pnx,did_pny/),vid_isu)
   ! define variables (r8)
   ierr = nf90_def_var(fileid,'alpha',nf90_real8,(/did_nvrts,did_pnx,did_pny/),vid_alpha)
   ierr = nf90_def_var(fileid,'beta',nf90_real8,(/did_nvrts,did_pnx,did_pny/),vid_beta)
   ierr = nf90_def_var(fileid,'lon',nf90_real8,(/did_nvrts,did_pnx,did_pny/),vid_lon)
   ierr = nf90_def_var(fileid,'lat',nf90_real8,(/did_nvrts,did_pnx,did_pny/),vid_lat)
   ! file close
   ierr = nf90_enddef(fileid)
!
   ! put variables (i4)
   call get_value_for_plot(cs,'iface',output_i4_2d)
   ierr = nf90_put_var(fileid,vid_iface,output_i4_2d)
   call get_value_for_plot(cs,'ie',output_i4_2d)
   ierr = nf90_put_var(fileid,vid_ie,output_i4_2d)
   call get_value_for_plot(cs,'je',output_i4_2d)
   ierr = nf90_put_var(fileid,vid_je,output_i4_2d)
   call get_value_for_plot(cs,'global',output_i4_2d)
   ierr = nf90_put_var(fileid,vid_global,output_i4_2d)
   call get_value_for_plot(cs,'isfc',output_i4_2d)
   ierr = nf90_put_var(fileid,vid_isfc,output_i4_2d)
   call get_value_for_plot(cs,'rank',output_i4_2d)
   ierr = nf90_put_var(fileid,vid_rank,output_i4_2d)
   call get_value_for_plot(cs,'local',output_i4_2d)
   ierr = nf90_put_var(fileid,vid_local,output_i4_2d)
   call get_value_for_plot(cs,'joiner',output_i4_2d)
   ierr = nf90_put_var(fileid,vid_joiner,output_i4_2d)
   call get_value_for_plot(cs,'nbrs',output_i4_4d)
   ierr = nf90_put_var(fileid,vid_nbrs,output_i4_4d)
   call get_value_for_plot(cs,output_i4_4d,output_l4_4d)
   ierr = nf90_put_var(fileid,vid_iup,output_i4_4d)
   ierr = nf90_put_var(fileid,vid_isu,output_l4_4d)
   ! put variables (r8)
   call get_value_for_plot(cs,'alpha',output_r8_3d)
   ierr = nf90_put_var(fileid,vid_alpha,output_r8_3d)
   call get_value_for_plot(cs,'beta',output_r8_3d)
   ierr = nf90_put_var(fileid,vid_beta,output_r8_3d)
   call get_value_for_plot(cs,'lon',output_r8_3d)
   ierr = nf90_put_var(fileid,vid_lon,output_r8_3d)
   call get_value_for_plot(cs,'lat',output_r8_3d)
   ierr = nf90_put_var(fileid,vid_lat,output_r8_3d)
!
   ! file close
   ierr = nf90_close(fileid)
!
   if (allocated(output_i4_2d)) deallocate(output_i4_2d)
   if (allocated(output_i4_4d)) deallocate(output_i4_4d)
   if (allocated(output_r8_2d)) deallocate(output_r8_2d)
   if (allocated(output_r8_3d)) deallocate(output_r8_3d)
!
   return
   end subroutine cs_write_cube_infos
! i4(2)   i4(2)     i4(3)         r8(3)           r8(3)        i4(2) i4(2)  i4(2)     l4(2)
! iface, global, ieje(x,y), ab(alpha,beta), lonlat(lon, lat), isfc, rank, local, isOut(logical)
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_value_for_plot_i4_2d(cs, vname, res)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t),                     intent(in   ) :: cs
   character(len=*),                         intent(in   ) :: vname
   integer(i4), dimension(:,:), allocatable, intent(  out) :: res
! local variables
   integer(i4) :: np, ne, nface, nprocs, pnx, pny, iface, ie, je, ic, jc
!
   np = cs%np
   ne = cs%ne
   nface = cs%nface
   nprocs = cs%nprocs
   pnx = cs%pnx
   pny = cs%pny
   ! value
   if (allocated(res)) deallocate(res)
   allocate(res(pnx,pny))
   res =-1
!
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
     call iface_ie_je_to_ix_iy(cs,iface,ie,je,ic,jc)
     !
     if (trim(vname).eq.'iface') then
       res(ic,jc) = cs%elem(ie,je,iface)%iface
     elseif (trim(vname).eq.'global') then
       res(ic,jc) = cs%elem(ie,je,iface)%global
     elseif (trim(vname).eq.'isfc') then
       res(ic,jc) = cs%elem(ie,je,iface)%isfc
     elseif (trim(vname).eq.'rank') then
       res(ic,jc) = cs%elem(ie,je,iface)%rank
     elseif (trim(vname).eq.'local') then
       res(ic,jc) = cs%elem(ie,je,iface)%local
     elseif (trim(vname).eq.'ie') then
       res(ic,jc) = cs%elem(ie,je,iface)%ie
     elseif (trim(vname).eq.'je') then
       res(ic,jc) = cs%elem(ie,je,iface)%je
     elseif (trim(vname).eq.'joiner') then
       res(ic,jc) = cs%elem(ie,je,iface)%joiner
     else
       print*,'check varname in get_value_for_plot_i4...'
       stop
     endif
   enddo
   enddo
   enddo
!
   return
   end subroutine get_value_for_plot_i4_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_value_for_plot_i4_4d(cs, vname, res)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t),                          intent(in   ) :: cs
   character(len=*),                              intent(in   ) :: vname
   integer(i4), dimension(:,:, :,:), allocatable, intent(  out) :: res
! local variables
   integer(i4) :: np, ne, nface, nprocs, pnx, pny, iface, ie, je, ic, jc, inbr
!
   np = cs%np
   ne = cs%ne
   nface = cs%nface
   nprocs = cs%nprocs
   pnx = cs%pnx
   pny = cs%pny
! value
   if (allocated(res)) deallocate(res)
   allocate(res(n_nbr,nnbrs,pnx,pny))
   res =-1
!
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
     call iface_ie_je_to_ix_iy(cs,iface,ie,je,ic,jc)
     if (trim(vname).eq.'nbrs') then
       do inbr = 1,nnbrs
       res(1,inbr,ic,jc) = cs%elem(ie,je,iface)%nbrs(inbr)%ie
       res(2,inbr,ic,jc) = cs%elem(ie,je,iface)%nbrs(inbr)%je
       res(3,inbr,ic,jc) = cs%elem(ie,je,iface)%nbrs(inbr)%global
       res(4,inbr,ic,jc) = cs%elem(ie,je,iface)%nbrs(inbr)%iface
       res(5,inbr,ic,jc) = cs%elem(ie,je,iface)%nbrs(inbr)%rank
       enddo
     else
       print*,'check varname in get_value_for_plot_i4_4d...'
       stop
     endif
   enddo
   enddo
   enddo
!
   return
   end subroutine get_value_for_plot_i4_4d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_value_for_plot_qp(cs, res, isu)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t),                         intent(in   ) :: cs
!   character(len=*),                             intent(in   ) :: vname
   integer(i4), dimension(:,:,:,:), allocatable, intent(  out) :: res
   integer(i4), dimension(:,:,:,:), allocatable, intent(  out) :: isu
! local variables
   integer(i4) :: np, ne, nface, nprocs, pnx, pny, iface, ip, jp, ie, je, ic, jc
!
   np     = cs%np
   ne     = cs%ne
   nface  = cs%nface
   nprocs = cs%nprocs
   pnx  = cs%pnx
   pny  = cs%pny
! value
   if (allocated(res)) deallocate(res)
   allocate(res(np,np,pnx,pny))
   if (allocated(isu)) deallocate(isu)
   allocate(isu(np,np,pnx,pny))
   res = -1
   isu = -1
!
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
     call iface_ie_je_to_ix_iy(cs,iface,ie,je,ic,jc)
     do jp = 1,np
     do ip = 1,np
       res(ip,jp,ic,jc) = cs%elem(ie,je,iface)%qp(ip,jp)%iup
       if (cs%elem(ie,je,iface)%qp(ip,jp)%isunique) then
         isu(ip,jp,ic,jc) = 1
       else
         isu(ip,jp,ic,jc) = 0
       endif
     enddo
     enddo
!     else
!       print*,'check varname in get_value_for_plot_r8...'
!       stop
!     endif
   enddo
   enddo
   enddo
!
   return
   end subroutine get_value_for_plot_qp
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_value_for_plot_r8_3d(cs, vname, res)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t),                    intent(in   ) :: cs
   character(len=*),                        intent(in   ) :: vname
   real(r8), dimension(:,:,:), allocatable, intent(  out) :: res
! local variables
   integer(i4) :: np, ne, nface, nprocs, pnx, pny, iface, ie, je, ic, jc
!
   np = cs%np
   ne = cs%ne
   nface = cs%nface
   nprocs = cs%nprocs
   pnx = cs%pnx
   pny = cs%pny
! value
   if (allocated(res)) deallocate(res)
   allocate(res(4,pnx,pny))
   res =-1
!
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
     call iface_ie_je_to_ix_iy(cs,iface,ie,je,ic,jc)
     if (trim(vname).eq.'alpha') then
       res(:,ic,jc) = cs%elem(ie,je,iface)%alpha
     elseif (trim(vname).eq.'beta') then
       res(:,ic,jc) = cs%elem(ie,je,iface)%beta
     elseif (trim(vname).eq.'lon') then
       res(:,ic,jc) = cs%elem(ie,je,iface)%lon
     elseif (trim(vname).eq.'lat') then
       res(:,ic,jc) = cs%elem(ie,je,iface)%lat
     else
       print*,'check varname in get_value_for_plot_r8...'
       stop
     endif
   enddo
   enddo
   enddo
!
   return
   end subroutine get_value_for_plot_r8_3d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ig_to_iface_ie_je(cs, ig, iface, ie, je)
   implicit none
   type(cubed_sphere_t), intent(in   ) :: cs
   integer(i4)         , intent(in   ) :: ig
   integer(i4)         , intent(  out) :: iface, ie, je
! local variables
   integer(i4) :: ne, idx_in_elem
!-------------------------------------------------------------------------------
!
   ne = cs%ne
   iface = (ig-1)/(ne*ne)+1
   idx_in_elem = mod(ig-1,ne*ne)
   ie = mod(idx_in_elem,ne)+1
   je = idx_in_elem/ne+1
!
   return
   end subroutine ig_to_iface_ie_je
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine iface_ie_je_to_ix_iy(cs, iface, ie, je, ix, iy)
   implicit none
   type(cubed_sphere_t), intent(in   ) :: cs
   integer(i4)         , intent(in   ) :: iface, ie, je
   integer(i4)         , intent(  out) :: ix, iy
! local variables
   integer(i4) :: ne, idx_in_elem
   integer(i4), dimension(2,6) :: base
!-------------------------------------------------------------------------------
!
   ne = cs%ne
   base(:,1) = (/ne,ne/)
   base(:,2) = (/2*ne,ne/)
   base(:,3) = (/3*ne,ne/)
   base(:,4) = (/0,ne/)
   base(:,5) = (/ne,0/)
   base(:,6) = (/ne,2*ne/)
!
   ix = base(1,iface)+ie
   iy = base(2,iface)+je
!
   return
   end subroutine iface_ie_je_to_ix_iy
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine create_plot_coordinate(cs, ifaces, ies, jes, global)
   implicit none
   type(cubed_sphere_t)                 , intent(in   ) :: cs
   integer(i4), dimension(cs%pnx,cs%pny), intent(  out) :: ifaces, ies, jes, global
! local variables
   integer(i4) :: nelem, ig, iface, ie, je, ix, iy
!-------------------------------------------------------------------------------
!
   nelem  = cs%nelem
   ifaces = -1
   ies    = -1
   jes    = -1
   global = -1
   do ig = 1,nelem
     call ig_to_iface_ie_je(cs,ig,iface,ie,je)
     call iface_ie_je_to_ix_iy(cs,iface,ie,je,ix,iy)
     ifaces(ix,iy) = iface
     ies   (ix,iy) = ie
     jes   (ix,iy) = je
     global(ix,iy) = ig
   enddo
!
   return
   end subroutine create_plot_coordinate
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine convert_global_elem_to_plot_0d_i4(cs, var,pvar)
   implicit none
   type(cubed_sphere_t)               , intent(in   ) :: cs
   integer(i4), dimension(cs%nelem)   , intent(in   ) :: var
   integer(i4), dimension(cs%pnx,cs%pny), intent(  out) :: pvar
! local variables
   integer(i4) :: nelem, ig, iface, ie, je, ix, iy
!-------------------------------------------------------------------------------
!
   nelem = cs%nelem
   pvar = -1
   do ig = 1,nelem
     call ig_to_iface_ie_je(cs,ig,iface,ie,je)
     call iface_ie_je_to_ix_iy(cs,iface,ie,je,ix,iy)
     pvar(ix,iy) = var(ig)
   enddo
!
   return
   end subroutine convert_global_elem_to_plot_0d_i4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine convert_global_elem_to_plot_0d_r8(cs, var, pvar)
   implicit none
   type(cubed_sphere_t)            , intent(in   ) :: cs
   real(r8), dimension(cs%nelem)   , intent(in   ) :: var
   real(r8), dimension(cs%pnx,cs%pny), intent(  out) :: pvar
! local variables
   integer(i4) :: nelem, ig, iface, ie, je, ix, iy
!-------------------------------------------------------------------------------
!
   nelem = cs%nelem
   pvar = -1.0_r8
   do ig = 1,nelem
     call ig_to_iface_ie_je(cs,ig,iface,ie,je)
     call iface_ie_je_to_ix_iy(cs,iface,ie,je,ix,iy)
     pvar(ix,iy) = var(ig)
   enddo
!
   return
   end subroutine convert_global_elem_to_plot_0d_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine convert_global_elem_to_plot_2d_i4(cs, var, pvar)
   implicit none
   type(cubed_sphere_t)                           , intent(in   ) :: cs
   integer(i4), dimension(cs%np,cs%np,cs%nelem)   , intent(in   ) :: var
   integer(i4), dimension(cs%np,cs%np,cs%pnx,cs%pny), intent(  out) :: pvar
! local variables
   integer(i4) :: nelem, ig, iface, ie, je, ix, iy
!-------------------------------------------------------------------------------
!
   nelem = cs%nelem
   pvar = -1
   do ig = 1,nelem
     call ig_to_iface_ie_je(cs,ig,iface,ie,je)
     call iface_ie_je_to_ix_iy(cs,iface,ie,je,ix,iy)
     pvar(:,:,ix,iy) = var(:,:,ig)
   enddo
!
   return
   end subroutine convert_global_elem_to_plot_2d_i4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine convert_global_elem_to_plot_2d_r8(cs,var,pvar)
   implicit none
   type(cubed_sphere_t)                        , intent(in   ) :: cs
   real(r8), dimension(cs%np,cs%np,cs%nelem)   , intent(in   ) :: var
   real(r8), dimension(cs%np,cs%np,cs%pnx,cs%pny), intent(  out) :: pvar
! local variables
   integer(i4) :: nelem, ig, iface, ie, je, ix, iy
!-------------------------------------------------------------------------------
!
   nelem = cs%nelem
   pvar = -1.0_r8
   do ig = 1,nelem
     call ig_to_iface_ie_je(cs,ig,iface,ie,je)
     call iface_ie_je_to_ix_iy(cs,iface,ie,je,ix,iy)
     pvar(:,:,ix,iy) = var(:,:,ig)
   enddo
!
   return
   end subroutine convert_global_elem_to_plot_2d_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine gen_lonslats(cs, lons, lats, nbr_lons, nbr_lats, cell_lons, cell_lats)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t),                  intent(in   ) :: cs
   real(r8), dimension(:),   allocatable, intent(inout) :: lons, lats
   real(r8), dimension(:,:), allocatable, intent(inout) :: nbr_lons, nbr_lats
   real(r8), dimension(:,:), allocatable, intent(inout) :: cell_lons, cell_lats
! local variables
   integer(i4) :: nups
   integer(i4) :: nface, ne, np, iface, ie, je, ip, jp, ii, jj
!
   nups = cs%nups
   nface = cs%nface
   ne = cs%ne
   np = cs%np
   allocate(lons(nups),lats(nups))
   allocate(nbr_lons(nnbrs,nups),nbr_lats(nnbrs,nups))
   allocate(cell_lons(ncell,nups),cell_lats(ncell,nups))
!
   ii = 0
   do iface = 1,nface
   do je = 1,ne
   do ie = 1,ne
   do jp = 1,np
   do ip = 1,np
     if (cs%elem(ie,je,iface)%qp(ip,jp)%isunique) then
       ii = ii+1
       if (cs%elem(ie,je,iface)%qp(ip,jp)%iup.ne.ii) then
         print*,iface,ie,je,ip,jp
         print*,'not match up...',cs%elem(ie,je,iface)%qp(ip,jp)%iup,ii
         stop
       endif
       lons(ii) = cs%elem(ie,je,iface)%qp(ip,jp)%lon
       lats(ii) = cs%elem(ie,je,iface)%qp(ip,jp)%lat
       do jj = 1,nnbrs
         nbr_lons(jj,ii) = cs%elem(ie,je,iface)%qp(ip,jp)%nbrs(jj)%lon
         nbr_lats(jj,ii) = cs%elem(ie,je,iface)%qp(ip,jp)%nbrs(jj)%lat
       enddo
       do jj = 1,ncell
         cell_lons(jj,ii) = cs%elem(ie,je,iface)%qp(ip,jp)%cell(jj)%lon
         cell_lats(jj,ii) = cs%elem(ie,je,iface)%qp(ip,jp)%cell(jj)%lat
       enddo
     endif
   enddo
   enddo
   enddo
   enddo
   enddo
!
   return
   end subroutine gen_lonslats
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_write_lonslats(cs, filename)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(in   ) :: cs
   character(len=*),     intent(in   ) :: filename
! local variables
   integer(i4) :: ierr
   integer(i4) :: fileid
   integer(i4) :: did_np, did_ne, did_nups, did_nnbrs, did_ncell
   integer(i4) :: vid_lons, vid_lats, vid_nbr_lons, vid_nbr_lats, vid_cell_lons, vid_cell_lats
   integer(i4) :: np, ne, nups!, nnbrs
   real(r8), dimension(:), allocatable :: lons, lats
   real(r8), dimension(:,:), allocatable :: nbr_lons, nbr_lats
   real(r8), dimension(:,:), allocatable :: cell_lons, cell_lats
!
   np = cs%np
   ne = cs%ne
   nups = cs%nups
!
   call gen_lonslats(cs,lons,lats,nbr_lons,nbr_lats,cell_lons,cell_lats)
!
   ! file create
   ierr = nf90_create(trim(filename),ior(nf90_clobber,nf90_64bit_offset),fileid)
   ! demensions
   ierr = nf90_def_dim(fileid,'np',np,did_np)
   ierr = nf90_def_dim(fileid,'ne',ne,did_ne)
   ierr = nf90_def_dim(fileid,'nups',nups,did_nups)
   ierr = nf90_def_dim(fileid,'nnbrs',nnbrs,did_nnbrs)
   ierr = nf90_def_dim(fileid,'ncell',ncell,did_ncell)
   ! define variable
   ierr = nf90_def_var(fileid,'lons',nf90_real8,did_nups,vid_lons)
   ierr = nf90_def_var(fileid,'lats',nf90_real8,did_nups,vid_lats)
   ierr = nf90_def_var(fileid,'nbr_lons',nf90_real8,(/did_nnbrs,did_nups/),vid_nbr_lons)
   ierr = nf90_def_var(fileid,'nbr_lats',nf90_real8,(/did_nnbrs,did_nups/),vid_nbr_lats)
   ierr = nf90_def_var(fileid,'cell_lons',nf90_real8,(/did_ncell,did_nups/),vid_cell_lons)
   ierr = nf90_def_var(fileid,'cell_lats',nf90_real8,(/did_ncell,did_nups/),vid_cell_lats)
   !
   ierr = nf90_enddef(fileid)
   ! put variables
   ierr = nf90_put_var(fileid,vid_lons,lons)
   ierr = nf90_put_var(fileid,vid_lats,lats)
   ierr = nf90_put_var(fileid,vid_nbr_lons,nbr_lons)
   ierr = nf90_put_var(fileid,vid_nbr_lats,nbr_lats)
   ierr = nf90_put_var(fileid,vid_cell_lons,cell_lons)
   ierr = nf90_put_var(fileid,vid_cell_lats,cell_lats)
   ! file close
   ierr = nf90_close(fileid)
!
   deallocate(lons,lats)
   deallocate(nbr_lons,nbr_lats)
!
   return
   end subroutine cs_write_lonslats
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_write_scrip(cs, filename)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t),       intent(in   ) :: cs
   character(len=*), optional, intent(in   ) :: filename
! local variables
   integer(i4)        :: nups, iup, ncorners, ico
   type(scrip_t)      :: scrip
   character(len=16)  :: str_rotated
   character(len=512) :: oufilename, description
   real(r8), dimension(:),   allocatable :: lons, lats
   real(r8), dimension(:,:), allocatable :: nbr_lons, nbr_lats
   real(r8), dimension(:,:), allocatable :: cell_lons, cell_lats
!
   if (cs%isrotated) then
     str_rotated = 'rotated'
   else
     str_rotated = 'unrotated'
   endif
!
   if (present(filename)) then
     write(oufilename,'(a)') trim(filename)
   else
     write(oufilename,'(a,i3.3,a)') './jhkim_ne',cs%ne,'_'//trim(str_rotated)//'.nc'
   endif
!
   write(description,'(a,i3.3,a)') 'Cubed-sphere on ne',cs%ne,' and np4('//trim(str_rotated)//') '
!
   call scrip_initialize(scrip,trim(oufilename),1,cs%nups,cs%ncell,trim(description))
!
   call gen_lonslats(cs,lons,lats,nbr_lons,nbr_lats,cell_lons,cell_lats)
!
   nups = cs%nups
   ! dimension
   scrip%grid_dims(:) = 1
   ! centers
   do iup = 1,nups
     scrip%grid_imask(iup) = 1
     !radian
     scrip%grid_center_lon(iup) = lons(iup)
     scrip%grid_center_lat(iup) = lats(iup)
   enddo
   ! corners
   ncorners = ncell
   do iup = 1,nups
   do ico = 1,ncorners
     scrip%grid_corner_lon(ico,iup) = cell_lons(ico,iup)
     scrip%grid_corner_lat(ico,iup) = cell_lats(ico,iup)
   enddo
   enddo
!
   call scrip_check(scrip)
!
   call scrip_write(scrip)
!
   call scrip_finalize(scrip)
!
   return
   end subroutine cs_write_scrip
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_get_distance(cs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t),                    intent(in   ) :: cs
! local variables
   integer(i4) :: np, ne, nface, npts
   integer(i4) :: ie, je, ip, jp, ii, jj
   real(r8)    :: rr, lon1, lat1, lon2, lat2, ang, dmin, dmax, dave
   real(r8), dimension(:,:,:), allocatable :: dist
   real(r8), dimension(:,:)  , allocatable :: dist_equ
!
   nface = cs%nface
   ne = cs%ne
   np = cs%np
   allocate(dist(2*np*(np-1),ne,ne))
   allocate(dist_equ(np-1,ne))
   dist     = 0.0_r8
   dist_equ = 0.0_r8
!
   ! cos(ang) = x dot y / |x||y|
   ! x dot y = r1r2(cosY1*cosY2*cos(X1-X2)+sinY1*sinY2)
   ! cos(ang) = cosY1*cosY2*cos(X2-X1)+sinY1*sinY2
   rr   = 6371220.0_r8
   do je = 1,ne
     do ie = 1,ne
       ii = 0
       jj = 0
       do jp = 1,np
       do ip = 1,np-1
         lon1 = cs%elem(ie,je,1)%qp(ip,jp)%lon
         lat1 = cs%elem(ie,je,1)%qp(ip,jp)%lat
         lon2 = cs%elem(ie,je,1)%qp(ip+1,jp)%lon
         lat2 = cs%elem(ie,je,1)%qp(ip+1,jp)%lat
         ii = ii+1
         ang = acos(cos(lat1)*cos(lat2)*cos(lon1-lon2)+sin(lat1)*sin(lat2))
         dist(ii,ie,je) = rr*ang
         if (je.eq.ne/2.and.jp.eq.np/2) then
           jj = jj+1
           dist_equ(jj,ie) = dist(ii,ie,je)
         endif
       enddo
       enddo
       do jp = 1,np-1
       do ip = 1,np
         lon1 = cs%elem(ie,je,1)%qp(ip,jp)%lon
         lat1 = cs%elem(ie,je,1)%qp(ip,jp)%lat
         lon2 = cs%elem(ie,je,1)%qp(ip,jp+1)%lon
         lat2 = cs%elem(ie,je,1)%qp(ip,jp+1)%lat
         ii = ii+1
         ang = acos(cos(lat1)*cos(lat2)*cos(lon1-lon2)+sin(lat1)*sin(lat2))
         dist(ii,ie,je) = rr*ang
       enddo
       enddo
#if 0
if (ie.eq.1.and.je.eq.1) then
  lon1 = cs%elem(ie,je,1)%qp(1,1)%lon
  lat1 = cs%elem(ie,je,1)%qp(1,1)%lat
  lon2 = cs%elem(ie,je,1)%qp(2,2)%lon
  lat2 = cs%elem(ie,je,1)%qp(2,2)%lat
  ang = acos(cos(lat1)*cos(lat2)*cos(lon1-lon2)+sin(lat1)*sin(lat2))
  print *,'DEBUG = ',rr*ang/1000.
endif
#endif
     enddo
   enddo
!
   dmin = minval(dist)
   dmax = maxval(dist)
   dave = sum(dist)/(2*np*(np-1)*ne*ne)
   print '(a,f7.2,x,f7.2,x,f7.2,x,f7.2)', 'min/max/max-min/average (entire) = ',dmin/1000.,dmax/1000.,(dmax-dmin)/1000.,dave/1000.
   dmin = minval(dist_equ)
   dmax = maxval(dist_equ)
   dave = sum(dist_equ)/((np-1)*ne)
   print '(a,f7.2,x,f7.2,x,f7.2,x,f7.2)', 'min/max/max-min/average (equator)= ',dmin/1000.,dmax/1000.,(dmax-dmin)/1000.,dave/1000.
!
   deallocate(dist,dist_equ)
!
   return
   end subroutine cs_get_distance
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_do_check(cs, np, ne, nprocs, filename1, filename2)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(inout) :: cs
   integer(i4),          intent(in   ) :: np, ne, nprocs
   character(len=*),     intent(in   ) :: filename1, filename2
! local variables
   real(r8), dimension(:)  , allocatable :: lons1, lats1, lons2, lats2
   real(r8), dimension(:,:), allocatable :: cell_lons1, cell_lats1, cell_lons2, cell_lats2
   integer(i4) :: fid, vid, ierr
   integer(i4) :: iface, ie, je, ip, jp, iup, ic
   real(r8) :: err
!
   call cs_initialize(cs,np,ne,nprocs,.false.)
   call set_coordinates(cs) ! iface,global index(global),cartesian index(ie,je),alpha,beta,lon,lat
   call set_sfc_index(cs)   ! isfc,joiner
   call set_rank(cs)        ! rank
   call set_local_index(cs) !
   ! local variables
   ! index(local)
   call set_neighbors(cs)
   call set_quadrature_ups_eps(cs)
!
   allocate(lons1(cs%nups))
   allocate(lats1(cs%nups))
   allocate(cell_lats1(cs%ncell,cs%nups))
   allocate(cell_lons1(cs%ncell,cs%nups))
   allocate(lons2(cs%nups))
   allocate(lats2(cs%nups))
   allocate(cell_lats2(cs%ncell,cs%nups))
   allocate(cell_lons2(cs%ncell,cs%nups))
!
   ! file 1
   ierr = nf90_open(trim(filename1),nf90_nowrite,fid)
   ierr = nf90_inq_varid(fid,'grid_center_lon',vid)
   ierr = nf90_get_var(fid,vid,lons1)
   ierr = nf90_inq_varid(fid,'grid_center_lat',vid)
   ierr = nf90_get_var(fid,vid,lats1)
   ierr = nf90_inq_varid(fid,'grid_corner_lon',vid)
   ierr = nf90_get_var(fid,vid,cell_lons1)
   ierr = nf90_inq_varid(fid,'grid_corner_lat',vid)
   ierr = nf90_get_var(fid,vid,cell_lats1)
   ierr = nf90_close(fid)
!
   ! file 2
   ierr = nf90_open(trim(filename2),nf90_nowrite,fid)
   ierr = nf90_inq_varid(fid,'grid_center_lon',vid)
   ierr = nf90_get_var(fid,vid,lons2)
   ierr = nf90_inq_varid(fid,'grid_center_lat',vid)
   ierr = nf90_get_var(fid,vid,lats2)
   ierr = nf90_inq_varid(fid,'grid_corner_lon',vid)
   ierr = nf90_get_var(fid,vid,cell_lons2)
   ierr = nf90_inq_varid(fid,'grid_corner_lat',vid)
   ierr = nf90_get_var(fid,vid,cell_lats2)
   ierr = nf90_close(fid)
!
   iup = 1
   do iface = 1,cs%nface
   do je = 1,cs%ne
   do ie = 1,cs%ne
   do jp = 1,cs%np
   do ip = 1,cs%np
     if (cs%elem(ie,je,iface)%qp(ip,jp)%isunique) then
       err =(lons1(iup)-lons2(iup))/lons1(iup)*100.0_r8
       if (lons1(iup).ne.lons2(iup)) then !.and.abs(err).gt.1.0d-15) then
         print*,'diff lons : ',iface,ie,je,ip,jp,iup,lons1(iup),lons2(iup),(lons1(iup)-lons2(iup))/lons1(iup)*100.0_r8
       endif
  !  
       err =(lats1(iup)-lats2(iup))/lats1(iup)*100.0_r8
       if (lats1(iup).ne.lats2(iup)) then !.and.abs(err).gt.1.0d-15) then
         print*,'diff lats : ',iface,ie,je,ip,jp,iup,lats1(iup),lats2(iup),(lats1(iup)-lats2(iup))/lats1(iup)*100.0_r8
       endif
  !
       do ic = 1,cs%ncell
         err =(cell_lons1(ic,iup)-cell_lons2(ic,iup))/cell_lons1(ic,iup)*100.0_r8
         if (cell_lons1(ic,iup).ne.cell_lons2(ic,iup)) then !.and.abs(err).gt.1.0d-15) then
           print*,'diff cell_lons : ',iface,ie,je,ip,jp,iup,cell_lons1(ic,iup),cell_lons2(ic,iup),(cell_lons1(ic,iup)-cell_lons2(ic,iup))/cell_lons1(ic,iup)*100.0_r8
         endif
         err =(cell_lats1(ic,iup)-cell_lats2(ic,iup))/cell_lats1(ic,iup)*100.0_r8
         if (cell_lats1(ic,iup).ne.cell_lats2(ic,iup)) then !.and.abs(err).gt.1.0d-15) then
           print*,'diff cell_lats : ',iface,ie,je,ip,jp,iup,cell_lats1(ic,iup),cell_lats2(ic,iup),(cell_lats1(ic,iup)-cell_lats2(ic,iup))/cell_lats1(ic,iup)*100.0_r8
         endif
       enddo
   ! 
       iup = iup+1
     endif
   enddo
   enddo
   enddo
   enddo
   enddo
!
   call cs_finalize(cs)
!
   deallocate(lons1)
   deallocate(lats1)
   deallocate(cell_lats1)
   deallocate(cell_lons1)
   deallocate(lons2)
   deallocate(lats2)
   deallocate(cell_lats2)
   deallocate(cell_lons2)
!
   return
   end subroutine cs_do_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_find_index_point(cs, lond, latd)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cubed_sphere_t), intent(in   ) :: cs
   real(r8),             intent(in   ) :: lond, latd
! local variables
   integer(i4) :: np, ne, nface, nups
   integer(i4) :: ie, je, iface, ip, jp, iup, pmin(1)
   real(r8)    :: rr, lonf, latf, lon1, lat1, ang, rmin, dave
   real(r8), dimension(:), allocatable :: dist
!
   nups  = cs%nups
   nface = cs%nface
   ne    = cs%ne
   np    = cs%np
   allocate(dist(nups))
   dist  = 0.0_r8
!
   lonf = lond/180.0*pi
   latf = latd/180.0*pi
   ! cos(ang) = x dot y / |x||y|
   ! x dot y = r1r2(cosY1*cosY2*cos(X1-X2)+sinY1*sinY2)
   ! cos(ang) = cosY1*cosY2*cos(X2-X1)+sinY1*sinY2
   rr   = 6371220.0_r8
   iup  = 0
   do iface =1,nface
    do je = 1,ne
     do ie = 1,ne
      do jp = 1,np
       do ip = 1,np
         if (cs%elem(ie,je,iface)%qp(ip,jp)%isunique) then
           iup = iup+1
           lon1 = cs%elem(ie,je,iface)%qp(ip,jp)%lon
           lat1 = cs%elem(ie,je,iface)%qp(ip,jp)%lat
           ang = acos(cos(lat1)*cos(latf)*cos(lon1-lonf)+sin(lat1)*sin(latf))
           dist(iup) = rr*ang
         endif
       enddo
      enddo
     enddo
    enddo
   enddo
!
   rmin = minval(dist)
   pmin = minloc(dist)
   print '(a,f6.1,a,f6.1,a,i8,a,f5.1,a)', &
         ' - lat/lon = ',latd,'/ ',lond,': iup = ',pmin(1),' (dist = ',rmin/1000.0,' km)'
!
   deallocate(dist)
!
   return
   end subroutine cs_find_index_point
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module cubed_sphere
!-------------------------------------------------------------------------------
