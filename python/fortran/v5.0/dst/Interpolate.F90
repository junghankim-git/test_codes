#include "KIM.h"
!-------------------------------------------------------------------------------
   module interpolate
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
   use kiapsbase, only: int_kind=>kim_int_kind, iulog=>kim_iu_log, real_kind=>kim_real8_kind
   use element, only: element_t
   use dimensions, only: np, ne, nelemd
   use quadrature, only: quadrature_t, legendre, quad_norm
   use coordinatesystems, only: spherical_polar_t, cartesian2d_t, cartesian3d_t, sphere2cubedsphere, spherical_to_cart, cubedsphere2cart, distance, ref2sphere_double
   use physicalconstants, only: dd_pi
   use quadrature, only: quadrature_t, gauss, gausslobatto
   use kiapsparallel, only: abortpar, barrierpar, parallel_t, par_double_precision, par_max, par_integer, par_sum, par_min
   use cube, only: convert_gbl_index, dmap
   use mesh, only: meshusemeshfile
   implicit none
!
   private
   integer, parameter, public :: max_vecvars = 25
   real(real_kind), parameter, private :: tol = 1e-10
   character(len=10), public :: vector_uvars(max_vecvars), vector_vvars(max_vecvars)
!
   type, public :: interpolate_t
     sequence
     real(real_kind), dimension(:,:), pointer :: imat ! p_k(xj)*wj/gamma(k)
     real(real_kind), dimension(:), pointer :: rk ! 1/k
     real(real_kind), dimension(:), pointer :: vtemp ! temp results
     real(real_kind), dimension(:), pointer :: glp ! gll pts(nair)
   end type interpolate_t
!
   type, public :: interpdata_t
! Output Interpolation points.  Used to output data on lat-lon (or other grid)
! with native element interpolation.  Each element keeps a list of points from the
! interpolation grid that are in this element
     type(cartesian2d_t), pointer, dimension(:) :: interp_xy ! element coordinate
     integer, pointer, dimension(:) :: ilat, ilon ! position of interpolation point in lat-lon grid
     integer :: n_interp
     integer :: nlat
     integer :: nlon
     logical :: first_entry = .true.
   end type interpdata_t
   real(real_kind), private :: delta = 1.0d-9 ! move tiny bit off center to
! avoid landing on element edges
   public :: interp_init
   public :: setup_latlon_interp
   public :: interpolate_scalar
   public :: interpolate_vector
   public :: set_interp_parameter
   public :: get_interp_parameter
   public :: get_interp_gweight
   public :: get_interp_lat
   public :: get_interp_lon
   public :: var_is_vector_uvar, var_is_vector_vvar
   public :: cube_facepoint_ne
   private :: parametric_coordinates
   private :: find_parametric_by_roots
   public :: interpolate_2d
   public :: interpolate_create
!
   interface interpolate_scalar
     module procedure interpolate_scalar2d
     module procedure interpolate_scalar3d
   end interface
!
   interface interpolate_vector
     module procedure interpolate_vector2d
     module procedure interpolate_vector3d
   end interface
   type(interpolate_t), target :: interp_p
! store the  lat-lon grid
! gridtype = 1       equally spaced, including poles (FV scalars output grid)
! gridtype = 2       Gauss grid (CAM Eulerian)
! gridtype = 3       equally spaced, no poles (FV staggered velocity)
! Seven possible history files, last one is inithist and should be native grid
   logical, public, save :: interpolate_analysis(8) =(/.true., .false., .false., .false., .false., .false., .false., .false./)
   integer :: nlat, nlon
   real(real_kind), pointer, dimension(:) :: lat(:), lon(:), gweight(:)
   integer :: gridtype = 1 !
   integer :: itype = 1 ! 0 = native high order
! 1 = bilinear
!
   contains
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine set_interp_parameter(parm_name, value)
!
   character*(*), intent(in   ) :: parm_name
   character(len=80) :: msg
   integer :: value, power
   real(real_kind) :: value_target
!
   if (parm_name.eq.'itype') then
     itype = value
   elseif (parm_name.eq.'nlon') then
     nlon = value
   elseif (parm_name.eq.'nlat') then
     nlat = value
   elseif (parm_name.eq.'gridtype') then
     gridtype = value
   elseif (parm_name.eq.'auto') then
! compute recommended nlat,nlon which has slightly higher
! resolution than the specifed number of points around equator given in "value"
! computed recommended lat-lon grid.
! nlon > peq   peq = points around equator cubed sphere grid
! take nlon power of 2, and at most 1 power of 3
     if (value.eq.0) then
! If reading in unstructured mesh, ne = 0
! This makes it hard to guess how many interpolation points to use
! So We'll set the default as 720 x 360
! BUT if you're running with an unstructured mesh, set interp_nlon and interp_nlat
       nlon = 720
       nlat = 360
     else
       value_target = value*1.25
       power = nint(.5+log(value_target)/log(2d0))
       power = max(power,7) ! min grid:64x128
       if (3*2**(power-2)>value_target) then
         nlon = 3*2**(power-2) ! use 1 power of 3
       else
         nlon = 2**power
       endif
     endif
       nlat = nlon/2
       if (gridtype==1) nlat = nlat+1
     else
       write(msg,*) 'Did not recognize parameter named ',parm_name,' in interpolate_mod:set_interp_parameter'
       call abortpar(message = msg)
   endif
!
   end subroutine set_interp_parameter
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function get_interp_parameter(parm_name) result(value)
!
   character*(*), intent(in   ) :: parm_name
   integer :: value
   character(len=80) :: msg
!
   if (parm_name.eq.'itype') then
     value = itype
   elseif (parm_name.eq.'nlon') then
     value = nlon
   elseif (parm_name.eq.'nlat') then
     value = nlat
   elseif (parm_name.eq.'gridtype') then
     value = gridtype
   else
     write(msg,*) 'Did not recognize parameter named ',parm_name,' in interpolate_mod:get_interp_parameter'
     value =-1
     call abortpar(message = msg)
   endif
!
   return
!
   end function get_interp_parameter
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function get_interp_gweight() result(gw)
!
   real(real_kind) :: gw(nlat)
!
   gw = gweight
   return
!
   end function get_interp_gweight
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function get_interp_lat() result(thislat)
!-------------------------------------------------------------------------------
   use physicalconstants, only: dd_pi
!
   real(real_kind) :: thislat(nlat)
!
   thislat = lat*180.0d0/dd_pi
   return
!
   end function get_interp_lat
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function get_interp_lon() result(thislon)
!-------------------------------------------------------------------------------
   use physicalconstants, only: dd_pi
!
   real(real_kind) :: thislon(nlon)
!
   thislon = lon*180.0d0/dd_pi
   return
!
   end function get_interp_lon
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine interpolate_create(gquad, interp)
!
   type(quadrature_t), intent(in   ) :: gquad
   type(interpolate_t), intent(  out) :: interp
! Local variables
   integer k, j
   integer npts
   real(real_kind), dimension(:), allocatable :: gamma
   real(real_kind), dimension(:), allocatable :: leg
!
   npts = size(gquad%points)
   allocate(interp%imat(npts,npts))
   allocate(interp%rk(npts))
   allocate(interp%vtemp(npts))
   allocate(interp%glp(npts))
   allocate(gamma(npts))
   allocate(leg(npts))
   gamma = quad_norm(gquad,npts)
   do k = 1,npts
     interp%rk(k) = 1.0d0/k
     interp%glp(k) = gquad%points(k) !nair
   enddo
   do j = 1,npts
     leg = legendre(gquad%points(j),npts-1)
     do k = 1,npts
       interp%imat(j,k) = leg(k)*gquad%weights(j)/gamma(k)
     enddo
   enddo
   deallocate(gamma)
   deallocate(leg)
!
   end subroutine interpolate_create
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function interpolate_2d(cart, f, interp, npts, fillvalue) result(fxy)
!
   integer, intent(in   ) :: npts
   type(cartesian2d_t), intent(in   ) :: cart
   real(real_kind), intent(in   ) :: f(npts, npts)
   type(interpolate_t) :: interp
   real(real_kind) :: fxy ! value of f interpolated to(x, y)
   real(real_kind), intent(in   ), optional :: fillvalue
! local variables
   real(real_kind) :: tmp_1, tmp_2
   real(real_kind) :: fk0, fk1
   real(real_kind) :: pk
   integer :: l, j, k
!
   if (present(fillvalue)) then
     if (any(f==fillvalue)) then
       fxy = fillvalue
       return
     endif
   endif
!
   do l = 1, npts, 2
! Compute Pk(cart%x) for Legendre order 0
     pk = 1.0d0
     fk0 = 0.0d0
     fk1 = 0.0d0
     do j = 1,npts
       fk0 = fk0+interp%imat(j,1)*f(j,l)
       fk1 = fk1+interp%imat(j,1)*f(j,l+1)
     enddo
     interp%vtemp(l) = pk*fk0
     interp%vtemp(l+1) = pk*fk1
! Compute Pk(cart%x) for Legendre order 1
     tmp_2 = pk
     pk = cart%x
     fk0 = 0.0d0
     fk1 = 0.0d0
     do j = 1,npts
       fk0 = fk0+interp%imat(j,2)*f(j,l)
       fk1 = fk1+interp%imat(j,2)*f(j,l+1)
     enddo
     interp%vtemp(l) = interp%vtemp(l)+pk*fk0
     interp%vtemp(l+1) = interp%vtemp(l+1)+pk*fk1
! Compute Pk(cart%x) for Legendre order 2 to npts-1
     do k = 2,npts-1
       tmp_1 = tmp_2
       tmp_2 = pk
       pk =((2*k-1)*cart%x*tmp_2-(k-1)*tmp_1)*interp%rk(k)
       fk0 = 0.0d0
       fk1 = 0.0d0
       do j = 1,npts
         fk0 = fk0+interp%imat(j,k+1)*f(j,l)
         fk1 = fk1+interp%imat(j,k+1)*f(j,l+1)
       enddo
       interp%vtemp(l) = interp%vtemp(l)+pk*fk0
       interp%vtemp(l+1) = interp%vtemp(l+1)+pk*fk1
     enddo
   enddo
! Compute Pk(cart%y) for Legendre order 0
!
   pk = 1.0
   fk0 = 0.0d0
   do j = 1,npts
     fk0 = fk0+interp%imat(j,1)*interp%vtemp(j)
   enddo
   fxy = pk*fk0
! Compute Pk(cart%y) for Legendre order 1
   tmp_2 = pk
   pk = cart%y
   fk0 = 0.0d0
   do j = 1,npts
     fk0 = fk0+interp%imat(j,2)*interp%vtemp(j)
   enddo
   fxy = fxy+pk*fk0
! Compute Pk(cart%y) for Legendre order 2,npts-1
   do k = 2,npts-1
     tmp_1 = tmp_2
     tmp_2 = pk
     pk =((2*k-1)*cart%y*tmp_2-(k-1)*tmp_1)*interp%rk(k)
     fk0 = 0.0d0
     do j = 1,npts
       fk0 = fk0+interp%imat(j,k+1)*interp%vtemp(j)
     enddo
     fxy = fxy+pk*fk0
   enddo
!
   end function interpolate_2d
!===============================
!(Nair) Bilinear interpolation for every GLL grid cell
!===============================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function interpol_bilinear(cart, f, interp, npts, fillvalue) result(fxy)
!
   integer, intent(in   ) :: npts
   type(cartesian2d_t), intent(in   ) :: cart
   real(real_kind), intent(in   ) :: f(npts, npts)
   type(interpolate_t) :: interp
   real(real_kind) :: fxy ! value of f interpolated to(x, y)
   real(real_kind), intent(in   ), optional :: fillvalue
! local variables
   real(real_kind) :: xoy(npts)
   real(real_kind) :: p, q, xp, yp, y4(4)
   integer :: l, j, k, ii, jj, na, nb, nm
!
   xp = cart%x
   yp = cart%y
   xoy(:) = interp%glp(:)
! Search index along "x"  (bisection method)
   na = 1
   nb = npts
   do
     if ((nb-na)<=1) exit
     nm =(nb+na)/2
     if (xp>xoy(nm)) then
       na = nm
     else
       nb = nm
     endif
   enddo
   ii = na
! Search index along "y"
   na = 1
   nb = npts
   do
     if ((nb-na)<=1) exit
     nm =(nb+na)/2
     if (yp>xoy(nm)) then
       na = nm
     else
       nb = nm
     endif
   enddo
   jj = na
! GLL cell containing (xp,yp)
   y4(1) = f(ii,jj)
   y4(2) = f(ii+1,jj)
   y4(3) = f(ii+1,jj+1)
   y4(4) = f(ii,jj+1)
   if (present(fillvalue)) then
     if (any(y4==fillvalue)) then
       fxy = fillvalue
       return
     endif
   endif
   p =(xp-xoy(ii))/(xoy(ii+1)-xoy(ii))
   q =(yp-xoy(jj))/(xoy(jj+1)-xoy(jj))
   fxy =(1.0d0-p)*(1.0d0-q)*y4(1)+p*(1.0d0-q)*y4(2)+p*q*y4(3)+(1.0d0-p)*q*y4(4)
!
   end function interpol_bilinear
!
! fast iterative search for bilinear elements on gnomonic cube face
! (replaced by more general version parametric_coordinates2)
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function parametric_coordinates(sphere, corners, face_no) result(par_coor)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian2d_t), intent(in   ) :: corners(4)
   type(spherical_polar_t), intent(in   ) :: sphere
   integer :: face_no
   type(cartesian2d_t) :: par_coor
   real(real_kind) :: isinelemconverged
   real(real_kind) :: x(4), y(4)
   real(real_kind) :: xp, yp
   real(real_kind) :: j(4), f(2), shapefct(4)
   real(real_kind) :: xinew, etanew
   real(real_kind) :: xicur, etacur
   real(real_kind) :: jdet
   real(real_kind) :: xidiff(2)
   type(cartesian2d_t) :: cube
   integer :: i
   integer :: max_nr_iter
!
   isinelemconverged = 1.0e-16
   max_nr_iter = 10 ! bilinear so should converge in 1 step
   cube = sphere2cubedsphere(sphere,face_no)
! Translate element so that (x,y) coordinates of the first node are (0,0)
   x(1) = 0
   x(2) = corners(2)%x-corners(1)%x
   x(3) = corners(3)%x-corners(1)%x
   x(4) = corners(4)%x-corners(1)%x
   y(1) = 0
   y(2) = corners(2)%y-corners(1)%y
   y(3) = corners(3)%y-corners(1)%y
   y(4) = corners(4)%y-corners(1)%y
! (xp,yp) is the point at which we're searching for (xi,eta)
! (must translate this also)
   xp = cube%x-corners(1)%x
   yp = cube%y-corners(1)%y
! Newton-Raphson iteration for (xi,eta)
   xinew = 0.5 ! initial guess
   etanew = 0.5
   xicur = 0.5
   etacur = 0.5
   xidiff(1) = 1.0
   xidiff(2) = 1.0
   i = 1
   do
     xicur = xinew
     etacur = etanew
     j(1) = 0.25*(1.00-etacur)*x(2)+0.25*(1.00+etacur)*x(3)-0.25*(1.00+etacur)*x(4)
     j(2) =-0.25*(1.00+xicur)*x(2)+0.25*(1.00+xicur)*x(3)+0.25*(1.00-xicur)*x(4)
     j(3) = 0.25*(1.00-etacur)*y(2)+0.25*(1.00+etacur)*y(3)-0.25*(1.00+etacur)*y(4)
     j(4) =-0.25*(1.00+xicur)*y(2)+0.25*(1.00+xicur)*y(3)+0.25*(1.00-xicur)*y(4)
     jdet = j(1)*j(4)-j(2)*j(3)
     shapefct(1) = 0.25*(1.00-etacur)*(1.00-xicur)
     shapefct(2) = 0.25*(1.00-etacur)*(1.00+xicur)
     shapefct(3) = 0.25*(1.00+etacur)*(1.00+xicur)
     shapefct(4) = 0.25*(1.00+etacur)*(1.00-xicur)
     f(1) =(shapefct(2)*x(2)+shapefct(3)*x(3)+shapefct(4)*x(4))-xp
     f(2) =(shapefct(2)*y(2)+shapefct(3)*y(3)+shapefct(4)*y(4))-yp
     xinew = xicur-(f(1)*j(4)-f(2)*j(2))/jdet
     etanew = etacur-(-f(1)*j(3)+f(2)*j(1))/jdet
     xidiff(1) = xinew-xicur
     xidiff(2) = etanew-etacur
     i = i+1
     if (dot_product(xidiff,xidiff)<isinelemconverged.or.max_nr_iter<i) exit
   enddo
   par_coor%x = huge(par_coor%x)
   par_coor%y = huge(par_coor%y)
   if (i<=max_nr_iter) then
     par_coor%x = xinew
     par_coor%y = etanew
   endif
!
   end function parametric_coordinates
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function parametric_coordinates2(sphere, elem) result(ref)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spherical_polar_t), intent(in   ) :: sphere
   type(element_t), intent(in   ) :: elem
   type(cartesian2d_t) :: ref
! local
   integer :: face_no, i, max_nr_iter = 10
   real(real_kind) :: d(2, 2), dinv(2, 2), detd, a, b, resa, resb, dela, delb, costh
   real(real_kind) :: tol_sq = 1e-26
   type(spherical_polar_t) :: sphere1, sphere_tmp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! newton iteration on: ref=ref - df^-1 (ref2sphere(ref) - sphere)
!
! Generic version written in terms of HOMME's 'ref2sphere' and 'Dmap' operaters,
! with no assumption as to the type of map (gnomonic, equi-angular, parametric)
!
! f = ref2sphere(xvec) - sphere
! df = d(ref2sphere)
!
! D = diag(cos(theta),1) * d(ref2sphere)       d(ref2sphere) = diag(1/cos(theta),1)*D
!
! df = diag(1/cos(theta),1)*D
! df^-1 =  D^-1 *  diag(cos(theta),1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   costh = cos(sphere%lat)
   a = 0
   b = 0
   i = 0
   do
     sphere1 = ref2sphere_double(a,b,elem%corners,elem%vertex%face_number)
     resa = sphere1%lon-sphere%lon
     if (resa>dd_pi) resa = resa-2*dd_pi
     if (resa<-dd_pi) resa = resa+2*dd_pi
     resb = sphere1%lat-sphere%lat
     call dmap(d,elem,a,b)
     detd = d(1,1)*d(2,2)-d(1,2)*d(2,1)
     dinv(1,1) = d(2,2)/detd
     dinv(1,2) =-d(1,2)/detd
     dinv(2,1) =-d(2,1)/detd
     dinv(2,2) = d(1,1)/detd
     dela = dinv(1,1)*costh*resa+dinv(1,2)*resb
     delb = dinv(2,1)*costh*resa+dinv(2,2)*resb
     a = a-dela
     b = b-delb
     i = i+1
     if ((costh*resa)**2+resb**2<tol_sq.or.max_nr_iter<i) exit
   enddo
   ref%x = a
   ref%y = b
!
   end function parametric_coordinates2
!
! analytic inverse of bilinear map from reference element to arbrirary quad
! on a cube face
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function find_parametric_by_roots(cart, elem, cube) result(found)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian2d_t), intent(in   ) :: cube
   type(element_t), intent(in   ) :: elem
   type(cartesian2d_t), intent(  out) :: cart
   logical :: found
   real(real_kind) :: xp, yp
! Store the two roots of inverse map
   real(real_kind) :: root1(2), root2(2)
   real(real_kind) :: aa, bb, cc, dd, ee, ff, gg, hh, tmp1, tmp2, det, det2
!
   found = .false.
   xp = cube%x
   yp = cube%y
! For each element,find the two values in reference element that map to
! the interpolation point. If either of these points is in [-1,1]^2,we are
! in the right element. Otherwise,go on to the next one.
   aa = elem%u2qmap(1,1)-xp
   bb = elem%u2qmap(2,1)
   cc = elem%u2qmap(3,1)
   dd = elem%u2qmap(4,1)
   ee = elem%u2qmap(1,2)-yp
   ff = elem%u2qmap(2,2)
   gg = elem%u2qmap(3,2)
   hh = elem%u2qmap(4,2)
! Need to divide by (dd*ff-bb*hh) and (dd*gg-cc*hh),so check for when either
! value is zero:
   if ((abs(dd).lt.tol).and.(abs(hh).lt.tol)) then
! dd=hh=0 => dd*ff - bb*hh = 0 AND dd*gg - cc*hh = 0
     root1(1) =(cc*ee-aa*gg)/(bb*gg-cc*ff)
     root1(2) =(bb*ee-aa*ff)/(cc*ff-bb*gg)
     root2 = root1
   elseif ((abs(ff).lt.tol).and.(abs(hh).lt.tol)) then
! ff=hh=0 => dd*ff - bb*hh = 0
     root1(1) =(cc*ee-aa*gg)/(bb*gg-dd*ee)
     root1(2) =-ee/gg
     root2 = root1
   elseif ((abs(dd).lt.tol).and.(abs(bb).lt.tol)) then
! dd=bb=0 => dd*ff - bb*hh = 0
     root1(1) =(cc*ee-aa*gg)/(aa*hh-cc*ff)
     root1(2) =-aa/cc
     root2 = root1
   elseif ((abs(ff).lt.tol).and.(abs(bb).lt.tol)) then
! ff=bb=0 => dd*ff - bb*hh = 0
     root1(1) =(cc*ee-aa*gg)/(aa*hh-dd*ee)
     root1(2) =(dd*ee-aa*hh)/(cc*hh-dd*gg)
     root2 = root1
   elseif ((abs(gg).lt.tol).and.(abs(hh).lt.tol)) then
     root1(1) =-ee/ff
     root1(2) =(bb*ee-aa*ff)/(cc*ff-dd*ee)
     root2 = root1
   elseif ((abs(dd).lt.tol).and.(abs(cc).lt.tol)) then
     root1(1) =-aa/bb
     root1(2) =(bb*ee-aa*ff)/(aa*hh-bb*gg)
     root2 = root1
   elseif ((abs(gg).lt.tol).and.(abs(cc).lt.tol)) then
     root1(1) =(bb*ee-aa*ff)/(aa*hh-dd*ee)
     root1(2) =(aa*hh-dd*ee)/(bb*hh-dd*ff)
     root2 = root1
   elseif ((abs(dd*ff-bb*hh).lt.tol*tol).or.(abs(dd*gg-cc*hh).lt.tol*tol)) then
     print*,' this should not be zero!'
     print*,"bb = ",bb
     print*,"cc = ",cc
     print*,"dd = ",dd
     print*,"ff = ",ff
     print*,"gg = ",gg
     print*,"hh = ",hh
     print*,"dd*ff-bb*hh = ",dd*ff-bb*hh
     print*,"dd*gg-cc*hh = ",dd*gg-cc*hh
   else
     tmp1 = dd*ee+cc*ff-bb*gg-aa*hh
     tmp2 = dd*ee-cc*ff+bb*gg-aa*hh
     det2 = tmp2*tmp2-4.0d0*(bb*ee-aa*ff)*(dd*gg-cc*hh)
     if (det2.lt.0.0d0) then
! no real roots
       root1 =(/9999.0d0,9999.0d0/)
       root2 = root1
     else
       det = sqrt(det2)
       root1(1) =(-tmp1+det)/(2.0d0*(dd*ff-bb*hh))
       root1(2) =(tmp2+det)/(2.0d0*(-dd*gg+cc*hh))
       root2(1) =(-tmp1-det)/(2.0d0*(dd*ff-bb*hh))
       root2(2) =(tmp2-det)/(2.0d0*(-dd*gg+cc*hh))
     endif
   endif
   if ((root1(1)+1.0d0.ge.-tol).and.(root1(1)-1.0d0.le.tol).and.(root1(2)+1.0d0.ge.-tol).and.(root1(2)-1.0d0.le.tol)) then
     found = .true.
     cart%x = root1(1)
     cart%y = root1(2)
   elseif ((root2(1)+1.0d0.ge.-tol).and.(root2(1)-1.0d0.le.tol).and.(root2(2)+1.0d0.ge.-tol).and.(root2(2)-1.0d0.le.tol)) then
     found = .true.
     cart%x = root2(1)
     cart%y = root2(2)
   endif
!
   end function find_parametric_by_roots
!
! find element containing given point, useing HOMME's standard
! equi-angular gnomonic map.
! note that with this map, only coordinate lines are great circle arcs
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function point_inside_equiangular(elem, sphere, sphere_xyz) result(inside)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spherical_polar_t), intent(in   ) :: sphere
   type(cartesian3d_t), intent(in   ) :: sphere_xyz
   type(element_t), intent(in   ) :: elem
   logical :: inside, inside2
   logical :: debug = .false.
   integer :: i, j
   type(cartesian2d_t) :: corners(4), sphere_xy, cart
   type(cartesian3d_t) :: corners_xyz(4), center, a, b, cross(4)
   real(real_kind) :: yp(4), y, elem_diam, dotprod
   real(real_kind) :: xp(4), x
   real(real_kind) :: tol_inside
   real(real_kind) :: d1, d2
!
   inside = .false.
! first check if point is near the element:
   corners_xyz(:) = elem%corners3d(:)
   elem_diam = max(distance(corners_xyz(1),corners_xyz(3)),distance(corners_xyz(2),corners_xyz(4)))
   center%x = sum(corners_xyz(1:4)%x)/4
   center%y = sum(corners_xyz(1:4)%y)/4
   center%z = sum(corners_xyz(1:4)%z)/4
   if (distance(center,sphere_xyz)>1.0*elem_diam) return
   tol_inside = 1e-10*elem_diam**2
! the point is close to the element,so project both to cubed sphere
! and perform contour integral
   sphere_xy = sphere2cubedsphere(sphere,elem%facenum)
   x = sphere_xy%x
   y = sphere_xy%y
   do i = 1,4
     xp(i) = elem%corners(i)%x
     yp(i) = elem%corners(i)%y
   enddo
   if (debug) then
     j = 4
     do i = 1,4
       print*,i,(x-xp(j))*(yp(i)-yp(j))-(y-yp(j))*(xp(i)-xp(j))
       j = i ! within this loopk j = i-1
     enddo
   endif
   j = 4
   do i = 1,4
! a = x-xp(j), y-yp(j)
! b = xp(i)-xp(j), yp(i)-yp(j)
! compute a cross b:
     if ((x-xp(j))*(yp(i)-yp(j))-(y-yp(j))*(xp(i)-xp(j))>tol_inside) then
       return
     endif
       j = i ! within this loopk j = i-1
   enddo
! all cross products were negative,must be inside:
   inside = .true.
!
   end function point_inside_equiangular
!
! find element containing given point, with element edges assumed to be great circle arcs
! this will work with any map where straight lines are mapped to great circle arcs.
! (thus it will fail on unstructured grids using the equi-angular gnomonic map)
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function point_inside_gc(elem, sphere, sphere_xyz) result(inside)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spherical_polar_t), intent(in   ) :: sphere
   type(cartesian3d_t), intent(in   ) :: sphere_xyz
   type(element_t), intent(in   ) :: elem
   logical :: inside, inside2
   logical :: debug
   integer :: i, j
   type(cartesian2d_t) :: corners(4), sphere_xy, cart
   type(cartesian3d_t) :: corners_xyz(4), center, a, b, cross(4)
   real(real_kind) :: yp(4), y, elem_diam, dotprod
   real(real_kind) :: xp(4), x
   real(real_kind) :: d1, d2
!
   inside = .false.
! first check if point is near the element:
   corners_xyz(:) = elem%corners3d(:)
   elem_diam = max(distance(corners_xyz(1),corners_xyz(3)),distance(corners_xyz(2),corners_xyz(4)))
   center%x = sum(corners_xyz(1:4)%x)/4
   center%y = sum(corners_xyz(1:4)%y)/4
   center%z = sum(corners_xyz(1:4)%z)/4
   if (distance(center,sphere_xyz)>1.0*elem_diam) return
   j = 4
   do i = 1,4
! outward normal to plane containing j->i edge:  corner(i) x corner(j)
! sphere dot (corner(i) x corner(j) ) = negative if outside
     cross(i)%x = corners_xyz(i)%y*corners_xyz(j)%z-corners_xyz(i)%z*corners_xyz(j)%y
     cross(i)%y =-(corners_xyz(i)%x*corners_xyz(j)%z-corners_xyz(i)%z*corners_xyz(j)%x)
     cross(i)%z = corners_xyz(i)%x*corners_xyz(j)%y-corners_xyz(i)%y*corners_xyz(j)%x
     dotprod = cross(i)%x*sphere_xyz%x+cross(i)%y*sphere_xyz%y+cross(i)%z*sphere_xyz%z
     j = i ! within this loopk j = i-1
     if (dotprod<0) return
   enddo
   inside = .true.
   return
!
   end function point_inside_gc
!================================================
!  (Nair) Cube face index and local coordinates
!================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine cube_facepoint_ne(sphere, ne, cart, number)
!-------------------------------------------------------------------------------
   use coordinatesystems, only: cube_face_number_from_sphere, sphere2cubedsphere
   implicit none
!
   type(spherical_polar_t), intent(in   ) :: sphere
   integer, intent(in   ) :: ne
   type(cartesian2d_t), intent(  out) :: cart
   integer, intent(  out) :: number
   real(real_kind) :: xp, yp
   type(cartesian2d_t) :: cube
   integer :: ie, je, face_no
   real(real_kind) :: x1, x2
   real(real_kind) :: dx
!
   face_no = cube_face_number_from_sphere(sphere)
   cube = sphere2cubedsphere(sphere,face_no)
   xp = cube%x
   yp = cube%y
! MNL: for uniform grids (on cube face),analytic solution is fine
   x1 = xp+0.25d0*dd_pi
   x2 = yp+0.25d0*dd_pi
   dx =(0.5d0*dd_pi)/ne
   ie = int(abs(x1)/dx)
   je = int(abs(x2)/dx)
! if we are exactly on an element edge,we can put the point in
! either the ie or ie+1 element,EXCEPT if ie==ne.
   if (abs(x1)<ne*dx) ie = ie+1
   if (abs(x2)<ne*dx) je = je+1
   if (ie>ne.or.je>ne) then
     write(iulog,*) 'ERROR:',ie,je,ne
     write(iulog,*) 'lat,lon = ',sphere%lat,sphere%lon
     write(iulog,*) 'face no = ',face_no
     write(iulog,*) x1,x2,x1/dx,x2/dx
     call abortpar(message = 'interpolate_mod:bad argument')
   endif
! bug fix MT 1/2009.  This was creating a plotting error at
! the row of elements in iface=2 at 50 degrees (NE=16 128x256 lat/lon grid)
! For point on element edge,we can have ie=2,but x1=dx
! but if ie>1,we must execute this statement.
! The only time we can skip this statement is if ie=1,but then
! the statement has no effect,so lets never skip it:
!    if (x1 > dx ) then
   x1 = x1-dble(ie-1)*dx
!    endif
   x1 = 2.0d0*(x1/dx)-1.0d0
!    if (x2 > dx ) then    ! removed MT 1/2009,see above
   x2 = x2-dble(je-1)*dx
!    endif
   x2 = 2.0d0*(x2/dx)-1.0d0
! coordinates within an element [-1,1]
   cart%x = x1
   cart%y = x2
   number = ie+(je-1)*ne+(face_no-1)*ne*ne
!
   end subroutine cube_facepoint_ne
!================================================
!  (Nair) Cube face index and local coordinates
!================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine cube_facepoint_unstructured(sphere, cart, number, elem)
!-------------------------------------------------------------------------------
   use coordinatesystems, only: cube_face_number_from_sphere, sphere2cubedsphere, change_coordinates, cube_face_number_from_cart
   implicit none
!
   type(element_t), intent(in   ), target :: elem(:)
   type(spherical_polar_t), intent(in   ) :: sphere
   type(cartesian2d_t), intent(  out) :: cart
   integer, intent(  out) :: number
   integer :: ii
   logical :: found, p
   type(cartesian3d_t) :: sphere_xyz
   type(cartesian2d_t) :: cube
!
   sphere_xyz = spherical_to_cart(sphere)
   p = .false.
   do ii = 1,nelemd
     found = point_inside_equiangular(elem(ii),sphere,sphere_xyz)
     if (found) then
       number = ii
!          cube = sphere2cubedsphere(sphere,elem(ii)%FaceNum)
!          found = find_parametric_by_roots(cart,elem(ii),cube)
!          if (.not. found) then
!             cart = parametric_coordinates(sphere,elem(ii)%corners,elem(ii)%FaceNum)
       cart = parametric_coordinates2(sphere,elem(ii))
!          endif
       exit
     endif
   enddo
!
   end subroutine cube_facepoint_unstructured
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine interp_init()
!
   type(quadrature_t) :: gp
!
   gp = gausslobatto(np)
   call interpolate_create(gp,interp_p)
!
   end subroutine interp_init
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine setup_latlon_interp(elem, interpdata, par)
!
! initialize interpolation data structures to interpolate to a lat-lon grid
!
!
!-------------------------------------------------------------------------------
   implicit none
!
   type(element_t), intent(in   ), target :: elem(:)
   type(parallel_t), intent(in   ) :: par
   type(interpdata_t), intent(  out) :: interpdata(:)
! local
   integer i, j, ii, count_total, n_interp, count_max
   integer ngrid, number, elem_num, plat
   integer countx, missing_pts, ierr
   integer :: npts_mult_claims, max_claims
   real(real_kind) :: dp, latdeg(nlat+1), clat(nlat+1), w(nlat+1), w_staggered(nlat)
   real(real_kind) :: clat_staggered(nlat), latdeg_st(nlat), err, err2
   type(spherical_polar_t) :: sphere
   type(cartesian2d_t) :: cart
   type(cartesian3d_t) :: sphere_xyz, sphere2_xyz
   type(quadrature_t) :: gp
   logical, save :: first_time = .true.
! Array to make sure each interp point is on exactly one process
   type(cartesian2d_t), allocatable :: cart_vec(:,:)
   integer :: k
   integer, allocatable :: global_elem_gid(:,:), local_elem_gid(:,:), local_elem_num(:,:)
! these arrays often are too large for stack, so lets make sure
! they go on the heap:
!
   allocate(local_elem_num(nlat,nlon))
   allocate(local_elem_gid(nlat,nlon))
   allocate(global_elem_gid(nlat,nlon))
   allocate(cart_vec(nlat,nlon))
   if (par%ismasterproc) then
     write(iulog,'(a,i4,a,i4,a) ') 'Initializing ',nlat,' x ',nlon,' lat-lon interpolation grid:'
   endif
   do ii = 1,nelemd
     interpdata(ii)%n_interp = 0 ! reset counter
   enddo
   if (first_time) then
     nullify(lat)
     nullify(gweight)
     nullify(lon)
     first_time = .false.
   endif
   if (associated(lat)) then
     if (size(lat)>0) deallocate(lat)
   endif
   if (associated(gweight)) then
     if (size(gweight)>0) deallocate(gweight)
   endif
   if (associated(lon)) then
     if (size(lon)>0) deallocate(lon)
   endif
   allocate(lat(nlat))
   allocate(gweight(nlat))
   allocate(lon(nlon))
   call interp_init()
   gweight = 0
   do i = 1,nlon
     lon(i) = 2*dd_pi*(i-1)/nlon
   enddo
   if (gridtype==1) then
     do j = 1, nlat
       lat(j) =-dd_pi/2+dd_pi*(j-1)/(nlat-1)
     enddo
       plat = nlat
   endif
   if (gridtype==2) then
     gp = gauss(nlat)
     do j = 1,nlat
       lat(j) = asin(gp%points(j))
       gweight(j) = gp%weights(j)
     enddo
   endif
   if (gridtype==3) then
     do j = 1, nlat
       lat(j) =-dd_pi/2+dd_pi*(j-.5d0)/nlat
     enddo
       plat = nlat+1
   endif
   if (gridtype==1.or.gridtype==3) then
! gridtype=1    plat=nlat    gweight(1:nlat)=w(1:plat)
! gridtype=3    plat=nlat+1  gweight(1:nlat)=w_staggered(1:plat-1)
! L-R dynamics uses a regular latitude distribution (not gausian).
! The algorithm below is a bastardized version of LSM: map.F.
     dp = 180d0/(plat-1)
     do j = 1,plat
       latdeg(j) =-90d0+(j-1)*dp
       clat(j) = latdeg(j)*dd_pi/180d0
     enddo
! Calculate latitudes for the staggered grid
     do j = 1,plat-1
       clat_staggered(j) =(clat(j)+clat(j+1))/2
       latdeg_st(j) = clat_staggered(j)*180d0/dd_pi
     enddo
! Weights are defined as cos(phi)*(delta-phi)
! For a sanity check,the sum of w across all lats should be 2,or 1 across
! half of the latitudes.
     do j = 2,plat-1
       w(j) = sin(clat_staggered(j))-sin(clat_staggered(j-1))
     enddo
     w(1) = sin(clat_staggered(1))+1
     w(plat) = w(1)
! with nlat=2048,this error was 4e-16
     if (abs(sum(w(1:plat))-2)>1e-8) then
       write(iulog,*) 'interpolate_mod:w weights do not sum to 2. sum = ',sum(w(1:plat))
       call abortpar(message = 'interpolate_mod:weights do not sum to 2.')
     endif
     dp = dd_pi/(plat-1)
     do j = 1,plat-1
       w_staggered(j) = sin(clat(j+1))-sin(clat(j))
     enddo
     if (abs(sum(w_staggered(1:plat-1))-2)>1e-8) then
       write(iulog,*) 'interpolate_mod:staggered weights do not sum to 2. sum = ',sum(w_staggered(1:plat-1))
       call abortpar(message = 'interpolate_mod:weights do not sum to 2.')
     endif
     if (gridtype==1) then
       gweight(1:nlat) = w(1:plat)
     endif
     if (gridtype==3) then
       gweight(1:nlat) = w_staggered(1:plat-1)
     endif
   endif
! go through once,counting the number of points on each element
   sphere%r = 1
   local_elem_num =-1
   local_elem_gid =-1
   global_elem_gid =-1
   err = 0
   do j = 1,nlat
     do i = 1, nlon
       sphere%lat = lat(j)
       sphere%lon = lon(i)
       number =-1
       if (meshusemeshfile) then
         call cube_facepoint_unstructured(sphere,cart,number,elem)
         if (number/=-1) then
! If points are outside element but within tolerance, move to boundary
           if (cart%x+1.0d0.le.0.0d0) cart%x =-1.0d0
           if (cart%x-1.0d0.ge.0.0d0) cart%x = 1.0d0
           if (cart%y+1.0d0.le.0.0d0) cart%y =-1.0d0
           if (cart%y-1.0d0.ge.0.0d0) cart%y = 1.0d0
           local_elem_num(j,i) = number
           local_elem_gid(j,i) = elem(number)%vertex%number
           cart_vec(j,i) = cart ! local element coordiante of interpolation point
         endif
       else
         call cube_facepoint_ne(sphere,ne,cart,number)
! the sphere point belongs to the element number on face = face_no.
! do I own this element?
         if (number/=-1) then
           do ii = 1, nelemd
             if (number==elem(ii)%vertex%number) then
               local_elem_gid(j,i) = number
               local_elem_num(j,i) = ii
               cart_vec(j,i) = cart ! local element coordinate found above
               exit
             endif
             enddo
         endif
       endif
       ii = local_elem_num(j,i)
       if (ii/=-1) then
! compute error: map 'cart' back to sphere and compare with original
! interpolation point:
         sphere2_xyz = spherical_to_cart(ref2sphere_double(cart%x,cart%y,elem(ii)%corners,elem(ii)%vertex%face_number))
         sphere_xyz = spherical_to_cart(sphere)
         err = max(err,distance(sphere2_xyz,sphere_xyz))
       endif
     enddo
       if (par%ismasterproc) then
         if ((mod(j, 64).eq.1).or.(j.eq.nlat)) then
           print*,'finished latitude ',j,' of ',nlat
         endif
         endif
   enddo
   err2 = err
#ifdef _MPI
   call mpi_allreduce(err,err2,1,par_double_precision,par_max,par%comm,ierr)
#endif
   if (par%ismasterproc) then
     write(iulog,'(a,e12.4) ') 'Max interpolation point search error:',err2
   endif
! if multile elements claim a interpolation point,take the one with largest gid:
   global_elem_gid = local_elem_gid
#ifdef _MPI
   call mpi_allreduce(local_elem_gid,global_elem_gid,nlat*nlon,par_integer,par_max,par%comm,ierr)
#endif
   missing_pts = 0
   do j = 1,nlat
     do i = 1, nlon
       if (global_elem_gid(j, i)==-1) then
         missing_pts = missing_pts+1
         if (par%ismasterproc) print*,'Error:point not claimed by any element j,i,lat(j),lon(i) = ',j,i,lat(j),lon(i)
       elseif (local_elem_gid(j,i)==global_elem_gid(j,i)) then
         ii = local_elem_num(j,i)
         interpdata(ii)%n_interp = interpdata(ii)%n_interp+1
       endif
       enddo
   enddo
   countx = maxval(interpdata(1:nelemd)%n_interp)
   count_max = countx
#ifdef _MPI
   call mpi_allreduce(countx,count_max,1,par_integer,par_max,par%comm,ierr)
#endif
   if (par%ismasterproc) then
     write(iulog,'(a,i6) ') 'Maximum number of interpolation points claimed by an element:',count_max
   endif
! allocate storage
   do ii = 1,nelemd
     ngrid = interpdata(ii)%n_interp
     if (interpdata(ii)%first_entry) then
       nullify(interpdata(ii)%interp_xy)
       nullify(interpdata(ii)%ilat)
       nullify(interpdata(ii)%ilon)
       interpdata(ii)%first_entry = .false.
     endif
     if (associated(interpdata(ii)%interp_xy)) then
       if (size(interpdata(ii)%interp_xy)>0) deallocate(interpdata(ii)%interp_xy)
     endif
     if (associated(interpdata(ii)%ilat)) then
       if (size(interpdata(ii)%ilat)>0) deallocate(interpdata(ii)%ilat)
     endif
     if (associated(interpdata(ii)%ilon)) then
       if (size(interpdata(ii)%ilon)>0) deallocate(interpdata(ii)%ilon)
     endif
     allocate(interpdata(ii)%interp_xy(ngrid))
     allocate(interpdata(ii)%ilat(ngrid))
     allocate(interpdata(ii)%ilon(ngrid))
     interpdata(ii)%n_interp = 0 ! reset counter
   enddo
   do j = 1,nlat
     do i = 1, nlon
       if (local_elem_gid(j, i)==global_elem_gid(j, i).and.local_elem_gid(j, i)/=-1) then
         ii = local_elem_num(j,i)
         ngrid = interpdata(ii)%n_interp+1
         interpdata(ii)%n_interp = ngrid
         interpdata(ii)%interp_xy(ngrid) = cart_vec(j,i)
         interpdata(ii)%ilon(ngrid) = i
         interpdata(ii)%ilat(ngrid) = j
       endif
       enddo
   enddo
! now lets compute the number of points that were claimed by
! more than one element:
   do j = 1,nlat
     do i = 1, nlon
       if (local_elem_gid(j, i)==-1) then
         local_elem_gid(j,i) = 0
       else
         local_elem_gid(j,i) = 1
       endif
       enddo
   enddo
   global_elem_gid = local_elem_gid
#ifdef _MPI
   call mpi_allreduce(local_elem_gid,global_elem_gid,nlat*nlon,par_integer,par_sum,par%comm,ierr)
#endif
   if (par%ismasterproc) then
     countx = 0
     do j = 1,nlat
       do i = 1, nlon
         if (global_elem_gid(j,i)>1) countx = countx+1
       enddo
     enddo
     npts_mult_claims = countx
     max_claims = maxval(global_elem_gid)
   endif
   if (par%ismasterproc) then
     print*,'Number of interpolation points claimed by more than one element:',npts_mult_claims
     print*,'max number of elements which claimed the same interpolation point:',max_claims
   endif
   deallocate(global_elem_gid)
   deallocate(local_elem_num)
   deallocate(local_elem_gid)
   deallocate(cart_vec)
! check if every point in interpolation grid was claimed by an element:
   if (missing_pts>0) then
     count_total = nlat*nlon
     if (par%ismasterproc) then
       write(iulog,"(3a,i4,a,i7,a,i5) ") "error:",__file__," ",__line__," count_total:",count_total," missing:",missing_pts
     endif
     call barrierpar(par = par)
     call abortpar(message = 'Error:interpolation points not claimed by any element')
   endif
!
   end subroutine setup_latlon_interp
! =======================================
! interpolate_scalar
!
! Interpolate a scalar field given in an element (fld_cube) to the points in
! interpdata%interp_xy(i), i=1 .. interpdata%n_interp.
!
! Note that it is possible the given element contains none of the interpolation points
! =======================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine interpolate_scalar2d(interpdata, fld_cube, npts, fld, fillvalue)
!
   integer :: npts
   real(real_kind) :: fld_cube(npts, npts) ! cube field
   real(real_kind) :: fld(:) ! field at new grid lat, lon coordinates
   type(interpdata_t) :: interpdata
   real(real_kind), intent(in   ), optional :: fillvalue
! Local variables
   type(interpolate_t), pointer :: interp ! interpolation structure
   integer :: ne
   integer :: i
   type(cartesian2d_t) :: cart
!
   if (npts==np) then
     interp=>interp_p
   else
     call abortpar(message = 'Error in interpolate_scalar():must be called with p or v grid data')
   endif
! Choice for Native (high-order) or Bilinear interpolations
!
   if (present(fillvalue)) then
     if (itype==0) then
       do i = 1, interpdata%n_interp
         fld(i) = interpolate_2d(interpdata%interp_xy(i),fld_cube,interp,npts,fillvalue)
       enddo
         elseif (itype==1) then
         do i = 1,interpdata%n_interp
           fld(i) = interpol_bilinear(interpdata%interp_xy(i),fld_cube,interp,npts,fillvalue)
         enddo
       endif
       else
         if (itype==0) then
           do i = 1, interpdata%n_interp
             fld(i) = interpolate_2d(interpdata%interp_xy(i),fld_cube,interp,npts)
           enddo
             elseif (itype==1) then
             do i = 1,interpdata%n_interp
               fld(i) = interpol_bilinear(interpdata%interp_xy(i),fld_cube,interp,npts)
             enddo
           endif
   endif
!
   end subroutine interpolate_scalar2d
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine interpolate_scalar3d(interpdata, fld_cube, npts, nlev, fld, fillvalue)
!
   integer, intent(in   ) :: npts, nlev
   real(real_kind) :: fld_cube(npts, npts, nlev) ! cube field
   real(real_kind) :: fld(:,:) ! field at new grid lat, lon coordinates
   type(interpdata_t) :: interpdata
   real(real_kind), intent(in   ), optional :: fillvalue
! Local variables
   type(interpolate_t), pointer :: interp ! interpolation structure
   integer :: ne
   integer :: i, k
   type(cartesian2d_t) :: cart
!
   if (npts==np) then
     interp=>interp_p
   else
     call abortpar(message = 'Error in interpolate_scalar():must be called with p or v grid data')
   endif
! Choice for Native (high-order) or Bilinear interpolations
!
   if (present(fillvalue)) then
     if (itype==0) then
       do k = 1, nlev
         do i = 1, interpdata%n_interp
           fld(i,k) = interpolate_2d(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts,fillvalue)
         enddo
         enddo
           elseif (itype==1) then
           do k = 1,nlev
             do i = 1, interpdata%n_interp
               fld(i,k) = interpol_bilinear(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts,fillvalue)
             enddo
           enddo
         endif
         else
           if (itype==0) then
             do k = 1, nlev
               do i = 1, interpdata%n_interp
                 fld(i,k) = interpolate_2d(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts)
               enddo
               enddo
                 elseif (itype==1) then
                 do k = 1,nlev
                   do i = 1, interpdata%n_interp
                     fld(i,k) = interpol_bilinear(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts)
                   enddo
                 enddo
               else
                 write(iulog,*) itype
                 call abortpar(message = "wrong interpolation type")
               endif
   endif
!
   end subroutine interpolate_scalar3d
! =======================================
! interpolate_vector
!
! Interpolate a vector field given in an element (fld_cube)
! to the points in interpdata%interp_xy(i), i=1 .. interpdata%n_interp.
!
! input_coords = 0    fld_cube given in lat-lon
! input_coords = 1    fld_cube given in contravariant
!
! Note that it is possible the given element contains none of the interpolation points
! =======================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine interpolate_vector2d(interpdata, elem, fld_cube, npts, fld, input_coords, fillvalue)
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: npts
   real(real_kind) :: fld_cube(npts, npts, 2) ! vector field
   real(real_kind) :: fld(:,:) ! field at new grid lat, lon coordinates
   type(interpdata_t) :: interpdata
   type(element_t), intent(in   ) :: elem
   real(real_kind), intent(in   ), optional :: fillvalue
   integer :: input_coords
! Local variables
   real(real_kind) :: fld_contra(npts, npts, 2) ! vector field
   type(interpolate_t), pointer :: interp ! interpolation structure
   real(real_kind) :: v1, v2
   real(real_kind) :: d(2, 2) ! derivative of gnomonic mapping
   real(real_kind) :: jj(2, 2), tmpd(2, 2) ! derivative of gnomonic mapping
   integer :: i, j
   type(cartesian2d_t) :: cart
!
   if (present(fillvalue)) then
     if (any(fld_cube==fillvalue)) then
       fld = fillvalue
       return
     endif
   endif
!
   if (input_coords==0) then
! convert to contra
     do j = 1, npts
       do i = 1, npts
! latlon->contra
         fld_contra(i,j,1) = elem%dinv(1,1,i,j)*fld_cube(i,j,1)+elem%dinv(1,2,i,j)*fld_cube(i,j,2)
         fld_contra(i,j,2) = elem%dinv(2,1,i,j)*fld_cube(i,j,1)+elem%dinv(2,2,i,j)*fld_cube(i,j,2)
       enddo
       enddo
       else
         fld_contra = fld_cube
   endif
!
   if (npts==np) then
     interp=>interp_p
   elseif (npts==np) then
     call abortpar(message = 'Error in interpolate_vector():input must be on velocity grid')
   endif
! Choice for Native (high-order) or Bilinear interpolations
!
   if (itype==0) then
     do i = 1, interpdata%n_interp
       fld(i,1) = interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,1),interp,npts)
       fld(i,2) = interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,2),interp,npts)
     enddo
       elseif (itype==1) then
       do i = 1,interpdata%n_interp
         fld(i,1) = interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,1),interp,npts)
         fld(i,2) = interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,2),interp,npts)
       enddo
     else
       write(iulog,*) itype
       call abortpar(message = "wrong interpolation type")
   endif
!
   do i = 1, interpdata%n_interp
! convert fld from contra->latlon
     call dmap(d,elem,interpdata%interp_xy(i)%x,interpdata%interp_xy(i)%y)
! convert fld from contra->latlon
     v1 = fld(i,1)
     v2 = fld(i,2)
     fld(i,1) = d(1,1)*v1+d(1,2)*v2
     fld(i,2) = d(2,1)*v1+d(2,2)*v2
   enddo
!
   end subroutine interpolate_vector2d
! =======================================
! interpolate_vector
!
! Interpolate a vector field given in an element (fld_cube)
! to the points in interpdata%interp_xy(i), i=1 .. interpdata%n_interp.
!
! input_coords = 0    fld_cube given in lat-lon
! input_coords = 1    fld_cube given in contravariant
!
! Note that it is possible the given element contains none of the interpolation points
! =======================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine interpolate_vector3d(interpdata, elem, fld_cube, npts, nlev, fld, input_coords, fillvalue)
!-------------------------------------------------------------------------------
   implicit none
!
   type(interpdata_t), intent(in   ) :: interpdata
   type(element_t), intent(in   ) :: elem
   integer, intent(in   ) :: npts, nlev
   real(real_kind), intent(in   ) :: fld_cube(npts, npts, 2, nlev) ! vector field
   real(real_kind), intent(  out) :: fld(:,:,:) ! field at new grid lat, lon coordinates
   real(real_kind), intent(in   ), optional :: fillvalue
   integer, intent(in   ) :: input_coords
! Local variables
   real(real_kind) :: fld_contra(npts, npts, 2, nlev) ! vector field
   type(interpolate_t), pointer :: interp ! interpolation structure
   real(real_kind) :: v1, v2
   real(real_kind) :: d(2, 2) ! derivative of gnomonic mapping
   real(real_kind) :: jj(2, 2), tmpd(2, 2) ! derivative of gnomonic mapping
   integer :: i, j, k
   type(cartesian2d_t) :: cart
!
   if (present(fillvalue)) then
     if (any(fld_cube==fillvalue)) then
       fld = fillvalue
       return
     endif
   endif
!
   if (input_coords==0) then
! convert to contra
     do k = 1, nlev
       do j = 1, npts
         do i = 1, npts
! latlon->contra
           fld_contra(i,j,1,k) = elem%dinv(1,1,i,j)*fld_cube(i,j,1,k)+elem%dinv(1,2,i,j)*fld_cube(i,j,2,k)
           fld_contra(i,j,2,k) = elem%dinv(2,1,i,j)*fld_cube(i,j,1,k)+elem%dinv(2,2,i,j)*fld_cube(i,j,2,k)
         enddo
         enddo
         enddo
         else
           fld_contra = fld_cube
   endif
!
   if (npts==np) then
     interp=>interp_p
   elseif (npts==np) then
     call abortpar(message = 'Error in interpolate_vector():input must be on velocity grid')
   endif
! Choice for Native (high-order) or Bilinear interpolations
!
   if (itype==0) then
     do k = 1, nlev
       do i = 1, interpdata%n_interp
         fld(i,k,1) = interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,1,k),interp,npts)
         fld(i,k,2) = interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,2,k),interp,npts)
       enddo
       enddo
         elseif (itype==1) then
         do k = 1,nlev
           do i = 1, interpdata%n_interp
             fld(i,k,1) = interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,1,k),interp,npts)
             fld(i,k,2) = interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,2,k),interp,npts)
           enddo
         enddo
       else
         call abortpar(message = "wrong interpolation type")
   endif
!
   do i = 1, interpdata%n_interp
! compute D(:,:) at the point elem%interp_cube(i)
     call dmap(d,elem,interpdata%interp_xy(i)%x,interpdata%interp_xy(i)%y)
     do k = 1,nlev
! convert fld from contra->latlon
       v1 = fld(i,k,1)
       v2 = fld(i,k,2)
       fld(i,k,1) = d(1,1)*v1+d(1,2)*v2
       fld(i,k,2) = d(2,1)*v1+d(2,2)*v2
     enddo
   enddo
!
   end subroutine interpolate_vector3d
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function var_is_vector_uvar(name)
!
   character(len=*), intent(in   ) :: name
   integer :: i, var_is_vector_uvar, null_index
!
   var_is_vector_uvar = 0
   null_index = 0
   do i = 1,max_vecvars
     if (trim(vector_uvars(i)).eq.'') then
       null_index = i
       exit
     endif
       if (trim(vector_uvars(i)).eq.name) then
         var_is_vector_uvar = i
         exit
       endif
   enddo
#if 0
   disabled:breaks in many cases,like uv
   if (var_is_vector_uvar==0) then
! default rules: if variable starts with U and was not found, add it to the list:
     if (name(1:1).eq.'u') then
       if (null_index==0) then
         call abortpar(message = "error:max_vecvars too small")
       endif
         vector_uvars(null_index) = name
         vector_vvars(null_index) = name
         vector_vvars(null_index)(1:1) = 'V'
         var_is_vector_uvar = null_index
       endif
   endif
#endif
!
   end function var_is_vector_uvar
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function var_is_vector_vvar(name)
!
   character(len=*), intent(in   ) :: name
   integer :: i, var_is_vector_vvar, null_index
!
   var_is_vector_vvar = 0
   null_index = 0
   do i = 1,max_vecvars
     if (trim(vector_vvars(i)).eq.'') then
       null_index = i
       exit
     endif
       if (trim(vector_vvars(i)).eq.name) then
         var_is_vector_vvar = i
         exit
       endif
   enddo
#if 0
   disabled:breaks in many cases,like uv
   if (var_is_vector_vvar==0) then
! default rules: if variable starts with V and was not found, add it to the list:
     if (name(1:1).eq.'v') then
       if (null_index==0) then
         call abortpar(message = "error:max_vecvars too small")
       endif
         vector_uvars(null_index) = name
         vector_uvars(null_index)(1:1) = 'U'
         vector_vvars(null_index) = name
         var_is_vector_vvar = null_index
       endif
   endif
#endif
!
   end function var_is_vector_vvar
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module interpolate
!-------------------------------------------------------------------------------
