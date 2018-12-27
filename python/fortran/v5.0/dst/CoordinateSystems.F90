#include "KIM.h"
!-------------------------------------------------------------------------------
   module coordinatesystems
! WARNING:  When using this class be sure that you know if the
! cubic coordinates are on the unit cube or the [-\pi/4,\pi/4] cube
! and if the spherical longitude is in [0,2\pi] or [-\pi,\pi]
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
   use kiapsbase, only: int_kind=>kim_int_kind, real_kind=>kim_real8_kind, longdouble_kind=>kim_longdouble_kind
   implicit none
!
   private
   real(real_kind), public, parameter :: dist_threshold = 1.0d-9
   real(real_kind), parameter :: one = 1.0d0, two = 2.0d0
!
   type, public :: cartesian2d_t
     sequence
     real(real_kind) :: x ! x coordinate
     real(real_kind) :: y ! y coordinate
   end type cartesian2d_t
!
   type, public :: cartesian3d_t
     sequence
     real(real_kind) :: x ! x coordinate
     real(real_kind) :: y ! y coordinate
     real(real_kind) :: z ! z coordinate
   end type cartesian3d_t
!
   type, public :: spherical_polar_t
     sequence
     real(real_kind) :: r ! radius
     real(real_kind) :: lon ! longitude
     real(real_kind) :: lat ! latitude
   end type spherical_polar_t
!
   interface assignment(=)
     module procedure copy_cart2d
   end interface
!
   interface operator(==)
     module procedure eq_cart2d
   end interface
!
   interface distance
     module procedure distance_cart2d
     module procedure distance_cart2d_v
     module procedure distance_cart3d
     module procedure distance_cart3d_v
   end interface
!
   interface change_coordinates
     module procedure spherical_to_cart_v
     module procedure spherical_to_cart
     module procedure cart_to_spherical_v
     module procedure cart_to_spherical
     module procedure aray_to_spherical
   end interface
! we cant make this an interfaced because of real_kind=longdouble_kind,
! these subroutine are the same
!  interface ref2sphere
!     module procedure ref2sphere_quad
!     module procedure ref2sphere_double
!  end interface
! ==========================================
! Public Interfaces
! ==========================================
   public :: sphere_tri_area
   public :: surfareaxy
   public :: distance
   public :: change_coordinates
   public :: ref2sphere
   public :: ref2sphere_double
   public :: cart2cubedsphere
   public :: spherical_to_cart !ce
! note: cant make these next two an interface since they only differ by return arg
   public :: cubedsphere2cart
   public :: sphere2cubedsphere
   public :: cube_face_number_from_cart
   public :: cube_face_number_from_sphere
! CE
   public :: cart2cubedspherexy
   public :: cart2spherical !ce
   private :: copy_cart2d
   private :: eq_cart2d
   private :: distance_cart2d
   private :: distance_cart2d_v
   private :: distance_cart3d
   private :: distance_cart3d_v
   private :: spherical_to_cart_v
!private :: spherical_to_cart
   private :: cart_to_spherical_v
   private :: cart_to_spherical
   private :: aray_to_spherical
!
   contains
! ============================================
! copy_cart2d:
!
! Overload assignment operator for cartesian2D_t
! ============================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine copy_cart2d(cart2, cart1)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian2d_t), intent(  out) :: cart2
   type(cartesian2d_t), intent(in   ) :: cart1
!
   cart2%x = cart1%x
   cart2%y = cart1%y
!
   end subroutine copy_cart2d
! ============================================
! eq_cart2d:
!
! Overload == operator for cartesian2D_t
! ============================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function eq_cart2d(cart2, cart1) result(is_same)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian2d_t), intent(in   ) :: cart2
   type(cartesian2d_t), intent(in   ) :: cart1
   logical :: is_same
!
   if (distance(cart1, cart2)<dist_threshold) then
     is_same = .true.
   else
     is_same = .false.
   endif
!
   end function eq_cart2d
! ===================================================
! distance_cart2D  : scalar version
! distance_cart2D_v: vector version
!
! computes distance between cartesian 2D coordinates
! ===================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function distance_cart2d(cart1, cart2) result(dist)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian2d_t), intent(in   ) :: cart1
   type(cartesian2d_t), intent(in   ), optional :: cart2
   real(real_kind) :: dist
!
   if (present(cart2)) then
     dist = sqrt((cart1%x-cart2%x)**2+&
                                                 (cart1%y-cart2%y)**2)
   else
     dist = sqrt(cart1%x*cart1%x+&
                                                            cart1%y*cart1%y)
   endif
!
   end function distance_cart2d
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function distance_cart2d_v(cart1, cart2) result(dist)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian2d_t), intent(in   ) :: cart1(:)
   type(cartesian2d_t), intent(in   ), optional :: cart2(:)
   real(real_kind) :: dist(size(cart1))
   integer :: i
!
   if (present(cart2)) then
     forall(i = 1:size(cart1)) dist(i) = distance_cart2d(cart1(i),cart2(i))
   else
     forall(i = 1:size(cart1)) dist(i) = distance_cart2d(cart1(i))
   endif
!
   end function distance_cart2d_v
! ===================================================
! distance_cart3D  : scalar version
! distance_cart3D_v: vector version
! ===================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function distance_cart3d(cart1, cart2) result(dist)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian3d_t), intent(in   ) :: cart1
   type(cartesian3d_t), intent(in   ), optional :: cart2
   real(real_kind) :: dist
!
   if (present(cart2)) then
     dist = sqrt((cart1%x-cart2%x)**2+&
                                                 (cart1%y-cart2%y)**2+&
                                                 (cart1%z-cart2%z)**2)
   else
     dist = sqrt(cart1%x*cart1%x+&
                                                            cart1%y*cart1%y+&
                                                            cart1%z*cart1%z)
   endif
!
   end function distance_cart3d
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function distance_cart3d_v(cart1, cart2) result(dist)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian3d_t), intent(in   ) :: cart1(:)
   type(cartesian3d_t), intent(in   ), optional :: cart2(:)
   real(real_kind) :: dist(size(cart1))
   integer :: i
!
   if (present(cart2)) then
     forall(i = 1:size(cart1)) dist(i) = distance_cart3d(cart1(i),cart2(i))
   else
     forall(i = 1:size(cart1)) dist(i) = distance_cart3d(cart1(i))
   endif
!
   end function distance_cart3d_v
! ===================================================================
! spherical_to_cart:
! converts spherical polar {lon,lat}  to 3D cartesian {x,y,z}
! on unit sphere.  Note: spherical longitude is [0,2\pi]
! ===================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function spherical_to_cart(sphere) result(cart)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spherical_polar_t), intent(in   ) :: sphere
   type(cartesian3d_t) :: cart
!
   cart%x = sphere%r*cos(sphere%lat)*cos(sphere%lon)
   cart%y = sphere%r*cos(sphere%lat)*sin(sphere%lon)
   cart%z = sphere%r*sin(sphere%lat)
!
   end function spherical_to_cart
! ===================================================================
! spherical_to_cart_v:
! converts spherical polar {lon,lat}  to 3D cartesian {x,y,z}
! on unit sphere.  Note: spherical longitude is [0,2\pi]
! ===================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function spherical_to_cart_v(sphere) result(cart)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spherical_polar_t), intent(in   ) :: sphere(:)
   type(cartesian3d_t) :: cart(size(sphere))
   integer :: i
!
   forall(i = 1:size(sphere)) cart(i) = spherical_to_cart(sphere(i))
!
   end function spherical_to_cart_v
! ==========================================================================
! cart_to_spherical:
!
! converts 3D cartesian {x,y,z} to spherical polar {lon,lat}
! on unit sphere. Note: spherical longitude is [0,2\pi]
! ==========================================================================
! scalar version
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function cart_to_spherical(cart) result(sphere)
!-------------------------------------------------------------------------------
   use physicalconstants, only: dd_pi
   implicit none
!
   type(cartesian3d_t), intent(in   ) :: cart
   type(spherical_polar_t) :: sphere
!
   sphere%r = distance(cart)
   sphere%lat = asin(cart%z/sphere%r)
   sphere%lon = 0
! ==========================================================
! enforce three facts:
!
! 1) lon at poles is defined to be zero
!
! 2) Grid points must be separated by about .01 Meter (on earth)
!    from pole to be considered "not the pole".
!
! 3) range of lon is { 0<= lon < 2*pi }
!
! ==========================================================
!   if point is away from the POLE.  distance(cart) = distance from center of earth,
!   so this was a bug:
!    if (distance(cart) >= DIST_THRESHOLD) then
   if (abs(abs(sphere%lat)-dd_pi/2)>= dist_threshold) then
     sphere%lon = atan2(cart%y,cart%x)
     if (sphere%lon<0) then
       sphere%lon = sphere%lon+2*dd_pi
     endif
   endif
!
   end function cart_to_spherical
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function aray_to_spherical(coordinates) result(sphere)
!-------------------------------------------------------------------------------
   implicit none
!
   real(real_kind), intent(in   ) :: coordinates(3)
   type(spherical_polar_t) :: sphere
   type(cartesian3d_t) :: cart
!
   cart%x = coordinates(1)
   cart%y = coordinates(2)
   cart%z = coordinates(3)
   sphere = cart_to_spherical(cart)
!
   end function aray_to_spherical
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function cart_to_spherical_v(cart) result(sphere)
!-------------------------------------------------------------------------------
   use physicalconstants, only: dd_pi
   implicit none
!
   type(cartesian3d_t), intent(in   ) :: cart(:)
   type(spherical_polar_t) :: sphere(size(cart))
   integer :: i
!
   forall(i = 1:size(cart)) sphere(i) = cart_to_spherical(cart(i))
!
   end function cart_to_spherical_v
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function unit_face_based_cube_to_unit_sphere(cart, face_no) result(sphere)
! Note: Output spherical longitude is [-pi,pi]
! Project from a UNIT cube to a UNIT sphere.  ie, the lenght of the cube edge is 2.
! Face 1 of the cube touches the sphere at longitude, latitude (0,0). The negative
! x axis is negative longitude (ie. going west is negative), the positive x axis
! is increasing longitude.  Face 1 maps the Face 1 to the lat,lon on the sphere:
!    [-1,1] x [-1,1] => [-\pi/4,\pi/4] x [-\pi/4, \pi/4]
! Face 2 continues with increasing longitude (ie to the east of Face 1).
! The left edge of Face 2 (negative x) is the right edge of Face 1 (positive x)
! The latitude is the same as Face 1, but the longitude increases:
!    [-1,1] x [-1,1] => [\pi/4, 3\pi/4] x [-\pi/4, \pi/4]
! Face 3 continues with increasing longitude (ie to the east of Face 2).
! Face 3 is like Face 1, but the x coordinates are reversed, ie. decreasing x
! is increasing longitude:
!    [-1,1] x [-1,1]  =    [-1,0] x [-1,1] U  [0,1] x [-1,1] =>
!            [3\pi/4,\pi] x [-\pi, -3\pi/4]
! Face 4 finally connects Face 3 to Face 1.  Like Face 2, but wtih opposite x
!    [-1,1] x [-1,1] => [-3\pi/4, -\pi/4] x [-\pi/4, \pi/4]
! Face 5 is along the bottom edges of Faces 1,2,3,and 4 so the latitude goes from
! -\pi/4 to -\pi/2.  The tricky part is lining up the longitude.  The zero longitude
! must line up with the center of Face 1. ATAN2(x,1) = 0 => x = 0.
! So the (0,1) point on Face 5 is the zero longitude on the sphere.  The top edge of
! Face 5 is the bottom edge of Face 1.
! ATAN(x,0) = \pi/2 => x = 1, so the right edge of Face 5 is the bottom of Face 2.
! Continueing, the bottom edge of 5 is the bottom of 3.  Left of 5 is bottom of 4.
! Face 6 is along the top edges of Faces 1,2,3 and 4 so the latitude goes from
! \pi/4 to \pi/2.   The zero longitude must line up with the center of Face 1.
! This is just like Face 5, but the y axis is reversed.  So the bottom edge of Face 6
! is the top edge of Face 1.  The right edge of Face 6 is the top of Face 2.  The
! top of 6 the top of 3 and the left of 6 the top of 4.
!-------------------------------------------------------------------------------
   use physicalconstants, only: dd_pi
!    use parallel_mod,       only : abortmp
   use kiapsparallel, only: abortpar
   implicit none
!
   type(cartesian2d_t), intent(in   ) :: cart ! on face_no of a unit cube
   integer, intent(in   ) :: face_no
   type(spherical_polar_t) :: sphere
   integer i, j
   real(real_kind) :: r!, l_inf
! MNL: removing check that points are on the unit cube because we allow spherical grids to map beyond the extent of the cube (though we probably should still have an upper bound for how far past the edge the element lies...)
!    l_inf = MAX(ABS(cart%x), ABS(cart%y))
!    if (1.01 < l_inf) then
!      call abortmp('unit_face_based_cube_to_unit_sphere: Input not on unit cube.')
!    end if
!
   sphere%r = one
   r = sqrt(one+(cart%x)**2+(cart%y)**2)
   select case(face_no)
   case(1)
     sphere%lat = asin((cart%y)/r)
     sphere%lon = atan2(cart%x,one)
   case(2)
     sphere%lat = asin((cart%y)/r)
     sphere%lon = atan2(one,-cart%x)
   case(3)
     sphere%lat = asin((cart%y)/r)
     sphere%lon = atan2(-cart%x,-one)
   case(4)
     sphere%lat = asin((cart%y)/r)
     sphere%lon = atan2(-one,cart%x)
   case(5)
     if (abs(cart%y)>dist_threshold.or.abs(cart%x)>dist_threshold) then
       sphere%lon = atan2(cart%x,cart%y)
     else
       sphere%lon = 0.0d0 ! longitude is meaningless at south pole set to 0.0
     endif
       sphere%lat = asin(-one/r)
     case(6)
       if (abs(cart%y)>dist_threshold.or.abs(cart%x)>dist_threshold) then
         sphere%lon = atan2(cart%x,-cart%y)
       else
         sphere%lon = 0.0d0 ! longitude is meaningless at north pole set to 0.0
       endif
         sphere%lat = asin(one/r)
       case default
!call abortmp('unit_face_based_cube_to_unit_sphere: Face number not 1 to 6.')
         call abortpar(message = 'unit_face_based_cube_to_unit_sphere:face number not 1 to 6.')
   end select
   if (sphere%lon<0.0d0) then
     sphere%lon = sphere%lon+two*dd_pi
   endif
!
   end function unit_face_based_cube_to_unit_sphere
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function cart2spherical(x, y, face_no) result(sphere)
! IMPORTANT: INPUT ARE the REAL cartesian from the cube sphere
! Note: Output spherical longitude is [-pi,pi]
! Project from a UNIT cube to a UNIT sphere.  ie, the lenght of the cube edge is 2.
! Face 1 of the cube touches the sphere at longitude, latitude (0,0). The negative
! x axis is negative longitude (ie. going west is negative), the positive x axis
! is increasing longitude.  Face 1 maps the Face 1 to the lat,lon on the sphere:
!    [-1,1] x [-1,1] => [-\pi/4,\pi/4] x [-\pi/4, \pi/4]
! Face 2 continues with increasing longitude (ie to the east of Face 1).
! The left edge of Face 2 (negative x) is the right edge of Face 1 (positive x)
! The latitude is the same as Face 1, but the longitude increases:
!    [-1,1] x [-1,1] => [\pi/4, 3\pi/4] x [-\pi/4, \pi/4]
! Face 3 continues with increasing longitude (ie to the east of Face 2).
! Face 3 is like Face 1, but the x coordinates are reversed, ie. decreasing x
! is increasing longitude:
!    [-1,1] x [-1,1]  =    [-1,0] x [-1,1] U  [0,1] x [-1,1] =>
!            [3\pi/4,\pi] x [-\pi, -3\pi/4]
! Face 4 finally connects Face 3 to Face 1.  Like Face 2, but wtih opposite x
!    [-1,1] x [-1,1] => [-3\pi/4, -\pi/4] x [-\pi/4, \pi/4]
! Face 5 is along the bottom edges of Faces 1,2,3,and 4 so the latitude goes from
! -\pi/4 to -\pi/2.  The tricky part is lining up the longitude.  The zero longitude
! must line up with the center of Face 1. ATAN2(x,1) = 0 => x = 0.
! So the (0,1) point on Face 5 is the zero longitude on the sphere.  The top edge of
! Face 5 is the bottom edge of Face 1.
! ATAN(x,0) = \pi/2 => x = 1, so the right edge of Face 5 is the bottom of Face 2.
! Continueing, the bottom edge of 5 is the bottom of 3.  Left of 5 is bottom of 4.
! Face 6 is along the top edges of Faces 1,2,3 and 4 so the latitude goes from
! \pi/4 to \pi/2.   The zero longitude must line up with the center of Face 1.
! This is just like Face 5, but the y axis is reversed.  So the bottom edge of Face 6
! is the top edge of Face 1.  The right edge of Face 6 is the top of Face 2.  The
! top of 6 the top of 3 and the left of 6 the top of 4.
!-------------------------------------------------------------------------------
   use physicalconstants, only: dd_pi
!    use parallel_mod,       only : abortmp
   use kiapsparallel, only: abortpar
   implicit none
!
   real(real_kind), intent(in   ) :: x, y ! on face_no of a unit cube
   integer, intent(in   ) :: face_no
   type(spherical_polar_t) :: sphere
   integer i, j
   real(real_kind) :: r!, l_inf
! MNL: removing check that points are on the unit cube because we allow spherical grids to map beyond the extent of the cube (though we probably should still have an upper bound for how far past the edge the element lies...)
!    l_inf = MAX(ABS(cart%x), ABS(cart%y))
!    if (1.01 < l_inf) then
!      call abortmp('unit_face_based_cube_to_unit_sphere: Input not on unit cube.')
!    end if
!
   sphere%r = one
   r = sqrt(one+x**2+y**2)
   select case(face_no)
   case(1)
     sphere%lat = asin(y/r)
     sphere%lon = atan2(x,one)
   case(2)
     sphere%lat = asin(y/r)
     sphere%lon = atan2(one,-x)
   case(3)
     sphere%lat = asin(y/r)
     sphere%lon = atan2(-x,-one)
   case(4)
     sphere%lat = asin(y/r)
     sphere%lon = atan2(-one,x)
   case(5)
     if (abs(y)>dist_threshold.or.abs(x)>dist_threshold) then
       sphere%lon = atan2(x,y)
     else
       sphere%lon = 0.0d0 ! longitude is meaningless at south pole set to 0.0
     endif
       sphere%lat = asin(-one/r)
     case(6)
       if (abs(y)>dist_threshold.or.abs(x)>dist_threshold) then
         sphere%lon = atan2(x,-y)
       else
         sphere%lon = 0.0d0 ! longitude is meaningless at north pole set to 0.0
       endif
         sphere%lat = asin(one/r)
       case default
!call abortmp('unit_face_based_cube_to_unit_sphere: Face number not 1 to 6.')
         call abortpar(message = 'unit_face_based_cube_to_unit_sphere:face number not 1 to 6.')
   end select
   if (sphere%lon<0.0d0) then
     sphere%lon = sphere%lon+two*dd_pi
   endif
!
   end function cart2spherical
#if 0
! Note: Output spherical longitude is [-pi,pi]
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine project(sphere, cart, face_no)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian2d_t), intent(in   ) :: cart(:,:) ! assumed to be cartesian coordinates of cube
   integer, intent(in   ) :: face_no
   type(spherical_polar_t), intent(  out) :: sphere(:,:)
   integer npts, i, j
!
   npts = size(cart,1)
   do j = 1,npts
     do i = 1, npts
       sphere(i,j) = unit_face_based_cube_to_unit_sphere(cart(i,j),face_no)
     enddo
   enddo
!
   end subroutine project
#endif
!
! map a point in the referece element to the sphere
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function ref2sphere_double(a, b, corners, face_no) result(sphere)
!-------------------------------------------------------------------------------
   implicit none
!
   real(real_kind) :: a, b
   integer, intent(in   ) :: face_no
   type(spherical_polar_t) :: sphere
   type(cartesian2d_t) :: corners(4)
! local
   real(real_kind) :: pi, pj, qi, qj
   type(cartesian2d_t) :: cart
! map (a,b) to the [-pi/2,pi/2] equi angular cube face:  x1,x2
! a = gp%points(i)
! b = gp%points(j)
!
   pi =(1-a)/2
   pj =(1-b)/2
   qi =(1+a)/2
   qj =(1+b)/2
   cart%x = pi*pj*corners(1)%x &
                                                 +qi*pj*corners(2)%x &
                                                 +qi*qj*corners(3)%x &
                                                 +pi*qj*corners(4)%x
   cart%y = pi*pj*corners(1)%y &
                                                 +qi*pj*corners(2)%y &
                                                 +qi*qj*corners(3)%y &
                                                 +pi*qj*corners(4)%y
! map from [pi/2,pi/2] equ angular cube face to sphere:
   sphere = projectpoint(cart,face_no)
!
   end function ref2sphere_double
!
! map a point in the referece element to the sphere
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function ref2sphere(a, b, corners, face_no) result(sphere)
!-------------------------------------------------------------------------------
   implicit none
!
   real(longdouble_kind) :: a, b
   integer, intent(in   ) :: face_no
   type(spherical_polar_t) :: sphere
   type(cartesian2d_t) :: corners(4)
! local
   real(real_kind) :: pi, pj, qi, qj
   type(cartesian2d_t) :: cart
! map (a,b) to the [-pi/2,pi/2] equi angular cube face:  x1,x2
! a = gp%points(i)
! b = gp%points(j)
!
   pi =(1-a)/2
   pj =(1-b)/2
   qi =(1+a)/2
   qj =(1+b)/2
   cart%x = pi*pj*corners(1)%x &
                                                 +qi*pj*corners(2)%x &
                                                 +qi*qj*corners(3)%x &
                                                 +pi*qj*corners(4)%x
   cart%y = pi*pj*corners(1)%y &
                                                 +qi*pj*corners(2)%y &
                                                 +qi*qj*corners(3)%y &
                                                 +pi*qj*corners(4)%y
! map from [pi/2,pi/2] equ angular cube face to sphere:
   sphere = projectpoint(cart,face_no)
!
   end function ref2sphere
! Note: Output spherical longitude is [-pi,pi]
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function projectpoint(cartin, face_no) result(sphere)
! Projection from a [-pi/4, \pi/4] sized cube.
! This will be checked because unit_face_based_cube_to_unit_sphere checks the ranges.
! See unit_face_based_cube_to_unit_sphere for documentation.
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian2d_t), intent(in   ) :: cartin
   integer, intent(in   ) :: face_no
   type(spherical_polar_t) :: sphere
   type(cartesian2d_t) :: cart
!ASC  This is X and Y and not xhi eta ...
!
   cart%x = tan(cartin%x)
   cart%y = tan(cartin%y)
   sphere = unit_face_based_cube_to_unit_sphere(cart,face_no)
!
   end function projectpoint
! takes a 2D point on a face of the cube of size [-\pi/4, \pi/4] and projects it
! onto a 3D point on a cube of size [-1,1] in R^3
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function cubedsphere2cart(cartin, face_no) result(cart)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian2d_t), intent(in   ) :: cartin ! assumed to be cartesian coordinates of cube
   integer, intent(in   ) :: face_no
   type(cartesian3d_t) :: cart
!
   cart = spherical_to_cart(projectpoint(cartin,face_no))
!
   end function cubedsphere2cart
! onto a cube of size [-\pi/2,\pi/2] in R^3
! the spherical longitude can be either in [0,2\pi] or [-\pi,\pi]
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function sphere2cubedsphere(sphere, face_no) result(cart)
!-------------------------------------------------------------------------------
   use physicalconstants, only: dd_pi
   implicit none
!
   type(spherical_polar_t), intent(in   ) :: sphere
   integer, intent(in   ) :: face_no
   type(cartesian2d_t) :: cart
   real(real_kind) :: xp, yp
   real(real_kind) :: lat, lon
   real(real_kind) :: pi, twopi, pi2, pi3, pi4
!
   lat = sphere%lat
   lon = sphere%lon
   pi = dd_pi
   twopi = 2.0d0*pi
   pi2 = pi*0.5d0
   pi3 = pi*1.5d0
   pi4 = pi*0.25d0
   select case(face_no)
   case(1)
     xp = lon
     if (pi<lon) xp = lon-twopi !if lon in [0,2\pi]
     yp = atan(tan(lat)/cos(xp))
   case(2)
     xp = lon-pi2
     yp = atan(tan(lat)/cos(xp))
   case(3)
     xp = lon-pi
     if (lon<0) xp = lon+pi !if lon in [0,2\pi]
     yp = atan(tan(lat)/cos(xp))
   case(4)
     xp = lon-pi3
     if (lon<0) xp = lon+pi2 !if lon in [0,2\pi]
     yp = atan(tan(lat)/cos(xp))
   case(5)
     xp = atan(-sin(lon)/tan(lat))
     yp = atan(-cos(lon)/tan(lat))
   case(6)
     xp = atan(sin(lon)/tan(lat))
     yp = atan(-cos(lon)/tan(lat))
   end select
! coordinates on the cube:
   cart%x = xp
   cart%y = yp
!
   end function sphere2cubedsphere
! Go from an arbitrary sized cube in 3D
! to a [-\pi/4,\pi/4] sized cube with (face,2d) coordinates.
!
!                        Z
!                        |
!                        |
!                        |
!                        |
!                        ---------------Y
!                       /
!                      /
!                     /
!                    /
!                   X
!
! NOTE: Face 1 =>  X positive constant face of cube
!       Face 2 =>  Y positive constant face of cube
!       Face 3 =>  X negative constant face of cube
!       Face 4 =>  Y negative constant face of cube
!       Face 5 =>  Z negative constant face of cube
!       Face 6 =>  Z positive constant face of cube
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function cart2cubedsphere(cart3d, face_no) result(cart)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian3d_t), intent(in   ) :: cart3d
   integer, intent(in   ) :: face_no
   type(cartesian2d_t) :: cart
   real(real_kind) :: x, y
!
   select case(face_no)
   case(1)
     x = cart3d%y/cart3d%x
     y = cart3d%z/cart3d%x
   case(2)
     x =-cart3d%x/cart3d%y
     y = cart3d%z/cart3d%y
   case(3)
     x = cart3d%y/cart3d%x
     y =-cart3d%z/cart3d%x
   case(4)
     x =-cart3d%x/cart3d%y
     y =-cart3d%z/cart3d%y
   case(5)
     x =-cart3d%y/cart3d%z
     y =-cart3d%x/cart3d%z
   case(6)
     x = cart3d%y/cart3d%z
     y =-cart3d%x/cart3d%z
   end select
!
   cart%x = atan(x)
   cart%y = atan(y)
!
   end function cart2cubedsphere
! This function divides three dimentional space up into
! six sectors.  These sectors are then considered as the
! faces of the cube.  It should work for any (x,y,z) coordinate
! if on a sphere or on a cube.
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function cube_face_number_from_cart(cart) result(face_no)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cartesian3d_t), intent(in   ) :: cart
   integer :: face_no
   real(real_kind) :: x, y, z
!
   x = cart%x
   y = cart%y
   z = cart%z
! Divide the X-Y plane into for quadrants of
! [-\pi/2,\pi/2],[\pi/2,3\pi/2],.....
! based on the lines X=Y and X=-Y.  This divides
! 3D space up into four sections.  Doing the same
! for the XZ and YZ planes divides space into six
! sections.  Can also be thought of as conic sections
! in the L_infinity norm.
   if (y<x.and.y>-x) then ! x>0,face 1,5 or 6
     if (z>x) then
       face_no = 6 ! north pole
     elseif (z<-x) then
       face_no = 5 ! south pole
     else
       face_no = 1
     endif
     elseif (y>x.and.y<-x) then ! x<0
       if (z>-x) then
         face_no = 6 ! north pole
       elseif (z<x) then
         face_no = 5 ! south pole
       else
         face_no = 3
       endif
       elseif (y>x.and.y>-x) then ! y>0
         if (z>y) then
           face_no = 6 ! north pole
         elseif (z<-y) then
           face_no = 5 ! south pole
         else
           face_no = 2
         endif
         elseif (y<x.and.y<-x) then ! y<0
           if (z>-y) then
             face_no = 6 ! north pole
           elseif (z<y) then
             face_no = 5 ! south pole
           else
             face_no = 4
           endif
           else
! abs(y) = abs(x).  point is on cube edge, or on face 5 or 6:
             if (abs(x)<z) then
               face_no = 6 ! north pole
             elseif (z<-abs(x)) then
               face_no = 5 ! south pole
             elseif (0<x.and.0<y) then
               face_no = 1
             elseif (x<0.and.0<y) then
               face_no = 2
             elseif (x<0.and.y<0) then
               face_no = 3
             else
               face_no = 4
             endif
   endif
!
   end function cube_face_number_from_cart
! This could be done directly by using the lon, lat coordinates,
! but call cube_face_number_from_cart just so that there is one place
! to do the conversions and they are all consistant.
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   pure function cube_face_number_from_sphere(sphere) result(face_no)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spherical_polar_t), intent(in   ) :: sphere
   type(cartesian3d_t) :: cart
   integer :: face_no
!
   cart = spherical_to_cart(sphere)
   face_no = cube_face_number_from_cart(cart)
!
   end function cube_face_number_from_sphere
! CE, need real (cartesian) xy coordinates on the cubed sphere
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine cart2cubedspherexy(cart3d, face_no, cartxy)
!
   type(cartesian3d_t), intent(in   ) :: cart3d
   integer, intent(in   ) :: face_no
   type(cartesian2d_t), intent(  out) :: cartxy
! a (half length of a cube side) is supposed to be 1
!
   select case(face_no)
   case(1)
     cartxy%x = cart3d%y/cart3d%x
     cartxy%y = cart3d%z/cart3d%x
   case(2)
     cartxy%x =-cart3d%x/cart3d%y
     cartxy%y = cart3d%z/cart3d%y
   case(3)
     cartxy%x = cart3d%y/cart3d%x
     cartxy%y =-cart3d%z/cart3d%x
   case(4)
     cartxy%x =-cart3d%x/cart3d%y
     cartxy%y =-cart3d%z/cart3d%y
   case(5) !bottom face
     cartxy%x =-cart3d%y/cart3d%z
     cartxy%y =-cart3d%x/cart3d%z
   case(6) !top face
     cartxy%x = cart3d%y/cart3d%z
     cartxy%y =-cart3d%x/cart3d%z
   end select
!
   end subroutine cart2cubedspherexy
! CE END
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine sphere_tri_area(v1, v2, v3, area)
!  input: v1(3),v2(3),v3(3)  cartesian coordinates of triangle
!  output: area
!  based on formulas in STRI_QUAD:
!  http://people.sc.fsu.edu/~burkardt/f_src/stri_quad/stri_quad.html
!-------------------------------------------------------------------------------
   use physicalconstants, only: dd_pi
   implicit none
!
   real(real_kind) area
   real(real_kind) a, b, c, al, bl, cl, sina, sinb, sinc, sins, a1, b1, c1
   type(cartesian3d_t) v1, v2, v3
! compute great circle lengths
!
   al = acos(v2%x*v3%x+v2%y*v3%y+v2%z*v3%z)
   bl = acos(v3%x*v1%x+v3%y*v1%y+v3%z*v1%z)
   cl = acos(v1%x*v2%x+v1%y*v2%y+v1%z*v2%z)
! compute angles
   sina = sin((bl+cl-al)/2) ! sin(sl-al)
   sinb = sin((al+cl-bl)/2) ! sin(sl-bl)
   sinc = sin((al+bl-cl)/2)
   sins = sin((al+bl+cl)/2)
#if 0
   a = 2*atan2(sqrt(sinb*sinc),sqrt(sins*sina))
   b = 2*atan2(sqrt(sina*sinc),sqrt(sins*sinb))
   c = 2*atan2(sqrt(sina*sinb),sqrt(sins*sinc))
! apply Girard's theorem
   area = a+b+c-dd_pi
#endif
! for small areas,formula above looses precision.
! 2atan(x) + 2atan(1/x) = pi
! 2atan(x) - pi = -2atan(1/x)
   a = sqrt((sinb*sinc)/(sins*sina))
   b = sqrt((sina*sinc)/(sins*sinb))
   c = sqrt((sina*sinb)/(sins*sinc))
   a1 = 2*atan(a)
   b1 = 2*atan(b)
   c1 = 2*atan(c)
   if (a.gt.b.and.a.gt.c) then
     a1 =-2*atan(1/a)
   elseif (b.gt.c) then
     b1 =-2*atan(1/b)
   else
     c1 =-2*atan(1/c)
   endif
! apply Girard's theorem
   area = a1+b1+c1
!
   end subroutine sphere_tri_area
!CE, 5.May 2011
!INPUT: Points in xy cubed sphere coordinates, counterclockwise
!OUTPUT: corresponding area on the sphere
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function surfareaxy(x1, x2, y1, y2) result(area)
!-------------------------------------------------------------------------------
   implicit none
!
   real(real_kind), intent(in   ) :: x1, x2, y1, y2
   real(real_kind) :: area
   real(real_kind) :: a1, a2, a3, a4
! cubed-sphere cell area, from Lauritzen & Nair MWR 2008
! central angles:
! cube face: -pi/4,-pi/4 -> pi/4,pi/4
! this formula gives 2   so normalize by 4pi/6 / 2 = pi/3
! use implementation where the nodes a counterclockwise (not as in the paper)
!
   a1 = acos(-sin(atan(x1))*sin(atan(y1)))
   a2 =-acos(-sin(atan(x2))*sin(atan(y1)))
   a3 = acos(-sin(atan(x2))*sin(atan(y2)))
   a4 =-acos(-sin(atan(x1))*sin(atan(y2)))
   area =(a1+a2+a3+a4)
   return
!
   end function surfareaxy
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module coordinatesystems
!-------------------------------------------------------------------------------
