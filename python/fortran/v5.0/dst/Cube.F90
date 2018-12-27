!-------------------------------------------------------------------------------
!
!  MODULE Cube
!
!> @brief
!>  - Cube module.
!>
!> @date ?????2012
!>  - Junghan Kim : First written from the HOMME and was modified for
!>                  KIAPSGM framework.
!> @date 30JAN2015
!>  - Junghan Kim : Added to some variables and subroutines. (CAM-SE 5.3)
!> @date 25FEB2015
!>  - In-Sun Song : Comments added and codes cleaned up
!
!-------------------------------------------------------------------------------
!
#include <KIM.h>
!
#define _BEGIN_FACE 1
#define _END_FACE   4
#undef _FACE_6
#undef _FACE_5
!
!-------------------------------------------------------------------------------
   module cube
!
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
   use coordinatesystems, only: spherical_polar_t, cartesian3d_t
   use physicalconstants, only: dd_pi
!
   implicit none
!
!
   private
!
   integer, parameter, public :: nfaces = 6 ! number of faces on the cube
   integer, parameter, public :: ninnerelemedge = 8 ! number of edges for an interior element
   integer, parameter, public :: ncornerelemedge = 4 ! number of corner elements
   real(real_kind), public, parameter :: cube_xstart =-0.25d0*dd_pi
   real(real_kind), public, parameter :: cube_xend = 0.25d0*dd_pi
   real(real_kind), public, parameter :: cube_ystart =-0.25d0*dd_pi
   real(real_kind), public, parameter :: cube_yend = 0.25d0*dd_pi
!
   type, public :: face_t
     sequence
     type(spherical_polar_t) :: sphere0 ! tangent point of face on sphere
     type(spherical_polar_t) :: sw ! sw corner of face on sphere
     type(spherical_polar_t) :: se ! se corner of face on sphere
     type(spherical_polar_t) :: ne ! ne corner of face on sphere
     type(spherical_polar_t) :: nw ! nw corner of face on sphere
     type(cartesian3d_t) :: p0
     type(cartesian3d_t) :: x0
     type(cartesian3d_t) :: y0
     integer :: number
     integer :: padding ! padd the struct
   end type face_t
!
   type, public :: cube_face_coord_t
     sequence
     real(real_kind) :: x ! x coordinate
     real(real_kind) :: y ! y coordinate
     type(face_t), pointer :: face ! face
   end type cube_face_coord_t
! ==========================================
! Public Interfaces
! ==========================================
   public :: cubetopology
! Rotate the North Pole:  used for JW baroclinic test case
! Settings this only changes Coriolis.
! User must also rotate initial condition
   real(real_kind), public :: rotate_grid = 0
! ===============================
! Public methods for cube
! ===============================
   public :: cube_init_atomic
   public :: convert_gbl_index
   public :: cube_assemble
   public :: vmap, dmap
   public :: covariant_rot
   public :: contravariant_rot
   public :: set_corner_coordinates
   public :: assign_node_numbers_to_elem
   public :: cubeedgecount
   public :: cubeelemcount
   public :: cubesetupedgeindex
   public :: rotation_init_atomic
! ===============================
! Private methods
! ===============================
   private :: coordinates_atomic
   private :: metric_atomic
   private :: coreolis_init, coreolis_init_atomic
   private :: getlatticespacing
!
!
   contains
!
! =======================================
!  cube_init_atomic:
!
! Initialize element descriptors for
! cube sphere case for each element ...
! =======================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine cube_init_atomic(elem, gll_points, alpha_in)
!-------------------------------------------------------------------------------
   use element, only: element_t
   use dimensions, only: np
!
   type(element_t), intent(inout) :: elem
   real(real_kind), optional :: alpha_in
   real(real_kind) :: alpha = 1
   real(longdouble_kind) :: gll_points(np)
!
   if (present(alpha_in)) alpha = alpha_in
   elem%facenum = elem%vertex%face_number
   call coordinates_atomic(elem,gll_points)
   call metric_atomic(elem,gll_points,alpha)
   call coreolis_init_atomic(elem)
   elem%desc%use_rotation = 0
   call solver_weights_atomic(elem)
!
   end subroutine cube_init_atomic
! =======================================
! coordinates_atomic:
!
! Initialize element coordinates for
! cube-sphere case ... (atomic)
!
! =======================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine coordinates_atomic(elem, gll_points)
!-------------------------------------------------------------------------------
   use element, only: element_t, element_var_coordinates
   use coordinatesystems, only: cartesian2d_t, ref2sphere, cubedsphere2cart, spherical_to_cart, sphere_tri_area
   use dimensions, only: np
   use kiapsparallel, only: abortpar
!
   type(element_t) :: elem
   real(longdouble_kind) :: gll_points(np)
   real(real_kind) :: area1, area2
   type(cartesian3d_t) :: quad(4)
   integer face_no, i, j
! =========================================
! compute coordinates of each GLL point
! =========================================
!
   face_no = elem%vertex%face_number
   do i = 1,np
     do j = 1, np
       elem%spherep(i,j) = ref2sphere(gll_points(i),gll_points(j),elem%corners,face_no)
     enddo
   enddo
! compute the corners in Cartesian coordinates
   do i = 1,4
     elem%corners3d(i) = cubedsphere2cart(elem%corners(i),face_no)
   enddo
! also compute the [-pi/2,pi/2] cubed sphere coordinates:
   elem%cartp = element_var_coordinates(elem%corners,gll_points)
! Matrix describing vector conversion to cartesian
! Zonal direction
   elem%vec_sphere2cart(:,:,1,1) =-sin(elem%spherep(:,:)%lon)
   elem%vec_sphere2cart(:,:,2,1) = cos(elem%spherep(:,:)%lon)
   elem%vec_sphere2cart(:,:,3,1) = 0.0_real_kind
! Meridional direction
   elem%vec_sphere2cart(:,:,1,2) =-sin(elem%spherep(:,:)%lat)*cos(elem%spherep(:,:)%lon)
   elem%vec_sphere2cart(:,:,2,2) =-sin(elem%spherep(:,:)%lat)*sin(elem%spherep(:,:)%lon)
   elem%vec_sphere2cart(:,:,3,2) = cos(elem%spherep(:,:)%lat)
!
   end subroutine coordinates_atomic
! elem_jacobians:
!
! Calculate Jacobian associated with mapping
! from arbitrary quadrilateral to [-1,1]^2
! along with its inverse and determinant
! ==========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine elem_jacobians(coords, unif2quadmap, nx)
!-------------------------------------------------------------------------------
   use coordinatesystems, only: cartesian2d_t
!    use Quadrature, only : quadrature_t
!
   integer, intent(in   ) :: nx
   type(cartesian2d_t), dimension(nx, nx), intent(in   ) :: coords
! unif2quadmap is the bilinear map from [-1,1]^2 -> arbitrary quadrilateral
   real(real_kind), dimension(4, 2), intent(  out) :: unif2quadmap
!    type (quadrature_t), intent(in) :: gauss
!    real (kind=real_kind), dimension(2,2,nx,nx), intent(out) :: J, Jinv
!    real (kind=real_kind), dimension(nx,nx), intent(out) :: Jdet
   integer :: ii, jj
!
   unif2quadmap(1,1) =(coords(1,1)%x+coords(nx,1)%x+coords(nx,nx)%x+coords(1,nx)%x)/4.0d0
   unif2quadmap(1,2) =(coords(1,1)%y+coords(nx,1)%y+coords(nx,nx)%y+coords(1,nx)%y)/4.0d0
   unif2quadmap(2,1) =(-coords(1,1)%x+coords(nx,1)%x+coords(nx,nx)%x-coords(1,nx)%x)/4.0d0
   unif2quadmap(2,2) =(-coords(1,1)%y+coords(nx,1)%y+coords(nx,nx)%y-coords(1,nx)%y)/4.0d0
   unif2quadmap(3,1) =(-coords(1,1)%x-coords(nx,1)%x+coords(nx,nx)%x+coords(1,nx)%x)/4.0d0
   unif2quadmap(3,2) =(-coords(1,1)%y-coords(nx,1)%y+coords(nx,nx)%y+coords(1,nx)%y)/4.0d0
   unif2quadmap(4,1) =(coords(1,1)%x-coords(nx,1)%x+coords(nx,nx)%x-coords(1,nx)%x)/4.0d0
   unif2quadmap(4,2) =(coords(1,1)%y-coords(nx,1)%y+coords(nx,nx)%y-coords(1,nx)%y)/4.0d0
#if 0
   do ii = 1,nx
     do jj = 1, nx
       j(1,1,ii,jj) = unif2quadmap(2,1)+unif2quadmap(4,1)*gauss%points(jj)
       j(1,2,ii,jj) = unif2quadmap(3,1)+unif2quadmap(4,1)*gauss%points(ii)
       j(2,1,ii,jj) = unif2quadmap(2,2)+unif2quadmap(4,2)*gauss%points(jj)
       j(2,2,ii,jj) = unif2quadmap(3,2)+unif2quadmap(4,2)*gauss%points(ii)
       jdet(ii,jj) = j(1,1,ii,jj)*j(2,2,ii,jj)-j(1,2,ii,jj)*j(2,1,ii,jj)
       jinv(1,1,ii,jj) = j(2,2,ii,jj)/jdet(ii,jj)
       jinv(1,2,ii,jj) =-j(1,2,ii,jj)/jdet(ii,jj)
       jinv(2,1,ii,jj) =-j(2,1,ii,jj)/jdet(ii,jj)
       jinv(2,2,ii,jj) = j(1,1,ii,jj)/jdet(ii,jj)
     enddo
   enddo
#endif
!
   end subroutine elem_jacobians
! =========================================
! metric_atomic:
!
! Initialize cube-sphere metric terms:
! equal angular elements (atomic)
! initialize:
!         metdet, rmetdet  (analytic)    = detD, 1/detD
!         met                (analytic)    D'D or DD' ?
!         metdet             (analytic)    = detD
!         metinv             (analytic)    Dinv'Dinv  or Dinv Dinv' ?
!         D     (from subroutine vmap)
!         Dinv  (computed directly from D)
!
! so if we want to tweak the mapping by a factor alpha (so he weights add up to 4pi, for example)
! we take:
!    NEW       OLD
!       D = sqrt(alpha) D  and then rederive all quantities.
!    detD = alpha detD
!
! where alpha = 4pi/SEMarea, SEMarea = global sum elem(ie)%mv(i,j)*elem(ie)%metdet(i,j)
!
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine metric_atomic(elem, gll_points, alpha)
!-------------------------------------------------------------------------------
   use element, only: element_t
   use dimensions, only: np
   use physicalconstants, only: rrearth=>kim_rrearth
   use kiapsparallel, only: abortpar
!
   type(element_t) :: elem
   real(real_kind) :: alpha
   real(longdouble_kind) :: gll_points(np)
! Local variables
   integer ii, face_no
   integer i, j
   integer iptr
   real(real_kind) :: r ! distance from origin for point on cube tangent to unit sphere
   real(real_kind) :: const, norm
   real(real_kind) :: detd ! determinant of vector field mapping matrix.
   real(real_kind) :: x1 ! 1st cube face coordinate
   real(real_kind) :: x2 ! 2nd cube face coordinate
   real(real_kind) :: tmpd(2, 2)
   real(real_kind) :: l1, l2 ! eigen values of met
!
   face_no = elem%vertex%face_number
! ==============================================
! Initialize differential mapping operator
! to and from vector fields on the sphere to
! contravariant vector fields on the cube
! i.e. dM/dx^i in Sadourney (1972) and it's
! inverse
! ==============================================
! MNL: Calculate Jacobians.  these must be computed before Dmap is used below
   call elem_jacobians(elem%cartp,elem%u2qmap,np)
   elem%max_eig = 0.0d0
   elem%min_eig = 1d99
   do j = 1,np
     do i = 1, np
       x1 = gll_points(i)
       x2 = gll_points(j)
       call dmap(elem%d(:,:,i,j),elem,x1,x2)
! Numerical metric tensor based on analytic D: met = D^T times D
! (D maps between sphere and reference element)
       elem%met(1,1,i,j) = elem%d(1,1,i,j)*elem%d(1,1,i,j)+elem%d(2,1,i,j)*elem%d(2,1,i,j)
       elem%met(1,2,i,j) = elem%d(1,1,i,j)*elem%d(1,2,i,j)+elem%d(2,1,i,j)*elem%d(2,2,i,j)
       elem%met(2,1,i,j) = elem%d(1,1,i,j)*elem%d(1,2,i,j)+elem%d(2,1,i,j)*elem%d(2,2,i,j)
       elem%met(2,2,i,j) = elem%d(1,2,i,j)*elem%d(1,2,i,j)+elem%d(2,2,i,j)*elem%d(2,2,i,j)
! compute D^-1...
! compute determinant of D mapping matrix... if not zero compute inverse
       detd = elem%d(1,1,i,j)*elem%d(2,2,i,j)-elem%d(1,2,i,j)*elem%d(2,1,i,j)
       elem%dinv(1,1,i,j) = elem%d(2,2,i,j)/detd
       elem%dinv(1,2,i,j) =-elem%d(1,2,i,j)/detd
       elem%dinv(2,1,i,j) =-elem%d(2,1,i,j)/detd
       elem%dinv(2,2,i,j) = elem%d(1,1,i,j)/detd
! L2 norm = sqrt max eigenvalue of metinv
!         = 1/sqrt(min eigenvalue of met)
! l1 and l2 are eigenvalues of met
! (should both be positive,l1 > l2)
       l1 =(elem%met(1,1,i,j)+elem%met(2,2,i,j)+sqrt(4.0d0*elem%met(1,2,i,j)*elem%met(2,1,i,j)+(elem%met(1,1,i,j)-elem%met(2,2,i,j))**2))/2.0d0
       l2 =(elem%met(1,1,i,j)+elem%met(2,2,i,j)-sqrt(4.0d0*elem%met(1,2,i,j)*elem%met(2,1,i,j)+(elem%met(1,1,i,j)-elem%met(2,2,i,j))**2))/2.0d0
! Max L2 norm of Dinv is sqrt of max eigenvalue of metinv
! max eigenvalue of metinv is 1/min eigenvalue of met
       norm = 1.0d0/sqrt(min(abs(l1),abs(l2)))
       elem%max_eig = max(norm,elem%max_eig)
! Min L2 norm of Dinv is sqrt of min eigenvalue of metinv
! min eigenvalue of metinv is 1/max eigenvalue of met
       norm = 1.0d0/sqrt(max(abs(l1),abs(l2)))
       elem%min_eig = min(norm,elem%min_eig)
! Need inverse of met if not calculated analytically
       elem%metdet(i,j) = abs(detd)
       elem%rmetdet(i,j) = 1.0d0/abs(detd)
       elem%metinv(1,1,i,j) = elem%met(2,2,i,j)/(detd*detd)
       elem%metinv(1,2,i,j) =-elem%met(1,2,i,j)/(detd*detd)
       elem%metinv(2,1,i,j) =-elem%met(2,1,i,j)/(detd*detd)
       elem%metinv(2,2,i,j) = elem%met(1,1,i,j)/(detd*detd)
     enddo
   enddo
   elem%dx_short = 1.0d0/(elem%max_eig*0.5d0*dble(np-1)*rrearth*1000.0d0)
   elem%dx_long = 1.0d0/(elem%min_eig*0.5d0*dble(np-1)*rrearth*1000.0d0)
! ===============================================
!
! Initialize equal angular metric tensor on each
! on velocity grid for unit sphere.
!
! Initialize gdet = SQRT(ABS(DET(gij)))
!
! These quantities are the same on every face
! of the cube.
!
! =================================================
! mt: better might be to compute all these quantities directly from D
! for consistency?
!
! MNL: done
   elem%d = elem%d*sqrt(alpha)
   elem%dinv = elem%dinv/sqrt(alpha)
   elem%metdet = elem%metdet*alpha
   elem%rmetdet = elem%rmetdet/alpha
   elem%met = elem%met*alpha
   elem%metinv = elem%metinv/alpha
!
   end subroutine metric_atomic
! =======================================
! solver_weights:
!
! For nonstaggered GaussLobatto elements,
! compute weights for redundant points in
! cg solver.
!
! =======================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine solver_weights_atomic(elem)
!-------------------------------------------------------------------------------
   use element, only: element_t
   use dimensions, only: np
   use kiapsparallel, only: abortpar
!
   type(element_t) :: elem
   real(real_kind) :: x
! Local variables
   integer :: i, j
! =========================================
! compute cube face coordinates of element
! =========================================
!
   do i = 1, np
     do j = 1, np
       if (i==1) then
         if (j==1) then
           x = 1.0_real_kind/elem%node_multiplicity(1)
         elseif (j==np) then
           x = 1.0_real_kind/elem%node_multiplicity(4)
         else
           x = 0.5_real_kind
         endif
         elseif (i==np) then
           if (j==1) then
             x = 1.0_real_kind/elem%node_multiplicity(2)
           elseif (j==np) then
             x = 1.0_real_kind/elem%node_multiplicity(3)
           else
             x = 0.5_real_kind
           endif
           elseif (j==1.or.j==np) then
             x = 0.5_real_kind
           else
             x = 1.0_real_kind
           endif
             elem%solver_wts(i,j) = x
           enddo
   enddo
!
   end subroutine solver_weights_atomic
#if 1
! ========================================
! covariant_rot:
!
! 2 x 2 matrix multiply:  Db^T * Da^-T
! for edge rotations: maps face a to face b
!
! ========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function covariant_rot(da, db) result(r)
!
   real(real_kind) :: da(2, 2)
   real(real_kind) :: db(2, 2)
   real(real_kind) :: r(2, 2)
   real(real_kind) :: detda
!
   detda = da(2,2)*da(1,1)-da(1,2)*da(2,1)
   r(1,1) =(da(2,2)*db(1,1)-da(1,2)*db(2,1))/detda
   r(1,2) =(da(1,1)*db(2,1)-da(2,1)*db(1,1))/detda
   r(2,1) =(da(2,2)*db(1,2)-da(1,2)*db(2,2))/detda
   r(2,2) =(da(1,1)*db(2,2)-da(2,1)*db(1,2))/detda
!
   end function covariant_rot
#else
! ========================================
! covariant_rot:
!
! 2 x 2 matrix multiply:  Db * Da^-1
! for edge rotations: maps face a to face b
!
! ========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function covariant_rot(da, db) result(r)
!
   real(real_kind) :: da(2, 2)
   real(real_kind) :: db(2, 2)
   real(real_kind) :: r(2, 2)
   real(real_kind) :: detda
!
   detda = da(2,2)*da(1,1)-da(1,2)*da(2,1)
   r(1,1) =(da(2,2)*db(1,1)-da(2,1)*db(1,2))/detda
   r(1,2) =(da(1,1)*db(1,2)-da(1,2)*db(1,1))/detda
   r(2,1) =(da(2,2)*db(2,1)-da(2,1)*db(2,2))/detda
   r(2,2) =(da(1,1)*db(2,2)-da(1,2)*db(2,1))/detda
!
   end function covariant_rot
#endif
! ========================================
! contravariant_rot:
!
! 2 x 2 matrix multiply:  Db^-1 * Da
! that maps a contravariant vector field
! from an edge of cube face a to a contiguous
! edge of cube face b.
!
! ========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function contravariant_rot(da, db) result(r)
!
   real(real_kind) :: da(2, 2)
   real(real_kind) :: db(2, 2)
   real(real_kind) :: r(2, 2)
   real(real_kind) :: detdb
!
   detdb = db(2,2)*db(1,1)-db(1,2)*db(2,1)
   r(1,1) =(da(1,1)*db(2,2)-da(2,1)*db(1,2))/detdb
   r(1,2) =(da(1,2)*db(2,2)-da(2,2)*db(1,2))/detdb
   r(2,1) =(da(2,1)*db(1,1)-da(1,1)*db(2,1))/detdb
   r(2,2) =(da(2,2)*db(1,1)-da(1,2)*db(2,1))/detdb
!
   end function contravariant_rot
! ========================================================
! Dmap:
!
! Initialize mapping that tranforms contravariant
! vector fields on the reference element onto vector fields on
! the sphere.
! For Gnomonic, followed by bilinear, this code uses the old vmap()
! for unstructured grids, this code uses the parametric map that
! maps quads on the sphere directly to the reference element
! ========================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine dmap(d, elem, a, b)
!-------------------------------------------------------------------------------
   use coordinatesystems, only: dist_threshold
   use element, only: element_t
!
   type(element_t) :: elem
   real(real_kind), intent(  out) :: d(2, 2)
   real(real_kind), intent(in   ) :: a, b
! local
   real(real_kind) :: tmpd(2, 2), jp(2, 2), x1, x2, pi, pj, qi, qj
! input (a,b) shold be a point in the reference element [-1,1]
! compute Jp(a,b)
!
   jp(1,1) = elem%u2qmap(2,1)+elem%u2qmap(4,1)*b
   jp(1,2) = elem%u2qmap(3,1)+elem%u2qmap(4,1)*a
   jp(2,1) = elem%u2qmap(2,2)+elem%u2qmap(4,2)*b
   jp(2,2) = elem%u2qmap(3,2)+elem%u2qmap(4,2)*a
! map (a,b) to the [-pi/2,pi/2] equi angular cube face:  x1,x2
! a = gp%points(i)
! b = gp%points(j)
   pi =(1-a)/2
   pj =(1-b)/2
   qi =(1+a)/2
   qj =(1+b)/2
   x1 = pi*pj*elem%corners(1)%x+qi*pj*elem%corners(2)%x+qi*qj*elem%corners(3)%x+pi*qj*elem%corners(4)%x
   x2 = pi*pj*elem%corners(1)%y+qi*pj*elem%corners(2)%y+qi*qj*elem%corners(3)%y+pi*qj*elem%corners(4)%y
   call vmap(tmpd,x1,x2,elem%vertex%face_number)
! Include map from element -> ref element in D
   d(1,1) = tmpd(1,1)*jp(1,1)+tmpd(1,2)*jp(2,1)
   d(1,2) = tmpd(1,1)*jp(1,2)+tmpd(1,2)*jp(2,2)
   d(2,1) = tmpd(2,1)*jp(1,1)+tmpd(2,2)*jp(2,1)
   d(2,2) = tmpd(2,1)*jp(1,2)+tmpd(2,2)*jp(2,2)
!
   end subroutine dmap
! ========================================================
! vmap:
!
! Initialize mapping that tranforms contravariant
! vector fields on the cube onto vector fields on
! the sphere. This follows Taylor's D matrix
!
!       | cos(theta)dlambda/dx1  cos(theta)dlambda/dx2 |
!   D = |                                              |
!       |     dtheta/dx1              dtheta/dx2       |
!
! ========================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine vmap(d, x1, x2, face_no)
!-------------------------------------------------------------------------------
   use coordinatesystems, only: dist_threshold
!
   real(real_kind), intent(inout) :: d(2, 2)
   real(real_kind), intent(in   ) :: x1
   real(real_kind), intent(in   ) :: x2
   integer, intent(in   ) :: face_no
! Local variables
   real(real_kind) :: poledist ! sqrt(tan(x1)**2+tan(x2)**2)
   real(real_kind) :: r ! distance from cube point to center of sphere
   real(real_kind) :: d11
   real(real_kind) :: d12
   real(real_kind) :: d21
   real(real_kind) :: d22
!
   r = sqrt(1.0d0+tan(x1)**2+tan(x2)**2)
   if (face_no>= 1.and.face_no<=4) then
     d11 = 1.0d0/(r*cos(x1))
     d12 = 0.0d0
     d21 =-tan(x1)*tan(x2)/(cos(x1)*r*r)
     d22 = 1.0d0/(r*r*cos(x1)*cos(x2)*cos(x2))
     d(1,1) = d11
     d(1,2) = d12
     d(2,1) = d21
     d(2,2) = d22
   elseif (face_no==6) then
     poledist = sqrt(tan(x1)**2+tan(x2)**2)
     if (poledist<=dist_threshold) then
! we set the D transform to the identity matrix
! which works ONLY for swtc1, phi starting at
! 3*PI/2... assumes lon at pole == 0
       d(1,1) = 1.0d0
       d(1,2) = 0.0d0
       d(2,1) = 0.0d0
       d(2,2) = 1.0d0
     else
       d11 =-tan(x2)/(poledist*cos(x1)*cos(x1)*r)
       d12 = tan(x1)/(poledist*cos(x2)*cos(x2)*r)
       d21 =-tan(x1)/(poledist*cos(x1)*cos(x1)*r*r)
       d22 =-tan(x2)/(poledist*cos(x2)*cos(x2)*r*r)
       d(1,1) = d11
       d(1,2) = d12
       d(2,1) = d21
       d(2,2) = d22
     endif
   elseif (face_no==5) then
     poledist = sqrt(tan(x1)**2+tan(x2)**2)
     if (poledist<=dist_threshold) then
! we set the D transform to the identity matrix
! which works ONLY for swtc1, phi starting at
! 3*PI/2... assumes lon at pole == 0, i.e. very specific
       d(1,1) = 1.0d0
       d(1,2) = 0.0d0
       d(2,1) = 0.0d0
       d(2,2) = 1.0d0
     else
       d11 = tan(x2)/(poledist*cos(x1)*cos(x1)*r)
       d12 =-tan(x1)/(poledist*cos(x2)*cos(x2)*r)
       d21 = tan(x1)/(poledist*cos(x1)*cos(x1)*r*r)
       d22 = tan(x2)/(poledist*cos(x2)*cos(x2)*r*r)
       d(1,1) = d11
       d(1,2) = d12
       d(2,1) = d21
       d(2,2) = d22
     endif
   endif
!
   end subroutine vmap
! ========================================
! coreolis_init:
!
! Initialize coreolis term ...
!
! ========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine coreolis_init(elem)
!-------------------------------------------------------------------------------
   use element, only: element_t
   use dimensions, only: np
   use physicalconstants, only: omega=>kim_omega
!
   type(element_t) :: elem(:)
! Local variables
   integer :: i, j
   integer :: ii
!
   do ii = 1, size(elem)
     call coreolis_init_atomic(elem(ii))
   enddo
!
   end subroutine coreolis_init
! ========================================
! coreolis_init_atomic:
!
! Initialize coreolis term ...
!
! ========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine coreolis_init_atomic(elem)
!-------------------------------------------------------------------------------
   use element, only: element_t
   use dimensions, only: np
   use physicalconstants, only: omega=>kim_omega
!
   type(element_t) :: elem
! Local variables
   integer :: i, j
   real(real_kind) :: lat, lon, rangle
!
   rangle = rotate_grid*dd_pi/180
   do j = 1,np
     do i = 1, np
       if (rotate_grid/= 0) then
         lat = elem%spherep(i,j)%lat
         lon = elem%spherep(i,j)%lon
         elem%fcor(i,j) = 2*omega*(-cos(lon)*cos(lat)*sin(rangle)+sin(lat)*cos(rangle))
       else
         elem%fcor(i,j) = 2.0d0*omega*sin(elem%spherep(i,j)%lat)
#if defined(CORE_SW)
!SJ.CHOI 14.03.17
         elem%ecor(i,j) = 2.0d0*omega*cos(elem%spherep(i,j)%lat)
#endif
       endif
       enddo
   enddo
!
   end subroutine coreolis_init_atomic
! =========================================
! rotation_init_atomic:
!
! Initialize cube rotation terms resulting
! from changing cube face coordinate systems
!
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine rotation_init_atomic(elem, rot_type)
!-------------------------------------------------------------------------------
   use element, only: element_t
   use dimensions, only: np
   use control, only: north, south, east, west, neast, seast, swest, nwest
   use kiapsparallel, only: abortpar
!
   type(element_t) :: elem
   character(len=*) rot_type
! =======================================
! Local variables
! =======================================
   integer :: myface_no ! current element face number
   integer :: nbrface_no ! neighbor element face number
   integer :: inbr
   integer :: nrot, irot
   integer :: ii, i, j, k
   integer :: ir, jr
   real(real_kind) :: dloc(2, 2, np)
   real(real_kind) :: drem(2, 2, np)
   real(real_kind) :: x1, x2
!
   myface_no = elem%vertex%face_number
   nrot = 0
   do inbr = 1,8
     if (associated(elem%vertex%nbrs(inbr)%n)) then
       do k = 1, size(elem%vertex%nbrs(inbr)%n)
         nbrface_no = elem%vertex%nbrs(inbr)%f(k)
         if (myface_no/= nbrface_no) nrot = nrot+1
       enddo
       endif
   enddo
   if (associated(elem%desc%rot)) then
     if (size(elem%desc%rot)>0) then
!         deallocate(elem%desc%rot)
       nullify(elem%desc%rot)
     endif
   endif
! =====================================================
! If there are neighbors on other cube faces,allocate
! an array of rotation matrix structs.
! =====================================================
   if (nrot>0) then
     allocate(elem%desc%rot(nrot))
     elem%desc%use_rotation = 1
     irot = 0
     do inbr = 1,8
       if (associated(elem%vertex%nbrs(inbr)%n)) then
         do k = 1, size(elem%vertex%nbrs(inbr)%n)
           nbrface_no = elem%vertex%nbrs(inbr)%f(k)
! The cube edge (myface_no,nbrface_no) and inbr defines
! a unique rotation given by (D^-1) on myface_no x (D) on nbrface_no
           if (myface_no/= nbrface_no.and.elem%vertex%nbrs(inbr)%n(k)/=-1) then
             irot = irot+1
             if (inbr<=4) then
               allocate(elem%desc%rot(irot)%r(2,2,np)) ! edge
             else
               allocate(elem%desc%rot(irot)%r(2,2,1)) ! corner
             endif
! must compute Dloc on my face,Drem on neighbor face,
! for each point on edge or corner.
! ====================================
! Equatorial belt east/west neighbors
! ====================================
             if (nbrface_no<=4.and.myface_no<=4) then
               if (inbr==west) then
                 do j = 1, np
                   x1 = elem%cartp(1,j)%x
                   x2 = elem%cartp(1,j)%y
                   call vmap(dloc(1,1,j),x1,x2,myface_no)
                   call vmap(drem(1,1,j),-x1,x2,nbrface_no)
                 enddo
                 elseif (inbr==east) then
                   do j = 1, np
                     x1 = elem%cartp(np,j)%x
                     x2 = elem%cartp(np,j)%y
                     call vmap(dloc(1,1,j),x1,x2,myface_no)
                     call vmap(drem(1,1,j),-x1,x2,nbrface_no)
                   enddo
                   elseif (inbr==swest) then
                     x1 = elem%cartp(1,1)%x
                     x2 = elem%cartp(1,1)%y
                     call vmap(dloc(1,1,1),x1,x2,myface_no)
                     call vmap(drem(1,1,1),-x1,x2,nbrface_no)
                   elseif (inbr==nwest) then
                     x1 = elem%cartp(1,np)%x
                     x2 = elem%cartp(1,np)%y
                     call vmap(dloc(1,1,1),x1,x2,myface_no)
                     call vmap(drem(1,1,1),-x1,x2,nbrface_no)
                   elseif (inbr==seast) then
                     x1 = elem%cartp(np,1)%x
                     x2 = elem%cartp(np,1)%y
                     call vmap(dloc(1,1,1),x1,x2,myface_no)
                     call vmap(drem(1,1,1),-x1,x2,nbrface_no)
                   elseif (inbr==neast) then
                     x1 = elem%cartp(np,np)%x
                     x2 = elem%cartp(np,np)%y
                     call vmap(dloc(1,1,1),x1,x2,myface_no)
                     call vmap(drem(1,1,1),-x1,x2,nbrface_no)
                   endif
             endif
! Northern Neighbors of Equatorial Belt
             if (myface_no<=4.and.nbrface_no==6) then
               if (inbr==north) then
                 do i = 1, np
                   ir = np+1-i
                   x1 = elem%cartp(i,np)%x
                   x2 = elem%cartp(i,np)%y
                   if (myface_no==1) then
                     call vmap(dloc(1,1,i),x1,x2,myface_no)
                     call vmap(drem(1,1,i),x1,-x2,nbrface_no)
                   endif
                   if (myface_no==2) then
                     call vmap(dloc(1,1,i),x1,x2,myface_no)
                     call vmap(drem(1,1,i),x2,x1,nbrface_no)
                   endif
                   if (myface_no==3) then
                     call vmap(dloc(1,1,ir),x1,x2,myface_no)
                     call vmap(drem(1,1,ir),-x1,x2,nbrface_no)
                   endif
                   if (myface_no==4) then
                     call vmap(dloc(1,1,ir),x1,x2,myface_no)
                     call vmap(drem(1,1,ir),-x2,-x1,nbrface_no)
                   endif
                 enddo
                 elseif (inbr==nwest) then
                   x1 = elem%cartp(1,np)%x
                   x2 = elem%cartp(1,np)%y
                   call vmap(dloc(1,1,1),x1,x2,myface_no)
                   if (myface_no==1) call vmap(drem(1,1,1),x1,-x2,nbrface_no)
                   if (myface_no==2) call vmap(drem(1,1,1),x2,x1,nbrface_no)
                   if (myface_no==3) call vmap(drem(1,1,1),-x1,x2,nbrface_no)
                   if (myface_no==4) call vmap(drem(1,1,1),-x2,-x1,nbrface_no)
                 elseif (inbr==neast) then
                   x1 = elem%cartp(np,np)%x
                   x2 = elem%cartp(np,np)%y
                   call vmap(dloc(1,1,1),x1,x2,myface_no)
                   if (myface_no==1) call vmap(drem(1,1,1),x1,-x2,nbrface_no)
                   if (myface_no==2) call vmap(drem(1,1,1),x2,x1,nbrface_no)
                   if (myface_no==3) call vmap(drem(1,1,1),-x1,x2,nbrface_no)
                   if (myface_no==4) call vmap(drem(1,1,1),-x2,-x1,nbrface_no)
                 endif
             endif
! Southern Neighbors of Equatorial Belt
             if (myface_no<=4.and.nbrface_no==5) then
               if (inbr==south) then
                 do i = 1, np
                   ir = np+1-i
                   x1 = elem%cartp(i,1)%x
                   x2 = elem%cartp(i,1)%y
                   if (myface_no==1) then
                     call vmap(dloc(1,1,i),x1,x2,myface_no)
                     call vmap(drem(1,1,i),x1,-x2,nbrface_no)
                   endif
                   if (myface_no==2) then
                     call vmap(dloc(1,1,ir),x1,x2,myface_no)
                     call vmap(drem(1,1,ir),-x2,-x1,nbrface_no)
                   endif
                   if (myface_no==3) then
                     call vmap(dloc(1,1,ir),x1,x2,myface_no)
                     call vmap(drem(1,1,ir),-x1,x2,nbrface_no)
                   endif
                   if (myface_no==4) then
                     call vmap(dloc(1,1,i),x1,x2,myface_no)
                     call vmap(drem(1,1,i),x2,x1,nbrface_no)
                   endif
                 enddo
                 elseif (inbr==swest) then
                   x1 = elem%cartp(1,1)%x
                   x2 = elem%cartp(1,1)%y
                   call vmap(dloc(1,1,1),x1,x2,myface_no)
                   if (myface_no==1) call vmap(drem(1,1,1),x1,-x2,nbrface_no)
                   if (myface_no==2) call vmap(drem(1,1,1),-x2,-x1,nbrface_no)
                   if (myface_no==3) call vmap(drem(1,1,1),-x1,x2,nbrface_no)
                   if (myface_no==4) call vmap(drem(1,1,1),x2,x1,nbrface_no)
                 elseif (inbr==seast) then
                   x1 = elem%cartp(np,1)%x
                   x2 = elem%cartp(np,1)%y
                   call vmap(dloc(1,1,1),x1,x2,myface_no)
                   if (myface_no==1) call vmap(drem(1,1,1),x1,-x2,nbrface_no)
                   if (myface_no==2) call vmap(drem(1,1,1),-x2,-x1,nbrface_no)
                   if (myface_no==3) call vmap(drem(1,1,1),-x1,x2,nbrface_no)
                   if (myface_no==4) call vmap(drem(1,1,1),x2,x1,nbrface_no)
                 endif
             endif
! Neighbors of Northern Capping Face Number 6
             if (myface_no==6) then
               if (nbrface_no==1) then
                 if (inbr==south) then
                   do i = 1, np
                     x1 = elem%cartp(i,1)%x
                     x2 = elem%cartp(i,1)%y
                     call vmap(dloc(1,1,i),x1,x2,myface_no)
                     call vmap(drem(1,1,i),x1,-x2,nbrface_no)
                   enddo
                   elseif (inbr==swest) then
                     x1 = elem%cartp(1,1)%x
                     x2 = elem%cartp(1,1)%y
                     call vmap(dloc(1,1,1),x1,x2,myface_no)
                     call vmap(drem(1,1,1),x1,-x2,nbrface_no)
                   elseif (inbr==seast) then
                     x1 = elem%cartp(np,1)%x
                     x2 = elem%cartp(np,1)%y
                     call vmap(dloc(1,1,1),x1,x2,myface_no)
                     call vmap(drem(1,1,1),x1,-x2,nbrface_no)
                   endif
                   elseif (nbrface_no==2) then
                     if (inbr==east) then
                       do j = 1, np
                         x1 = elem%cartp(np,j)%x
                         x2 = elem%cartp(np,j)%y
                         call vmap(dloc(1,1,j),x1,x2,myface_no)
                         call vmap(drem(1,1,j),x2,x1,nbrface_no)
                       enddo
                       elseif (inbr==seast) then
                         x1 = elem%cartp(np,1)%x
                         x2 = elem%cartp(np,1)%y
                         call vmap(dloc(1,1,1),x1,x2,myface_no)
                         call vmap(drem(1,1,1),x2,x1,nbrface_no)
                       elseif (inbr==neast) then
                         x1 = elem%cartp(np,np)%x
                         x2 = elem%cartp(np,np)%y
                         call vmap(dloc(1,1,1),x1,x2,myface_no)
                         call vmap(drem(1,1,1),x2,x1,nbrface_no)
                       endif
                       elseif (nbrface_no==3) then
                         if (inbr==north) then
                           do i = 1, np
                             ir = np+1-i
                             x1 = elem%cartp(i,np)%x
                             x2 = elem%cartp(i,np)%y
                             call vmap(dloc(1,1,ir),x1,x2,myface_no)
                             call vmap(drem(1,1,ir),-x1,x2,nbrface_no)
                           enddo
                           elseif (inbr==nwest) then
                             x1 = elem%cartp(1,np)%x
                             x2 = elem%cartp(1,np)%y
                             call vmap(dloc(1,1,1),x1,x2,myface_no)
                             call vmap(drem(1,1,1),-x1,x2,nbrface_no)
                           elseif (inbr==neast) then
                             x1 = elem%cartp(np,np)%x
                             x2 = elem%cartp(np,np)%y
                             call vmap(dloc(1,1,1),x1,x2,myface_no)
                             call vmap(drem(1,1,1),-x1,x2,nbrface_no)
                           endif
                           elseif (nbrface_no==4) then
                             if (inbr==west) then
                               do j = 1, np
                                 jr = np+1-j
                                 x1 = elem%cartp(1,j)%x
                                 x2 = elem%cartp(1,j)%y
                                 call vmap(dloc(1,1,jr),x1,x2,myface_no)
                                 call vmap(drem(1,1,jr),-x2,-x1,nbrface_no)
                               enddo
                               elseif (inbr==swest) then
                                 x1 = elem%cartp(1,1)%x
                                 x2 = elem%cartp(1,1)%y
                                 call vmap(dloc(1,1,1),x1,x2,myface_no)
                                 call vmap(drem(1,1,1),-x2,-x1,nbrface_no)
                               elseif (inbr==nwest) then
                                 x1 = elem%cartp(1,np)%x
                                 x2 = elem%cartp(1,np)%y
                                 call vmap(dloc(1,1,1),x1,x2,myface_no)
                                 call vmap(drem(1,1,1),-x2,-x1,nbrface_no)
                               endif
                               endif
             endif
! Neighbors of South Capping Face Number 5
             if (myface_no==5) then
               if (nbrface_no==1) then
                 if (inbr==north) then
                   do i = 1, np
                     x1 = elem%cartp(i,np)%x
                     x2 = elem%cartp(i,np)%y
                     call vmap(dloc(1,1,i),x1,x2,myface_no)
                     call vmap(drem(1,1,i),x1,-x2,nbrface_no)
                   enddo
                   elseif (inbr==nwest) then
                     x1 = elem%cartp(1,np)%x
                     x2 = elem%cartp(1,np)%y
                     call vmap(dloc(:,:,1),x1,x2,myface_no)
                     call vmap(drem(:,:,1),x1,-x2,nbrface_no)
                   elseif (inbr==neast) then
                     x1 = elem%cartp(np,np)%x
                     x2 = elem%cartp(np,np)%y
                     call vmap(dloc(1,1,1),x1,x2,myface_no)
                     call vmap(drem(1,1,1),x1,-x2,nbrface_no)
                   endif
                   elseif (nbrface_no==2) then
                     if (inbr==east) then
                       do j = 1, np
                         jr = np+1-j
                         x1 = elem%cartp(np,j)%x
                         x2 = elem%cartp(np,j)%y
                         call vmap(dloc(1,1,jr),x1,x2,myface_no)
                         call vmap(drem(1,1,jr),-x2,-x1,nbrface_no)
                       enddo
                       elseif (inbr==seast) then
                         x1 = elem%cartp(np,1)%x
                         x2 = elem%cartp(np,1)%y
                         call vmap(dloc(1,1,1),x1,x2,myface_no)
                         call vmap(drem(1,1,1),-x2,-x1,nbrface_no)
                       elseif (inbr==neast) then
                         x1 = elem%cartp(np,np)%x
                         x2 = elem%cartp(np,np)%y
                         call vmap(dloc(1,1,1),x1,x2,myface_no)
                         call vmap(drem(1,1,1),-x2,-x1,nbrface_no)
                       endif
                       elseif (nbrface_no==3) then
                         if (inbr==south) then
                           do i = 1, np
                             ir = np+1-i
                             x1 = elem%cartp(i,1)%x
                             x2 = elem%cartp(i,1)%y
                             call vmap(dloc(1,1,ir),x1,x2,myface_no)
                             call vmap(drem(1,1,ir),-x1,x2,nbrface_no)
                           enddo
                           elseif (inbr==swest) then
                             x1 = elem%cartp(1,1)%x
                             x2 = elem%cartp(1,1)%y
                             call vmap(dloc(1,1,1),x1,x2,myface_no)
                             call vmap(drem(1,1,1),-x1,x2,nbrface_no)
                           elseif (inbr==seast) then
                             x1 = elem%cartp(np,1)%x
                             x2 = elem%cartp(np,1)%y
                             call vmap(dloc(1,1,1),x1,x2,myface_no)
                             call vmap(drem(1,1,1),-x1,x2,nbrface_no)
                           endif
                           elseif (nbrface_no==4) then
                             if (inbr==west) then
                               do j = 1, np
                                 x1 = elem%cartp(1,j)%x
                                 x2 = elem%cartp(1,j)%y
                                 call vmap(dloc(1,1,j),x1,x2,myface_no)
                                 call vmap(drem(1,1,j),x2,x1,nbrface_no)
                               enddo
                               elseif (inbr==swest) then
                                 x1 = elem%cartp(1,1)%x
                                 x2 = elem%cartp(1,1)%y
                                 call vmap(dloc(1,1,1),x1,x2,myface_no)
                                 call vmap(drem(1,1,1),x2,x1,nbrface_no)
                               elseif (inbr==nwest) then
                                 x1 = elem%cartp(1,np)%x
                                 x2 = elem%cartp(1,np)%y
                                 call vmap(dloc(1,1,1),x1,x2,myface_no)
                                 call vmap(drem(1,1,1),x2,x1,nbrface_no)
                               endif
                               endif
             endif
             elem%desc%rot(irot)%nbr = inbr
             if (rot_type=="covariant") then
               do i = 1, size(elem%desc%rot(irot)%r(:,:,:), 3)
                 elem%desc%rot(irot)%r(:,:,i) = covariant_rot(dloc(:,:,i),drem(:,:,i))
               enddo
               elseif (rot_type=="contravariant") then
                 do i = 1, size(elem%desc%rot(irot)%r(:,:,:), 3)
                   elem%desc%rot(irot)%r(:,:,i) = contravariant_rot(dloc(:,:,i),drem(:,:,i))
                 enddo
             endif
           endif
         enddo
         endif
     enddo
   endif
!
   end subroutine rotation_init_atomic
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine set_corner_coordinates(elem)
!-------------------------------------------------------------------------------
   use element, only: element_t
   use dimensions, only: ne
   use kiapsparallel, only: abortpar
!
   type(element_t) :: elem
! Local variables
   integer i, ie, je, face_no, nn
   real(real_kind) :: dx, dy, startx, starty
!if (0==ne) call AbortPar(message='Error in set_corner_coordinates: ne is zero')
!
   if (0==ne) call abortpar(message = 'Error in set_corner_coordinates:ne is zero')
! ========================================
! compute cube face coordinates of element
! =========================================
   call convert_gbl_index(elem%vertex%number,ie,je,face_no)
   elem%vertex%face_number = face_no
   dx =(cube_xend-cube_xstart)/ne
   dy =(cube_yend-cube_ystart)/ne
   startx = cube_xstart+ie*dx
   starty = cube_ystart+je*dy
   elem%corners(1)%x = startx
   elem%corners(1)%y = starty
   elem%corners(2)%x = startx+dx
   elem%corners(2)%y = starty
   elem%corners(3)%x = startx+dx
   elem%corners(3)%y = starty+dy
   elem%corners(4)%x = startx
   elem%corners(4)%y = starty+dy
   do i = 1,4
     elem%node_multiplicity(i) = 4
   enddo
   ie = ie+1
   je = je+1
   if (ie==1.and.je==1) then
     elem%node_multiplicity(1) = 3
   elseif (ie==ne.and.je==1) then
     elem%node_multiplicity(2) = 3
   elseif (ie==ne.and.je==ne) then
     elem%node_multiplicity(3) = 3
   elseif (ie==1.and.je==ne) then
     elem%node_multiplicity(4) = 3
   endif
!
   end subroutine set_corner_coordinates
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine assign_node_numbers_to_elem(elements, gridvertex)
!-------------------------------------------------------------------------------
   use dimensions, only: ne
   use element, only: element_t
   use control, only: north, south, east, west, neast, seast, swest, nwest
   use kiapsparallel, only: abortpar
   use gridgraph, only: gridvertex_t
   implicit none
!
   type(element_t), intent(inout) :: elements(:)
   type(gridvertex_t), intent(in   ) :: gridvertex(:)
   type(gridvertex_t) :: vertex
   integer :: connectivity(6*ne*ne, 4)
   integer :: nn(4), en(4)
   integer el, i, n, direction
   integer current_node_num, tot_ne
!
   current_node_num = 0
   tot_ne = 6*ne*ne
   if (0==ne) call abortpar(message = 'Error in assign_node_numbers_to_elem:ne is zero')
   if (tot_ne/= size(gridvertex)) call abortpar(message = 'Error in assign_node_numbers_to_elem:gridvertex not correct length')
   connectivity = 0
   do el = 1,tot_ne
     vertex = gridvertex(el)
     en = 0
     do direction = 1,8
       if (associated(vertex%nbrs(direction)%n)) then
         do i = 1, size(vertex%nbrs(direction)%n)
           n = vertex%nbrs(direction)%n(i)
           if (n/=-1) then
             nn = connectivity(n,:)
             select case(direction)
             case(north)
               if (nn(1)/= 0) en(4) = nn(1)
               if (nn(2)/= 0) en(3) = nn(2)
             case(south)
               if (nn(4)/= 0) en(1) = nn(4)
               if (nn(3)/= 0) en(2) = nn(3)
             case(east)
               if (nn(1)/= 0) en(2) = nn(1)
               if (nn(4)/= 0) en(3) = nn(4)
             case(west)
               if (nn(2)/= 0) en(1) = nn(2)
               if (nn(3)/= 0) en(4) = nn(3)
             case(neast)
               if (nn(1)/= 0) en(3) = nn(1)
             case(seast)
               if (nn(4)/= 0) en(2) = nn(4)
             case(swest)
               if (nn(3)/= 0) en(1) = nn(3)
             case(nwest)
               if (nn(2)/= 0) en(4) = nn(2)
             end select
           endif
         enddo
         endif
     enddo
     do i = 1,4
       if (en(i)==0) then
         current_node_num = current_node_num+1
         en(i) = current_node_num
       endif
     enddo
     connectivity(el,:) = en
   enddo
   if (current_node_num/=(6*ne*ne+2)) then
     call abortpar(message = 'Error in assignment of node numbers:failed euler test')
   endif
   do el = 1,size(elements)
     elements(el)%node_numbers = connectivity(elements(el)%vertex%number,:)
   enddo
!
   end subroutine assign_node_numbers_to_elem
! ================================================
! convert_gbl_index:
!
! Convert global element index to cube index
! ================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine convert_gbl_index(number, ie, je, face_no)
!-------------------------------------------------------------------------------
   use dimensions, only: ne
   use kiapsparallel, only: abortpar
!
   integer, intent(in   ) :: number
   integer, intent(  out) :: ie, je, face_no
!
   if (0==ne) call abortpar(message = 'Error in cube:convert_gbl_index:ne is zero')
!  inverse of the function:      number = 1 + ie + ne*je + ne*ne*(face_no-1)
   face_no =((number-1)/(ne*ne))+1
   ie = modulo(number-1,ne)
   je =(number-1)/ne-(face_no-1)*ne
!
   end subroutine convert_gbl_index
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine cubetopology(gridedge, gridvertex)
!-------------------------------------------------------------------------------
   use params, only: recursive, sfcurve
   use control, only: partmethod
   use gridgraph, only: gridedge_t, gridvertex_t, initgridedge
   use dimensions, only: np, ne
   use spacecurve, only: isfactorable, genspacecurve
   use control, only: north, south, east, west, neast, seast, swest, nwest
   use kiapsparallel, only: abortpar
!-----------------------
   implicit none
!
   type(gridedge_t), intent(inout), target :: gridedge(:)
   type(gridvertex_t), intent(inout), target :: gridvertex(:)
   integer, allocatable :: mesh(:,:)
   integer, allocatable :: mesh2(:,:), mesh2_map(:,:,:), sfcij(:,:)
   type(gridvertex_t), allocatable :: gridelem(:,:,:)
   integer :: i, j, k, l, number, irev, ne2, i2, j2, sfc_index
   integer :: edgewgtp, cornerwgt
   integer :: ielem, nedge
   integer :: offset, ierr, loc
   logical, allocatable :: nbrs_used(:,:,:,:)
!
   if (0==ne) call abortpar(message = 'Error in cubetopology:ne is zero')
   allocate(gridelem(ne,ne,nfaces),stat = ierr)
   if (ierr/= 0) then
     call abortpar(message = 'error in allocation of gridelem structure')
   endif
   allocate(nbrs_used(ne,ne,nfaces,8))
   nbrs_used = .false.
   number = 1
   edgewgtp = np
   cornerwgt = 1
   do k = 1,nfaces
     do j = 1, ne
       do i = 1, ne
! ====================================
! Number elements
! ====================================
! Do some initalization here
         do l = 1, 8
           nullify(gridelem(i,j,k)%nbrs(l)%n)
           nullify(gridelem(i,j,k)%nbrs(l)%f)
         enddo
           gridelem(i,j,k)%nbrs_ptr(:) = 0
           allocate(gridelem(i,j,k)%nbrs_face(8))
           gridelem(i,j,k)%wgtp(:) = 0
           gridelem(i,j,k)%wgtp_ghost(:) = 1 ! always this value
           gridelem(i,j,k)%spacecurve = 0
           gridelem(i,j,k)%number = number
           number = number+1
         enddo
         enddo
   enddo
!    print *,'CubeTopology: Ne,IsFactorable,IsLoadBalanced : ',ne,IsFactorable(ne),IsLoadBalanced(nelem,npart)
   allocate(mesh(ne,ne))
   if (isfactorable(ne)) then
     call genspacecurve(mesh)
!      call PrintCurve(Mesh)
   else
! find the smallest ne2 which is a power of 2 and ne2>ne
     ne2 = 2**ceiling(log(real(ne))/log(2d0))
     if (ne2<ne) call abortpar(message = 'Fatel sfc error')
     allocate(mesh2(ne2,ne2))
     allocate(mesh2_map(ne2,ne2,2))
     allocate(sfcij(0:ne2*ne2,2))
     call genspacecurve(mesh2) ! sfc partition for ne2
! associate every element on the ne x ne mesh (Mesh)
! with its closest element on the ne2 x ne2 mesh (Mesh2)
! Store this as a map from Mesh2 -> Mesh in Mesh2_map.
! elements in Mesh2 which are not mapped get assigned a value of 0
     mesh2_map = 0
     do j = 1,ne
       do i = 1, ne
! map this element to an (i2,j2) element
! [ (i-.5)/ne , (j-.5)/ne ]  = [ (i2-.5)/ne2 , (j2-.5)/ne2 ]
         i2 = nint(((i-.5)/ne)*ne2+.5)
         j2 = nint(((j-.5)/ne)*ne2+.5)
         if (i2<1) i2 = 1
         if (i2>ne2) i2 = ne2
         if (j2<1) j2 = 1
         if (j2>ne2) j2 = ne2
         mesh2_map(i2,j2,1) = i
         mesh2_map(i2,j2,2) = j
       enddo
     enddo
! create a reverse index array for Mesh2
! k = Mesh2(i,j)
! (i,j) = (sfcij(k,1),sfci(k,2))
     do j = 1,ne2
       do i = 1, ne2
         k = mesh2(i,j)
         sfcij(k,1) = i
         sfcij(k,2) = j
       enddo
     enddo
! generate a SFC for Mesh with the same ordering as the
! elements in Mesh2 which map to Mesh.
     sfc_index = 0
     do k = 0,ne2*ne2-1
       i2 = sfcij(k,1)
       j2 = sfcij(k,2)
       i = mesh2_map(i2,j2,1)
       j = mesh2_map(i2,j2,2)
       if (i/= 0) then
! (i2,j2) element maps to (i,j) element
         mesh(i,j) = sfc_index
         sfc_index = sfc_index+1
       endif
     enddo
#if 0
     print*,'SFC mapping to non powers of 2,3 used. mesh:'
     do j = 1,ne
       write(*,'(99i3) ')(mesh(i,j),i = 1,ne)
     enddo
     call printcurve(mesh2)
#endif
     deallocate(mesh2)
     deallocate(mesh2_map)
     deallocate(sfcij)
   endif
! -------------------------------------------
!  Setup the space-filling curve for face 1
! -------------------------------------------
   offset = 0
   do j = 1,ne
     do i = 1, ne
       gridelem(i,j,1)%spacecurve = offset+mesh(i,ne-j+1)
     enddo
   enddo
! -------------------------------------------
!  Setup the space-filling curve for face 2
! -------------------------------------------
   offset = offset+ne*ne
   do j = 1,ne
     do i = 1, ne
       gridelem(i,j,2)%spacecurve = offset+mesh(i,ne-j+1)
     enddo
   enddo
! -------------------------------------------
!  Setup the space-filling curve for face 6
! -------------------------------------------
   offset = offset+ne*ne
   do j = 1,ne
     do i = 1, ne
       gridelem(i,j,6)%spacecurve = offset+mesh(ne-i+1,ne-j+1)
     enddo
   enddo
! -------------------------------------------
!  Setup the space-filling curve for face 4
! -------------------------------------------
   offset = offset+ne*ne
   do j = 1,ne
     do i = 1, ne
       gridelem(i,j,4)%spacecurve = offset+mesh(ne-j+1,i)
     enddo
   enddo
! -------------------------------------------
!  Setup the space-filling curve for face 5
! -------------------------------------------
   offset = offset+ne*ne
   do j = 1,ne
     do i = 1, ne
       gridelem(i,j,5)%spacecurve = offset+mesh(i,j)
     enddo
   enddo
! -------------------------------------------
!  Setup the space-filling curve for face 3
! -------------------------------------------
   offset = offset+ne*ne
   do j = 1,ne
     do i = 1, ne
       gridelem(i,j,3)%spacecurve = offset+mesh(i,j)
     enddo
   enddo
! ==================
! face interiors
! ==================
   do k = 1,6
! setup  SOUTH, WEST, SW neighbors
     do j = 2, ne
       do i = 2, ne
         nbrs_used(i,j,k,west) = .true.
         nbrs_used(i,j,k,south) = .true.
         nbrs_used(i,j,k,swest) = .true.
         allocate(gridelem(i,j,k)%nbrs(west)%n(1))
         allocate(gridelem(i,j,k)%nbrs(south)%n(1))
         allocate(gridelem(i,j,k)%nbrs(swest)%n(1))
         allocate(gridelem(i,j,k)%nbrs(west)%f(1))
         allocate(gridelem(i,j,k)%nbrs(south)%f(1))
         allocate(gridelem(i,j,k)%nbrs(swest)%f(1))
         gridelem(i,j,k)%nbrs(west)%n(1) = gridelem(i-1,j,k)%number
         gridelem(i,j,k)%nbrs(west)%f(1) = k
         gridelem(i,j,k)%wgtp(west) = edgewgtp
         gridelem(i,j,k)%nbrs(south)%n(1) = gridelem(i,j-1,k)%number
         gridelem(i,j,k)%nbrs(south)%f(1) = k
         gridelem(i,j,k)%wgtp(south) = edgewgtp
         gridelem(i,j,k)%nbrs(swest)%n(1) = gridelem(i-1,j-1,k)%number
         gridelem(i,j,k)%nbrs(swest)%f(1) = k
         gridelem(i,j,k)%wgtp(swest) = cornerwgt
       enddo
       enddo
!  setup EAST, NORTH, NE neighbors
         do j = 1, ne-1
           do i = 1, ne-1
             nbrs_used(i,j,k,east) = .true.
             nbrs_used(i,j,k,north) = .true.
             nbrs_used(i,j,k,neast) = .true.
             allocate(gridelem(i,j,k)%nbrs(east)%n(1))
             allocate(gridelem(i,j,k)%nbrs(north)%n(1))
             allocate(gridelem(i,j,k)%nbrs(neast)%n(1))
             allocate(gridelem(i,j,k)%nbrs(east)%f(1))
             allocate(gridelem(i,j,k)%nbrs(north)%f(1))
             allocate(gridelem(i,j,k)%nbrs(neast)%f(1))
             gridelem(i,j,k)%nbrs(east)%n(1) = gridelem(i+1,j,k)%number
             gridelem(i,j,k)%nbrs(east)%f(1) = k
             gridelem(i,j,k)%wgtp(east) = edgewgtp
             gridelem(i,j,k)%nbrs(north)%n(1) = gridelem(i,j+1,k)%number
             gridelem(i,j,k)%nbrs(north)%f(1) = k
             gridelem(i,j,k)%wgtp(north) = edgewgtp
             gridelem(i,j,k)%nbrs(neast)%n(1) = gridelem(i+1,j+1,k)%number
             gridelem(i,j,k)%nbrs(neast)%f(1) = k
             gridelem(i,j,k)%wgtp(neast) = cornerwgt
           enddo
           enddo
! Setup the remaining SOUTH, EAST, and SE neighbors
             do j = 2, ne
               do i = 1, ne-1
                 nbrs_used(i,j,k,south) = .true.
                 nbrs_used(i,j,k,east) = .true.
                 nbrs_used(i,j,k,seast) = .true.
                 allocate(gridelem(i,j,k)%nbrs(south)%n(1))
                 allocate(gridelem(i,j,k)%nbrs(east)%n(1))
                 allocate(gridelem(i,j,k)%nbrs(seast)%n(1))
                 allocate(gridelem(i,j,k)%nbrs(south)%f(1))
                 allocate(gridelem(i,j,k)%nbrs(east)%f(1))
                 allocate(gridelem(i,j,k)%nbrs(seast)%f(1))
                 gridelem(i,j,k)%nbrs(south)%n(1) = gridelem(i,j-1,k)%number
                 gridelem(i,j,k)%nbrs(south)%f(1) = k
                 gridelem(i,j,k)%wgtp(south) = edgewgtp
                 gridelem(i,j,k)%nbrs(east)%n(1) = gridelem(i+1,j,k)%number
                 gridelem(i,j,k)%nbrs(east)%f(1) = k
                 gridelem(i,j,k)%wgtp(east) = edgewgtp
                 gridelem(i,j,k)%nbrs(seast)%n(1) = gridelem(i+1,j-1,k)%number
                 gridelem(i,j,k)%nbrs(seast)%f(1) = k
                 gridelem(i,j,k)%wgtp(seast) = cornerwgt
               enddo
               enddo
! Setup the remaining NORTH, WEST, and NW neighbors
                 do j = 1, ne-1
                   do i = 2, ne
                     nbrs_used(i,j,k,north) = .true.
                     nbrs_used(i,j,k,west) = .true.
                     nbrs_used(i,j,k,nwest) = .true.
                     allocate(gridelem(i,j,k)%nbrs(north)%n(1))
                     allocate(gridelem(i,j,k)%nbrs(west)%n(1))
                     allocate(gridelem(i,j,k)%nbrs(nwest)%n(1))
                     allocate(gridelem(i,j,k)%nbrs(north)%f(1))
                     allocate(gridelem(i,j,k)%nbrs(west)%f(1))
                     allocate(gridelem(i,j,k)%nbrs(nwest)%f(1))
                     gridelem(i,j,k)%nbrs(north)%n(1) = gridelem(i,j+1,k)%number
                     gridelem(i,j,k)%nbrs(north)%f(1) = k
                     gridelem(i,j,k)%wgtp(north) = edgewgtp
                     gridelem(i,j,k)%nbrs(west)%n(1) = gridelem(i-1,j,k)%number
                     gridelem(i,j,k)%nbrs(west)%f(1) = k
                     gridelem(i,j,k)%wgtp(west) = edgewgtp
                     gridelem(i,j,k)%nbrs(nwest)%n(1) = gridelem(i-1,j+1,k)%number
                     gridelem(i,j,k)%nbrs(nwest)%f(1) = k
                     gridelem(i,j,k)%wgtp(nwest) = cornerwgt
                   enddo
                   enddo
   enddo
! ======================
! west/east "belt" edges
! ======================
   do k = 1,4
     do j = 1, ne
       nbrs_used(1,j,k,west) = .true.
       nbrs_used(ne,j,k,east) = .true.
       allocate(gridelem(1,j,k)%nbrs(west)%n(1))
       allocate(gridelem(ne,j,k)%nbrs(east)%n(1))
       allocate(gridelem(1,j,k)%nbrs(west)%f(1))
       allocate(gridelem(ne,j,k)%nbrs(east)%f(1))
       gridelem(1,j,k)%nbrs(west)%n(1) = gridelem(ne,j,modulo(2+k,4)+1)%number
       gridelem(1,j,k)%nbrs(west)%f(1) = modulo(2+k,4)+1
       gridelem(1,j,k)%wgtp(west) = edgewgtp
       gridelem(ne,j,k)%nbrs(east)%n(1) = gridelem(1,j,modulo(k,4)+1)%number
       gridelem(ne,j,k)%nbrs(east)%f(1) = modulo(k,4)+1
       gridelem(ne,j,k)%wgtp(east) = edgewgtp
!  Special rules for corner 'edges'
       if (j/= 1) then
         nbrs_used(1,j,k,swest) = .true.
         nbrs_used(ne,j,k,seast) = .true.
         allocate(gridelem(1,j,k)%nbrs(swest)%n(1))
         allocate(gridelem(ne,j,k)%nbrs(seast)%n(1))
         allocate(gridelem(1,j,k)%nbrs(swest)%f(1))
         allocate(gridelem(ne,j,k)%nbrs(seast)%f(1))
         gridelem(1,j,k)%nbrs(swest)%n(1) = gridelem(ne,j-1,modulo(2+k,4)+1)%number
         gridelem(1,j,k)%nbrs(swest)%f(1) = modulo(2+k,4)+1
         gridelem(1,j,k)%wgtp(swest) = cornerwgt
         gridelem(ne,j,k)%nbrs(seast)%n(1) = gridelem(1,j-1,modulo(k,4)+1)%number
         gridelem(ne,j,k)%nbrs(seast)%f(1) = modulo(k,4)+1
         gridelem(ne,j,k)%wgtp(seast) = cornerwgt
       endif
       if (j/= ne) then
         nbrs_used(1,j,k,nwest) = .true.
         nbrs_used(ne,j,k,neast) = .true.
         allocate(gridelem(1,j,k)%nbrs(nwest)%n(1))
         allocate(gridelem(ne,j,k)%nbrs(neast)%n(1))
         allocate(gridelem(1,j,k)%nbrs(nwest)%f(1))
         allocate(gridelem(ne,j,k)%nbrs(neast)%f(1))
         gridelem(1,j,k)%nbrs(nwest)%n(1) = gridelem(ne,j+1,modulo(2+k,4)+1)%number
         gridelem(1,j,k)%nbrs(nwest)%f(1) = modulo(2+k,4)+1
         gridelem(1,j,k)%wgtp(nwest) = cornerwgt
         gridelem(ne,j,k)%nbrs(neast)%n(1) = gridelem(1,j+1,modulo(k,4)+1)%number
         gridelem(ne,j,k)%nbrs(neast)%f(1) = modulo(k,4)+1
         gridelem(ne,j,k)%wgtp(neast) = cornerwgt
       endif
     enddo
   enddo
! ==================================
! south edge of 1 / north edge of 5
! ==================================
   do i = 1,ne
     nbrs_used(i,1,1,south) = .true.
     nbrs_used(i,ne,5,north) = .true.
     allocate(gridelem(i,1,1)%nbrs(south)%n(1))
     allocate(gridelem(i,ne,5)%nbrs(north)%n(1))
     allocate(gridelem(i,1,1)%nbrs(south)%f(1))
     allocate(gridelem(i,ne,5)%nbrs(north)%f(1))
     gridelem(i,1,1)%nbrs(south)%n(1) = gridelem(i,ne,5)%number
     gridelem(i,1,1)%nbrs(south)%f(1) = 5
     gridelem(i,1,1)%wgtp(south) = edgewgtp
     gridelem(i,ne,5)%nbrs(north)%n(1) = gridelem(i,1,1)%number
     gridelem(i,ne,5)%nbrs(north)%f(1) = 1
     gridelem(i,ne,5)%wgtp(north) = edgewgtp
!  Special rules for corner 'edges'
     if (i/= 1) then
       nbrs_used(i,1,1,swest) = .true.
       nbrs_used(i,ne,5,nwest) = .true.
       allocate(gridelem(i,1,1)%nbrs(swest)%n(1))
       allocate(gridelem(i,ne,5)%nbrs(nwest)%n(1))
       allocate(gridelem(i,1,1)%nbrs(swest)%f(1))
       allocate(gridelem(i,ne,5)%nbrs(nwest)%f(1))
       gridelem(i,1,1)%nbrs(swest)%n(1) = gridelem(i-1,ne,5)%number
       gridelem(i,1,1)%nbrs(swest)%f(1) = 5
       gridelem(i,1,1)%wgtp(swest) = cornerwgt
       gridelem(i,ne,5)%nbrs(nwest)%n(1) = gridelem(i-1,1,1)%number
       gridelem(i,ne,5)%nbrs(nwest)%f(1) = 1
       gridelem(i,ne,5)%wgtp(nwest) = cornerwgt
     endif
     if (i/= ne) then
       nbrs_used(i,1,1,seast) = .true.
       nbrs_used(i,ne,5,neast) = .true.
       allocate(gridelem(i,1,1)%nbrs(seast)%n(1))
       allocate(gridelem(i,ne,5)%nbrs(neast)%n(1))
       allocate(gridelem(i,1,1)%nbrs(seast)%f(1))
       allocate(gridelem(i,ne,5)%nbrs(neast)%f(1))
       gridelem(i,1,1)%nbrs(seast)%n(1) = gridelem(i+1,ne,5)%number
       gridelem(i,1,1)%nbrs(seast)%f(1) = 5
       gridelem(i,1,1)%wgtp(seast) = cornerwgt
       gridelem(i,ne,5)%nbrs(neast)%n(1) = gridelem(i+1,1,1)%number
       gridelem(i,ne,5)%nbrs(neast)%f(1) = 1
       gridelem(i,ne,5)%wgtp(neast) = cornerwgt
     endif
   enddo
! ==================================
! south edge of 2 / east edge of 5
! ==================================
   do i = 1,ne
     irev = ne+1-i
     nbrs_used(i,1,2,south) = .true.
     nbrs_used(ne,i,5,east) = .true.
     allocate(gridelem(i,1,2)%nbrs(south)%n(1))
     allocate(gridelem(ne,i,5)%nbrs(east)%n(1))
     allocate(gridelem(i,1,2)%nbrs(south)%f(1))
     allocate(gridelem(ne,i,5)%nbrs(east)%f(1))
     gridelem(i,1,2)%nbrs(south)%n(1) = gridelem(ne,irev,5)%number
     gridelem(i,1,2)%nbrs(south)%f(1) = 5
     gridelem(i,1,2)%wgtp(south) = edgewgtp
     gridelem(ne,i,5)%nbrs(east)%n(1) = gridelem(irev,1,2)%number
     gridelem(ne,i,5)%nbrs(east)%f(1) = 2
     gridelem(ne,i,5)%wgtp(east) = edgewgtp
!  Special rules for corner 'edges'
     if (i/= 1) then
       nbrs_used(i,1,2,swest) = .true.
       nbrs_used(ne,i,5,seast) = .true.
       allocate(gridelem(i,1,2)%nbrs(swest)%n(1))
       allocate(gridelem(ne,i,5)%nbrs(seast)%n(1))
       allocate(gridelem(i,1,2)%nbrs(swest)%f(1))
       allocate(gridelem(ne,i,5)%nbrs(seast)%f(1))
       gridelem(i,1,2)%nbrs(swest)%n(1) = gridelem(ne,irev+1,5)%number
       gridelem(i,1,2)%nbrs(swest)%f(1) = 5
       gridelem(i,1,2)%wgtp(swest) = cornerwgt
       gridelem(ne,i,5)%nbrs(seast)%n(1) = gridelem(irev+1,1,2)%number
       gridelem(ne,i,5)%nbrs(seast)%f(1) = 2
       gridelem(ne,i,5)%wgtp(seast) = cornerwgt
     endif
     if (i/= ne) then
       nbrs_used(i,1,2,seast) = .true.
       nbrs_used(ne,i,5,neast) = .true.
       allocate(gridelem(i,1,2)%nbrs(seast)%n(1))
       allocate(gridelem(ne,i,5)%nbrs(neast)%n(1))
       allocate(gridelem(i,1,2)%nbrs(seast)%f(1))
       allocate(gridelem(ne,i,5)%nbrs(neast)%f(1))
       gridelem(i,1,2)%nbrs(seast)%n(1) = gridelem(ne,irev-1,5)%number
       gridelem(i,1,2)%nbrs(seast)%f(1) = 5
       gridelem(i,1,2)%wgtp(seast) = cornerwgt
       gridelem(ne,i,5)%nbrs(neast)%n(1) = gridelem(irev-1,1,2)%number
       gridelem(ne,i,5)%nbrs(neast)%f(1) = 2
       gridelem(ne,i,5)%wgtp(neast) = cornerwgt
     endif
   enddo
! ==================================
! south edge of 3 / south edge of 5
! ==================================
   do i = 1,ne
     irev = ne+1-i
     nbrs_used(i,1,3,south) = .true.
     nbrs_used(i,1,5,south) = .true.
     allocate(gridelem(i,1,3)%nbrs(south)%n(1))
     allocate(gridelem(i,1,5)%nbrs(south)%n(1))
     allocate(gridelem(i,1,3)%nbrs(south)%f(1))
     allocate(gridelem(i,1,5)%nbrs(south)%f(1))
     gridelem(i,1,3)%nbrs(south)%n(1) = gridelem(irev,1,5)%number
     gridelem(i,1,3)%nbrs(south)%f(1) = 5
     gridelem(i,1,3)%wgtp(south) = edgewgtp
     gridelem(i,1,5)%nbrs(south)%n(1) = gridelem(irev,1,3)%number
     gridelem(i,1,5)%nbrs(south)%f(1) = 3
     gridelem(i,1,5)%wgtp(south) = edgewgtp
!  Special rules for corner 'edges'
     if (i/= 1) then
       nbrs_used(i,1,3,swest) = .true.
       nbrs_used(i,1,5,swest) = .true.
       allocate(gridelem(i,1,3)%nbrs(swest)%n(1))
       allocate(gridelem(i,1,5)%nbrs(swest)%n(1))
       allocate(gridelem(i,1,3)%nbrs(swest)%f(1))
       allocate(gridelem(i,1,5)%nbrs(swest)%f(1))
       gridelem(i,1,3)%nbrs(swest)%n(1) = gridelem(irev+1,1,5)%number
       gridelem(i,1,3)%nbrs(swest)%f(1) = 5
       gridelem(i,1,3)%wgtp(swest) = cornerwgt
       gridelem(i,1,5)%nbrs(swest)%n(1) = gridelem(irev+1,1,3)%number
       gridelem(i,1,5)%nbrs(swest)%f(1) = 3
       gridelem(i,1,5)%wgtp(swest) = cornerwgt
     endif
     if (i/= ne) then
       nbrs_used(i,1,3,seast) = .true.
       nbrs_used(i,1,5,seast) = .true.
       allocate(gridelem(i,1,3)%nbrs(seast)%n(1))
       allocate(gridelem(i,1,5)%nbrs(seast)%n(1))
       allocate(gridelem(i,1,3)%nbrs(seast)%f(1))
       allocate(gridelem(i,1,5)%nbrs(seast)%f(1))
       gridelem(i,1,3)%nbrs(seast)%n(1) = gridelem(irev-1,1,5)%number
       gridelem(i,1,3)%nbrs(seast)%f(1) = 5
       gridelem(i,1,3)%wgtp(seast) = cornerwgt
       gridelem(i,1,5)%nbrs(seast)%n(1) = gridelem(irev-1,1,3)%number
       gridelem(i,1,5)%nbrs(seast)%f(1) = 3
       gridelem(i,1,5)%wgtp(seast) = cornerwgt
     endif
   enddo
! ==================================
! south edge of 4 / west edge of 5
! ==================================
   do i = 1,ne
     irev = ne+1-i
     nbrs_used(i,1,4,south) = .true.
     nbrs_used(1,i,5,west) = .true.
     allocate(gridelem(i,1,4)%nbrs(south)%n(1))
     allocate(gridelem(1,i,5)%nbrs(west)%n(1))
     allocate(gridelem(i,1,4)%nbrs(south)%f(1))
     allocate(gridelem(1,i,5)%nbrs(west)%f(1))
     gridelem(i,1,4)%nbrs(south)%n(1) = gridelem(1,i,5)%number
     gridelem(i,1,4)%nbrs(south)%f(1) = 5
     gridelem(i,1,4)%wgtp(south) = edgewgtp
     gridelem(1,i,5)%nbrs(west)%n(1) = gridelem(i,1,4)%number
     gridelem(1,i,5)%nbrs(west)%f(1) = 4
     gridelem(1,i,5)%wgtp(west) = edgewgtp
!  Special rules for corner 'edges'
     if (i/= 1) then
       nbrs_used(i,1,4,swest) = .true.
       nbrs_used(1,i,5,swest) = .true.
       allocate(gridelem(i,1,4)%nbrs(swest)%n(1))
       allocate(gridelem(1,i,5)%nbrs(swest)%n(1))
       allocate(gridelem(i,1,4)%nbrs(swest)%f(1))
       allocate(gridelem(1,i,5)%nbrs(swest)%f(1))
       gridelem(i,1,4)%nbrs(swest)%n(1) = gridelem(1,i-1,5)%number
       gridelem(i,1,4)%nbrs(swest)%f(1) = 5
       gridelem(i,1,4)%wgtp(swest) = cornerwgt
       gridelem(1,i,5)%nbrs(swest)%n(1) = gridelem(i-1,1,4)%number
       gridelem(1,i,5)%nbrs(swest)%f(1) = 4
       gridelem(1,i,5)%wgtp(swest) = cornerwgt
     endif
     if (i/= ne) then
       nbrs_used(i,1,4,seast) = .true.
       nbrs_used(1,i,5,nwest) = .true.
       allocate(gridelem(i,1,4)%nbrs(seast)%n(1))
       allocate(gridelem(1,i,5)%nbrs(nwest)%n(1))
       allocate(gridelem(i,1,4)%nbrs(seast)%f(1))
       allocate(gridelem(1,i,5)%nbrs(nwest)%f(1))
       gridelem(i,1,4)%nbrs(seast)%n(1) = gridelem(1,i+1,5)%number
       gridelem(i,1,4)%nbrs(seast)%f(1) = 5
       gridelem(i,1,4)%wgtp(seast) = cornerwgt
       gridelem(1,i,5)%nbrs(nwest)%n(1) = gridelem(i+1,1,4)%number
       gridelem(1,i,5)%nbrs(nwest)%f(1) = 4
       gridelem(1,i,5)%wgtp(nwest) = cornerwgt
     endif
   enddo
! ==================================
! north edge of 1 / south edge of 6
! ==================================
   do i = 1,ne
     nbrs_used(i,ne,1,north) = .true.
     nbrs_used(i,1,6,south) = .true.
     allocate(gridelem(i,ne,1)%nbrs(north)%n(1))
     allocate(gridelem(i,1,6)%nbrs(south)%n(1))
     allocate(gridelem(i,ne,1)%nbrs(north)%f(1))
     allocate(gridelem(i,1,6)%nbrs(south)%f(1))
     gridelem(i,ne,1)%nbrs(north)%n(1) = gridelem(i,1,6)%number
     gridelem(i,ne,1)%nbrs(north)%f(1) = 6
     gridelem(i,ne,1)%wgtp(north) = edgewgtp
     gridelem(i,1,6)%nbrs(south)%n(1) = gridelem(i,ne,1)%number
     gridelem(i,1,6)%nbrs(south)%f(1) = 1
     gridelem(i,1,6)%wgtp(south) = edgewgtp
!  Special rules for corner 'edges'
     if (i/= 1) then
       nbrs_used(i,ne,1,nwest) = .true.
       nbrs_used(i,1,6,swest) = .true.
       allocate(gridelem(i,ne,1)%nbrs(nwest)%n(1))
       allocate(gridelem(i,1,6)%nbrs(swest)%n(1))
       allocate(gridelem(i,ne,1)%nbrs(nwest)%f(1))
       allocate(gridelem(i,1,6)%nbrs(swest)%f(1))
       gridelem(i,ne,1)%nbrs(nwest)%n(1) = gridelem(i-1,1,6)%number
       gridelem(i,ne,1)%nbrs(nwest)%f(1) = 6
       gridelem(i,ne,1)%wgtp(nwest) = cornerwgt
       gridelem(i,1,6)%nbrs(swest)%n(1) = gridelem(i-1,ne,1)%number
       gridelem(i,1,6)%nbrs(swest)%f(1) = 1
       gridelem(i,1,6)%wgtp(swest) = cornerwgt
     endif
     if (i/= ne) then
       nbrs_used(i,ne,1,neast) = .true.
       nbrs_used(i,1,6,seast) = .true.
       allocate(gridelem(i,ne,1)%nbrs(neast)%n(1))
       allocate(gridelem(i,1,6)%nbrs(seast)%n(1))
       allocate(gridelem(i,ne,1)%nbrs(neast)%f(1))
       allocate(gridelem(i,1,6)%nbrs(seast)%f(1))
       gridelem(i,ne,1)%nbrs(neast)%n(1) = gridelem(i+1,1,6)%number
       gridelem(i,ne,1)%nbrs(neast)%f(1) = 6
       gridelem(i,ne,1)%wgtp(neast) = cornerwgt
       gridelem(i,1,6)%nbrs(seast)%n(1) = gridelem(i+1,ne,1)%number
       gridelem(i,1,6)%nbrs(seast)%f(1) = 1
       gridelem(i,1,6)%wgtp(seast) = cornerwgt
     endif
   enddo
! ==================================
! north edge of 2 / east edge of 6
! ==================================
   do i = 1,ne
     nbrs_used(i,ne,2,north) = .true.
     nbrs_used(ne,i,6,east) = .true.
     allocate(gridelem(i,ne,2)%nbrs(north)%n(1))
     allocate(gridelem(ne,i,6)%nbrs(east)%n(1))
     allocate(gridelem(i,ne,2)%nbrs(north)%f(1))
     allocate(gridelem(ne,i,6)%nbrs(east)%f(1))
     gridelem(i,ne,2)%nbrs(north)%n(1) = gridelem(ne,i,6)%number
     gridelem(i,ne,2)%nbrs(north)%f(1) = 6
     gridelem(i,ne,2)%wgtp(north) = edgewgtp
     gridelem(ne,i,6)%nbrs(east)%n(1) = gridelem(i,ne,2)%number
     gridelem(ne,i,6)%nbrs(east)%f(1) = 2
     gridelem(ne,i,6)%wgtp(east) = edgewgtp
!  Special rules for corner 'edges'
     if (i/= 1) then
       nbrs_used(i,ne,2,nwest) = .true.
       nbrs_used(ne,i,6,seast) = .true.
       allocate(gridelem(i,ne,2)%nbrs(nwest)%n(1))
       allocate(gridelem(ne,i,6)%nbrs(seast)%n(1))
       allocate(gridelem(i,ne,2)%nbrs(nwest)%f(1))
       allocate(gridelem(ne,i,6)%nbrs(seast)%f(1))
       gridelem(i,ne,2)%nbrs(nwest)%n(1) = gridelem(ne,i-1,6)%number
       gridelem(i,ne,2)%nbrs(nwest)%f(1) = 6
       gridelem(i,ne,2)%wgtp(nwest) = cornerwgt
       gridelem(ne,i,6)%nbrs(seast)%n(1) = gridelem(i-1,ne,2)%number
       gridelem(ne,i,6)%nbrs(seast)%f(1) = 2
       gridelem(ne,i,6)%wgtp(seast) = cornerwgt
     endif
     if (i/= ne) then
       nbrs_used(i,ne,2,neast) = .true.
       nbrs_used(ne,i,6,neast) = .true.
       allocate(gridelem(i,ne,2)%nbrs(neast)%n(1))
       allocate(gridelem(ne,i,6)%nbrs(neast)%n(1))
       allocate(gridelem(i,ne,2)%nbrs(neast)%f(1))
       allocate(gridelem(ne,i,6)%nbrs(neast)%f(1))
       gridelem(i,ne,2)%nbrs(neast)%n(1) = gridelem(ne,i+1,6)%number
       gridelem(i,ne,2)%nbrs(neast)%f(1) = 6
       gridelem(i,ne,2)%wgtp(neast) = cornerwgt
       gridelem(ne,i,6)%nbrs(neast)%n(1) = gridelem(i+1,ne,2)%number
       gridelem(ne,i,6)%nbrs(neast)%f(1) = 2
       gridelem(ne,i,6)%wgtp(neast) = cornerwgt
     endif
   enddo
! ===================================
! north edge of 3 / north edge of 6
! ===================================
   do i = 1,ne
     irev = ne+1-i
     nbrs_used(i,ne,3,north) = .true.
     nbrs_used(i,ne,6,north) = .true.
     allocate(gridelem(i,ne,3)%nbrs(north)%n(1))
     allocate(gridelem(i,ne,6)%nbrs(north)%n(1))
     allocate(gridelem(i,ne,3)%nbrs(north)%f(1))
     allocate(gridelem(i,ne,6)%nbrs(north)%f(1))
     gridelem(i,ne,3)%nbrs(north)%n(1) = gridelem(irev,ne,6)%number
     gridelem(i,ne,3)%nbrs(north)%f(1) = 6
     gridelem(i,ne,3)%wgtp(north) = edgewgtp
     gridelem(i,ne,6)%nbrs(north)%n(1) = gridelem(irev,ne,3)%number
     gridelem(i,ne,6)%nbrs(north)%f(1) = 3
     gridelem(i,ne,6)%wgtp(north) = edgewgtp
!  Special rules for corner 'edges'
     if (i/= 1) then
       nbrs_used(i,ne,3,nwest) = .true.
       nbrs_used(i,ne,6,nwest) = .true.
       allocate(gridelem(i,ne,3)%nbrs(nwest)%n(1))
       allocate(gridelem(i,ne,6)%nbrs(nwest)%n(1))
       allocate(gridelem(i,ne,3)%nbrs(nwest)%f(1))
       allocate(gridelem(i,ne,6)%nbrs(nwest)%f(1))
       gridelem(i,ne,3)%nbrs(nwest)%n(1) = gridelem(irev+1,ne,6)%number
       gridelem(i,ne,3)%nbrs(nwest)%f(1) = 6
       gridelem(i,ne,3)%wgtp(nwest) = cornerwgt
       gridelem(i,ne,6)%nbrs(nwest)%n(1) = gridelem(irev+1,ne,3)%number
       gridelem(i,ne,6)%nbrs(nwest)%f(1) = 3
       gridelem(i,ne,6)%wgtp(nwest) = cornerwgt
     endif
     if (i/= ne) then
       nbrs_used(i,ne,3,neast) = .true.
       nbrs_used(i,ne,6,neast) = .true.
       allocate(gridelem(i,ne,3)%nbrs(neast)%n(1))
       allocate(gridelem(i,ne,6)%nbrs(neast)%n(1))
       allocate(gridelem(i,ne,3)%nbrs(neast)%f(1))
       allocate(gridelem(i,ne,6)%nbrs(neast)%f(1))
       gridelem(i,ne,3)%nbrs(neast)%n(1) = gridelem(irev-1,ne,6)%number
       gridelem(i,ne,3)%nbrs(neast)%f(1) = 6
       gridelem(i,ne,3)%wgtp(neast) = cornerwgt
       gridelem(i,ne,6)%nbrs(neast)%n(1) = gridelem(irev-1,ne,3)%number
       gridelem(i,ne,6)%nbrs(neast)%f(1) = 3
       gridelem(i,ne,6)%wgtp(neast) = cornerwgt
     endif
   enddo
! ===================================
! north edge of 4 / west edge of 6
! ===================================
   do i = 1,ne
     irev = ne+1-i
     nbrs_used(i,ne,4,north) = .true.
     nbrs_used(1,i,6,west) = .true.
     allocate(gridelem(i,ne,4)%nbrs(north)%n(1))
     allocate(gridelem(1,i,6)%nbrs(west)%n(1))
     allocate(gridelem(i,ne,4)%nbrs(north)%f(1))
     allocate(gridelem(1,i,6)%nbrs(west)%f(1))
     gridelem(i,ne,4)%nbrs(north)%n(1) = gridelem(1,irev,6)%number
     gridelem(i,ne,4)%nbrs(north)%f(1) = 6
     gridelem(i,ne,4)%wgtp(north) = edgewgtp
     gridelem(1,i,6)%nbrs(west)%n(1) = gridelem(irev,ne,4)%number
     gridelem(1,i,6)%nbrs(west)%f(1) = 4
     gridelem(1,i,6)%wgtp(west) = edgewgtp
!  Special rules for corner 'edges'
     if (i/= 1) then
       nbrs_used(i,ne,4,nwest) = .true.
       nbrs_used(1,i,6,swest) = .true.
       allocate(gridelem(i,ne,4)%nbrs(nwest)%n(1))
       allocate(gridelem(1,i,6)%nbrs(swest)%n(1))
       allocate(gridelem(i,ne,4)%nbrs(nwest)%f(1))
       allocate(gridelem(1,i,6)%nbrs(swest)%f(1))
       gridelem(i,ne,4)%nbrs(nwest)%n(1) = gridelem(1,irev+1,6)%number
       gridelem(i,ne,4)%nbrs(nwest)%f(1) = 6
       gridelem(i,ne,4)%wgtp(nwest) = cornerwgt
       gridelem(1,i,6)%nbrs(swest)%n(1) = gridelem(irev+1,ne,4)%number
       gridelem(1,i,6)%nbrs(swest)%f(1) = 4
       gridelem(1,i,6)%wgtp(swest) = cornerwgt
     endif
     if (i/= ne) then
       nbrs_used(i,ne,4,neast) = .true.
       nbrs_used(1,i,6,nwest) = .true.
       allocate(gridelem(i,ne,4)%nbrs(neast)%n(1))
       allocate(gridelem(1,i,6)%nbrs(nwest)%n(1))
       allocate(gridelem(i,ne,4)%nbrs(neast)%f(1))
       allocate(gridelem(1,i,6)%nbrs(nwest)%f(1))
       gridelem(i,ne,4)%nbrs(neast)%n(1) = gridelem(1,irev-1,6)%number
       gridelem(i,ne,4)%nbrs(neast)%f(1) = 6
       gridelem(i,ne,4)%wgtp(neast) = cornerwgt
       gridelem(1,i,6)%nbrs(nwest)%n(1) = gridelem(irev-1,ne,4)%number
       gridelem(1,i,6)%nbrs(nwest)%f(1) = 4
       gridelem(1,i,6)%wgtp(nwest) = cornerwgt
     endif
   enddo
   ielem = 1 ! element counter
   do k = 1,6
     do j = 1, ne
       do i = 1, ne
         gridvertex(ielem)%nbrs_ptr(1) = 1
         do l = 1,8
           if (associated(gridelem(i, j, k)%nbrs(l)%n)) then
             allocate(gridvertex(ielem)%nbrs(l)%n(size(gridelem(i,j,k)%nbrs(l)%n)))
             allocate(gridvertex(ielem)%nbrs(l)%f(size(gridelem(i,j,k)%nbrs(l)%n)))
             gridvertex(ielem)%nbrs(l)%n(:) = gridelem(i,j,k)%nbrs(l)%n(:)
             gridvertex(ielem)%nbrs(l)%f(:) = gridelem(i,j,k)%nbrs(l)%f(:)
           endif
! for VLAG
             loc = gridvertex(ielem)%nbrs_ptr(l)
             if (nbrs_used(i,j,k,l)) then
               gridvertex(ielem)%nbrs_ptr(l+1) = gridvertex(ielem)%nbrs_ptr(l)+1
             else
               gridvertex(ielem)%nbrs_ptr(l+1) = gridvertex(ielem)%nbrs_ptr(l)
             endif
         enddo
         gridvertex(ielem)%wgtp = gridelem(i,j,k)%wgtp
         gridvertex(ielem)%wgtp_ghost = gridelem(i,j,k)%wgtp_ghost
         gridvertex(ielem)%number = gridelem(i,j,k)%number
         gridvertex(ielem)%processor_number = 0
         gridvertex(ielem)%spacecurve = gridelem(i,j,k)%spacecurve
         ielem = ielem+1
       enddo
       enddo
   enddo
   deallocate(mesh)
   do k = 1,nfaces
     do j = 1, ne
       do i = 1, ne
         deallocate(gridelem(i,j,k)%nbrs_face)
       enddo
       enddo
   enddo
   deallocate(gridelem)
   deallocate(nbrs_used)
#if 0
   if (outputfiles) then
     close(7)
     close(8)
   endif
#endif
! =======================================
! Generate cube graph...
! =======================================
#if 0
   if (outputfiles) then
     write(9,*) nelem,2*nelem ! metis requires this first line
   endif
#endif
! ============================================
!  Setup the Grid edges (topology independent)
! ============================================
   call initgridedge(gridedge,gridvertex)
! ============================================
!  Setup the Grid edge Indirect addresses
!          (topology dependent)
! ============================================
   nedge = size(gridedge)
   do i = 1,nedge
     call cubesetupedgeindex(gridedge(i))
   enddo
!
   end subroutine cubetopology
! =======================================
! cube_assemble:
!
! Assemble the cube field element by element
! this routine is assumed to be single
! threaded...
! =======================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function cube_assemble(gbl, fld, elem, par, nelemd, nelem, ielem) result(ierr)
!-------------------------------------------------------------------------------
   use element, only: element_t
   use kiapsparallel, only: parallel_t, abortpar, par_any_source, par_any_tag, par_status_size, par_tag, par_double_precision
   use kiapsparmpi, only: par_send, par_recv
!
   real(real_kind) :: gbl(:,:,:,:) ! global output field
   real(real_kind) :: fld(:,:,:) ! local model field
   type(element_t) :: elem ! element to assemble
   type(parallel_t) :: par ! parallel structure
   integer :: nelemd ! number of elements on the node
   integer :: nelem ! number of elements on the node
   integer :: ielem ! local element ctr
   integer :: ierr ! returned error code
! Local variables
   integer :: ie, je, face_no
   integer :: ibase, jbase
   integer :: i, j, k
   integer :: elem_number
   integer :: ne1, ne2 ! element dimensions
   integer :: n1, n2 ! gbl face dimensions
   integer :: nface ! number of faces(must be 6)
   integer :: nlyr ! number of layers
#if defined(_MPI)
   integer :: ectr ! global element counter
   integer tag
   integer :: count ! w/o " :: ", triggers pgi 3.1 f90 bug
   integer pe
   integer status(par_status_size)
   integer mpi_err
#endif
!
   call abortpar(message = 'Because convert_gbl_index is not used cube_assemble is broken. ')
   ne1 = size(fld,1)
   ne2 = size(fld,2)
   nlyr = size(fld,3)
   n1 = size(gbl,1)
   n2 = size(gbl,2)
   nface = size(gbl,3)
! =========================
! Enforce certain rules...
! =========================
   ierr = 0
   if (modulo(n1,ne1)/= 0) then
     ierr =-1
     return
   endif
   if (modulo(n2,ne2)/= 0) then
     ierr =-2
     return
   endif
   if (nface/= 6) then
     ierr =-3
     return
   endif
! =========================================================
! Perform global assembly procedure element by element ...
! =========================================================
   if (par%iproc==par%root) then
     if (ielem<=nelemd) then
       elem_number = elem%vertex%number
       call convert_gbl_index(elem_number,ie,je,face_no)
       if (face_no/= elem%vertex%face_number) call abortpar(message = 'Error in getting face number')
       ibase = ie*ne1
       jbase = je*ne2
       do k = 1,nlyr
         do j = 1, ne2
           do i = 1, ne1
             gbl(i+ibase,j+jbase,face_no,k) = fld(i,j,k)
           enddo
           enddo
       enddo
     endif
#if defined(_MPI)
       if (ielem==nelemd) then
         ectr = nelemd
         do while(ectr<nelem)
           pe = par_any_source
           tag = par_any_tag
           count = ne1*ne2*nlyr
!             call MPI_RECV(fld(1,1,1),  !                  count,       !                  PAR_DOUBLE_PRECISION,   !                  pe,          !                  tag,         &
!                  par%comm,    !                  status,      !                  mpi_err)
           call par_recv(fld(1,1,1),count,pe,ierr = mpi_err)
           elem_number = status(par_tag)
! call convert_gbl_index(elem_number,ie,je,face_no)
           call abortpar(message = 'Because convert_gbl_index is not used for neghbors,the _mpi version needs to be fixed')
           ibase = ie*ne1
           jbase = je*ne2
           do k = 1,nlyr
             do j = 1, ne2
               do i = 1, ne1
                 gbl(i+ibase,j+jbase,face_no,k) = fld(i,j,k)
               enddo
               enddo
           enddo
           ectr = ectr+1
         enddo
       endif
       else
         pe = par%root
         tag = elem%vertex%number
         count = ne1*ne2*nlyr
!       call MPI_SEND(fld(1,1,1),   !            count,        !            PAR_DOUBLE_PRECISION,    !            pe,           !            tag,          !            par%comm,     !            mpi_err)
         call par_send(fld(1,1,1),count,pe,ierr = mpi_err)
#endif
   endif
!
   end function cube_assemble
! ===================================================================
! CubeEdgeCount:
!
!  Determine the number of Grid Edges
!
! ===================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function cubeedgecount() result(nedge)
!-------------------------------------------------------------------------------
   use dimensions, only: ne
   use kiapsparallel, only: abortpar
   implicit none
!
   integer :: nedge
!
   if (0==ne) call abortpar(message = 'Error in cubeedgecount:ne is zero')
   nedge = nfaces*(ne*ne*ninnerelemedge-ncornerelemedge)
!
   end function cubeedgecount
!-------------------------------------------------------------------------------
!
!  SUBROUTINE CubeElemCount
!
!> @brief
!>  - Determines the total number of elements on cubed sphere grid
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function cubeelemcount() result(nelem)
!
!-------------------------------------------------------------------------------
   use dimensions, only: ne
   use kiapsparallel, only: abortpar
!
   implicit none
!
!
   integer :: nelem
!
! Begins
! ------
!
!
   if (ne==0) call abortpar(message = 'Error in cubeelemcount:ne is zero')
!
   nelem = nfaces*ne*ne
!
!
   end function cubeelemcount
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine cubesetupedgeindex(edge)
!-------------------------------------------------------------------------------
   use gridgraph, only: gridedge_t
   use dimensions, only: np
   use control, only: north, south, east, west, neast, seast, swest, nwest
!
   type(gridedge_t), target :: edge
   integer :: np0, sface, dface
   logical :: reverse
   integer, allocatable :: forwardv(:), forwardp(:)
   integer, allocatable :: backwardv(:), backwardp(:)
   integer :: i, ii
!
   ii = edge%tail_face
   np0 = edge%tail%wgtp(ii)
#ifdef TESTGRID
   allocate(forwardp(np0))
   allocate(backwardp(np0))
   do i = 1,np0
     forwardp(i) = i
     backwardp(i) = np0-i+1
   enddo
#endif
   sface = edge%tail_face
   dface = edge%head_face
! Do not reverse the indices
   reverse = .false.
! Under special conditions use index reversal
   if ((sface==south.and.dface==east).or.(sface==east.and.dface==south).or.(sface==north.and.dface==west).or.(sface==west.and.dface==north).or.(sface==south.and.dface==south).or.(sface==north.and.dface==north).or.(sface==east.and.dface==east).or.(sface==west.and.dface==west)) then
     reverse = .true.
     edge%reverse = .true.
   endif
#ifdef TESTGRID
!  Setup the destination indices
   select case(dface)
   case(east)
     edge%headindex%ixv = nv
     edge%headindex%iyv = forwardv
     edge%headindex%ixp = np
     edge%headindex%iyp = forwardp
   case(west)
     edge%headindex%ixv = 1
     edge%headindex%iyv = forwardv
     edge%headindex%ixp = 1
     edge%headindex%iyp = forwardp
   case(north)
     edge%headindex%ixv = forwardv
     edge%headindex%iyv = nv
     edge%headindex%ixp = forwardp
     edge%headindex%iyp = np
   case(south)
     edge%headindex%ixv = forwardv
     edge%headindex%iyv = 1
     edge%headindex%ixp = forwardp
     edge%headindex%iyp = 1
   case(swest)
     edge%headindex%ixv = 1
     edge%headindex%iyv = 1
     edge%headindex%ixp = 1
     edge%headindex%iyp = 1
   case(seast)
     edge%headindex%ixv = nv
     edge%headindex%iyv = 1
     edge%headindex%ixp = np
     edge%headindex%iyp = 1
   case(nwest)
     edge%headindex%ixv = 1
     edge%headindex%iyv = nv
     edge%headindex%ixp = 1
     edge%headindex%iyp = np
   case(neast)
     edge%headindex%ixv = nv
     edge%headindex%iyv = nv
     edge%headindex%ixp = np
     edge%headindex%iyp = np
   case default
     write(*,*) 'SetupEdgeIndex:error in dface select statement'
   end select
! Setup the source indices
   select case(sface)
   case(north)
     edge%tailindex%ixv = forwardv
     if (reverse) edge%tailindex%ixv = backwardv
     edge%tailindex%iyv = nv
     edge%tailindex%ixp = forwardp
     if (reverse) edge%tailindex%ixp = backwardp
     edge%tailindex%iyp = np
   case(south)
     edge%tailindex%ixv = forwardv
     if (reverse) edge%tailindex%ixv = backwardv
     edge%tailindex%iyv = 1
     edge%tailindex%ixp = forwardp
     if (reverse) edge%tailindex%ixp = backwardp
     edge%tailindex%iyp = 1
   case(east)
     edge%tailindex%ixv = nv
     edge%tailindex%iyv = forwardv
     if (reverse) edge%tailindex%iyv = backwardv
     edge%tailindex%ixp = np
     edge%tailindex%iyp = forwardp
     if (reverse) edge%tailindex%iyp = backwardp
   case(west)
     edge%tailindex%ixv = 1
     edge%tailindex%iyv = forwardv
     if (reverse) edge%tailindex%iyv = backwardv
     edge%tailindex%ixp = 1
     edge%tailindex%iyp = forwardp
     if (reverse) edge%tailindex%iyp = backwardp
   case(swest)
     edge%tailindex%ixv = 1
     edge%tailindex%iyv = 1
     edge%tailindex%ixp = 1
     edge%tailindex%iyp = 1
   case(seast)
     edge%tailindex%ixv = nv
     edge%tailindex%iyv = 1
     edge%tailindex%ixp = np
     edge%tailindex%iyp = 1
   case(nwest)
     edge%tailindex%ixv = 1
     edge%tailindex%iyv = nv
     edge%tailindex%ixp = 1
     edge%tailindex%iyp = np
   case(neast)
     edge%tailindex%ixv = nv
     edge%tailindex%iyv = nv
     edge%tailindex%ixp = np
     edge%tailindex%iyp = np
   case default
     write(*,*) 'SetupEdgeIndex:error in sface select statement'
   end select
   deallocate(forwardv)
   deallocate(forwardp)
   deallocate(backwardv)
   deallocate(backwardp)
#endif
!
   end subroutine cubesetupedgeindex
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine getlatticespacing(spherev, dxv, spherep, dxp)
!-------------------------------------------------------------------------------
   use physicalconstants, only: rearth=>kim_rearth
   use dimensions, only: np
!
   type(spherical_polar_t), intent(in   ) :: spherev(np, np)
   real(real_kind) :: dxv
   type(spherical_polar_t), intent(in   ) :: spherep(np, np)
   real(real_kind) :: dxp
   real(real_kind) xcorner, ycorner, zcorner
   real(real_kind) x, y, z
   real(real_kind) chord
   real(real_kind) theta
!
   xcorner = cos(spherev(1,1)%lat)*cos(spherev(1,1)%lon)
   ycorner = cos(spherev(1,1)%lat)*sin(spherev(1,1)%lon)
   zcorner = sin(spherev(1,1)%lat)
   x = cos(spherev(2,1)%lat)*cos(spherev(2,1)%lon)
   y = cos(spherev(2,1)%lat)*sin(spherev(2,1)%lon)
   z = sin(spherev(2,1)%lat)
   chord = sqrt((xcorner-x)**2+(ycorner-y)**2+(zcorner-z)**2)
   theta = 2.0d0*asin(0.50d0*chord)
   dxv = theta*rearth
   x = cos(spherep(1,1)%lat)*cos(spherep(1,1)%lon)
   y = cos(spherep(1,1)%lat)*sin(spherep(1,1)%lon)
   z = sin(spherep(1,1)%lat)
   chord = sqrt((xcorner-x)**2+(ycorner-y)**2+(zcorner-z)**2)
   theta = 2.0d0*asin(0.50d0*chord)
   dxp = theta*rearth
!
   end subroutine getlatticespacing
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module cube
!-------------------------------------------------------------------------------
