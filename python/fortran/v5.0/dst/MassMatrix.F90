#include "KIM.h"
!-------------------------------------------------------------------------------
   module massmatrix
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
   use kiapsbase, only: real_kind=>kim_real8_kind
   use dimensions, only: np, nelemd
   use quadrature, only: quadrature_t, gauss, gausslobatto
   use element, only: element_t
   use kiapsparallel, only: parallel_t
   use edge, only: edgebuffer_t, edgevpack, edgevunpack, &
                                                freeedgebuffer, initedgebuffer
   use bndry, only: bndry_exchangev
   implicit none
!
   private
   public :: mass_matrix
!
   contains
! ===========================================
! mass_matrix:
!
! Compute the mass matrix for each element...
! ===========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine mass_matrix(par, elem)
!
   type(parallel_t), intent(in   ) :: par
   type(element_t) :: elem(:)
   type(edgebuffer_t) :: edge
   real(real_kind) da ! area element
   type(quadrature_t) :: gp
   integer ii
   integer i, j
   integer kptr
   integer iptr
! ===================
! begin code
! ===================
!
   call initedgebuffer(edge,1)
! =================================================
! mass matrix on the velocity grid
! =================================================
   gp = gausslobatto(np)
   do ii = 1,nelemd
     do j = 1, np
       do i = 1, np
! MNL: metric term for map to reference element is now in metdet!
         elem(ii)%mp(i,j) = gp%weights(i)*gp%weights(j)
         elem(ii)%rmp(i,j) = elem(ii)%mp(i,j)
       enddo
       enddo
         kptr = 0
         call edgevpack(edge,elem(ii)%rmp,1,kptr,elem(ii)%desc)
   enddo
! ==============================
! Insert boundary exchange here
! ==============================
   call bndry_exchangev(par,edge)
   do ii = 1,nelemd
     kptr = 0
     call edgevunpack(edge,elem(ii)%rmp,1,kptr,elem(ii)%desc)
     do j = 1,np
       do i = 1, np
         elem(ii)%rmp(i,j) = 1.0d0/elem(ii)%rmp(i,j)
       enddo
     enddo
   enddo
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
   deallocate(gp%points)
   deallocate(gp%weights)
! =============================================
! compute spherical element mass matrix
! =============================================
   do ii = 1,nelemd
     do j = 1, np
       do i = 1, np
         elem(ii)%spheremp(i,j) = elem(ii)%mp(i,j)*elem(ii)%metdet(i,j)
         elem(ii)%rspheremp(i,j) = elem(ii)%spheremp(i,j)
       enddo
       enddo
         kptr = 0
         call edgevpack(edge,elem(ii)%rspheremp,1,kptr,elem(ii)%desc)
   enddo
   call bndry_exchangev(par,edge)
   do ii = 1,nelemd
     kptr = 0
     call edgevunpack(edge,elem(ii)%rspheremp,1,kptr,elem(ii)%desc)
     do j = 1,np
       do i = 1, np
         elem(ii)%rspheremp(i,j) = 1.0d0/elem(ii)%rspheremp(i,j)
       enddo
     enddo
   enddo
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
! =============================================
! compute the mass matrix
! =============================================
! Jose Garcia: Not sure but I think this code is just dead code
!do ii=1,nelemd
!   iptr=1
!   do j=1,np
!      do i=1,np
!         elem(ii)%mp(i,j)=elem(ii)%mp(i,j)
!         iptr=iptr+1
!      end do
!   end do
!end do
   call freeedgebuffer(edge)
!
   end subroutine mass_matrix
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module massmatrix
!-------------------------------------------------------------------------------
