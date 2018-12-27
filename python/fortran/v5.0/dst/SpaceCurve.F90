#include "KIM.h"
!-------------------------------------------------------------------------------
   module spacecurve
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
   use kiapsbase, only: iulog=>kim_iu_log
   implicit none
!
   private
!
   type, public :: factor_t
     integer :: numfact
     integer, dimension(:), pointer :: factors
   end type factor_t
   integer, public, dimension(:,:), allocatable :: ordered
   integer, public, dimension(:,:), allocatable :: dir ! direction to move along each level
   integer, public, dimension(:), allocatable :: pos ! position along each of the axes
   integer, public :: maxdim ! dimensionality of entire space
   integer, public :: vcnt ! visitation count
   logical, private :: verbose = .false.
   type(factor_t), public :: fact
!JMD new addition
!
   save :: fact
   public :: map
   public :: hilbert_old
   public :: peanom, hilbert, cinco
   public :: gencurve
   public :: genspacecurve
   public :: log2, factor
   public :: printcurve
   public :: isfactorable, isloadbalanced
   public :: genspacepart
!
   contains
!---------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   recursive function cinco(l, type, ma, md, ja, jd) result(ierr)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: l, type, ma, md, ja, jd
   integer :: lma, lmd, lja, ljd, ltype
   integer :: ll
   integer :: ierr
   logical :: debug = .false.
!
   ll = l
   if (ll.gt.1) ltype = fact%factors(ll-1) ! set the next type of space curve
!--------------------------------------------------------------
!  Position [0,0]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,21) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'Cinco:after position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [1,0]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,22) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [1,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [2,0]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,23) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [2,1]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,24) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [2,2]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd = md
   lja = ma
   ljd =-md
   if (ll.gt.1) then
     if (debug) write(iulog,25) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [1,2]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd =-md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,26) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [1,1]
!--------------------------------------------------------------
   lma = ma
   lmd =-md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,27) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [0,1]
!--------------------------------------------------------------
   lma = ma
   lmd =-md
   lja = mod(ma+1,maxdim)
   ljd = md
   if (ll.gt.1) then
     if (debug) write(iulog,28) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [0,2]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,29) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [0,3]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,30) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [0,4]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,31) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [1,4]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = mod(ma+1,maxdim)
   ljd =-md
   if (ll.gt.1) then
     if (debug) write(iulog,32) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [1,3]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd =-md
   lja = ma
   ljd = md
   if (ll.gt.1) then
     if (debug) write(iulog,33) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [2,3]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,34) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [2,4]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,35) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [3,4]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,36) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [4,4]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = mod(ma+1,maxdim)
   ljd =-md
   if (ll.gt.1) then
     if (debug) write(iulog,37) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [4,3]
!--------------------------------------------------------------
   lma = ma
   lmd =-md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,38) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [3,3]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd =-md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,39) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [3,2]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd =-md
   lja = ma
   ljd = md
   if (ll.gt.1) then
     if (debug) write(iulog,40) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [4,2]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = mod(ma+1,maxdim)
   ljd =-md
   if (ll.gt.1) then
     if (debug) write(iulog,41) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [4,1]
!--------------------------------------------------------------
   lma = ma
   lmd =-md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,42) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [3,1]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd =-md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,43) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [3,0]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd =-md
   lja = ma
   ljd = md
   if (ll.gt.1) then
     if (debug) write(iulog,44) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
!  Position [4,0]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = ja
   ljd = jd
   if (ll.gt.1) then
     if (debug) write(iulog,45) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
   21 format('Call cinco pos [0,0] level ',i1,' at(',i2,',',i2,') ',4(i3))
   22 format('Call cinco pos [1,0] level ',i1,' at(',i2,',',i2,') ',4(i3))
   23 format('Call cinco pos [2,0] level ',i1,' at(',i2,',',i2,') ',4(i3))
   24 format('Call cinco pos [2,1] level ',i1,' at(',i2,',',i2,') ',4(i3))
   25 format('Call cinco pos [2,2] level ',i1,' at(',i2,',',i2,') ',4(i3))
   26 format('Call cinco pos [1,2] level ',i1,' at(',i2,',',i2,') ',4(i3))
   27 format('Call cinco pos [1,1] level ',i1,' at(',i2,',',i2,') ',4(i3))
   28 format('Call cinco pos [0,1] level ',i1,' at(',i2,',',i2,') ',4(i3))
   29 format('Call cinco pos [0,2] level ',i1,' at(',i2,',',i2,') ',4(i3))
   30 format('Call cinco pos [0,3] level ',i1,' at(',i2,',',i2,') ',4(i3))
   31 format('Call cinco pos [0,4] level ',i1,' at(',i2,',',i2,') ',4(i3))
   32 format('Call cinco pos [1,4] level ',i1,' at(',i2,',',i2,') ',4(i3))
   33 format('Call cinco pos [1,3] level ',i1,' at(',i2,',',i2,') ',4(i3))
   34 format('Call cinco pos [2,3] level ',i1,' at(',i2,',',i2,') ',4(i3))
   35 format('Call cinco pos [2,4] level ',i1,' at(',i2,',',i2,') ',4(i3))
   36 format('Call cinco pos [3,4] level ',i1,' at(',i2,',',i2,') ',4(i3))
   37 format('Call cinco pos [4,4] level ',i1,' at(',i2,',',i2,') ',4(i3))
   38 format('Call cinco pos [4,3] level ',i1,' at(',i2,',',i2,') ',4(i3))
   39 format('Call cinco pos [3,3] level ',i1,' at(',i2,',',i2,') ',4(i3))
   40 format('Call cinco pos [3,2] level ',i1,' at(',i2,',',i2,') ',4(i3))
   41 format('Call cinco pos [4,2] level ',i1,' at(',i2,',',i2,') ',4(i3))
   42 format('Call cinco pos [4,1] level ',i1,' at(',i2,',',i2,') ',4(i3))
   43 format('Call cinco pos [3,1] level ',i1,' at(',i2,',',i2,') ',4(i3))
   44 format('Call cinco pos [3,0] level ',i1,' at(',i2,',',i2,') ',4(i3))
   45 format('Call cinco pos [4,0] level ',i1,' at(',i2,',',i2,') ',4(i3))
!
   end function cinco
!---------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   recursive function peanom(l, type, ma, md, ja, jd) result(ierr)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: l, type, ma, md, ja, jd
   integer :: lma, lmd, lja, ljd, ltype
   integer :: ll
   integer :: ierr
   logical :: debug = .false.
!
   ll = l
   if (ll.gt.1) ltype = fact%factors(ll-1) ! set the next type of space curve
!--------------------------------------------------------------
!  Position [0,0]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,21) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
! Position [0,1]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,22) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,1] ',pos
   endif
!--------------------------------------------------------------
! Position [0,2]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,23) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,2] ',pos
   endif
!--------------------------------------------------------------
! Position [1,2]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,24) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [1,2] ',pos
   endif
!--------------------------------------------------------------
! Position [2,2]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = mod(lma+1,maxdim)
   ljd =-lmd
   if (ll.gt.1) then
     if (debug) write(iulog,25) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [2,2] ',pos
   endif
!--------------------------------------------------------------
! Position [2,1]
!--------------------------------------------------------------
   lma = ma
   lmd =-md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,26) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [2,1] ',pos
   endif
!--------------------------------------------------------------
! Position [1,1]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd =-md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,27) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [1,1] ',pos
   endif
!--------------------------------------------------------------
! Position [1,0]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd =-md
   lja = mod(lma+1,maxdim)
   ljd =-lmd
   if (ll.gt.1) then
     if (debug) write(iulog,28) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [1,0] ',pos
   endif
!--------------------------------------------------------------
! Position [2,0]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = ja
   ljd = jd
   if (ll.gt.1) then
     if (debug) write(iulog,29) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [2,0] ',pos
   endif
   21 format('Call peanom pos [0,0] level ',i1,' at(',i2,',',i2,') ',4(i3))
   22 format('Call peanom pos [0,1] level ',i1,' at(',i2,',',i2,') ',4(i3))
   23 format('Call peanom pos [0,2] level ',i1,' at(',i2,',',i2,') ',4(i3))
   24 format('Call peanom pos [1,2] level ',i1,' at(',i2,',',i2,') ',4(i3))
   25 format('Call peanom pos [2,2] level ',i1,' at(',i2,',',i2,') ',4(i3))
   26 format('Call peanom pos [2,1] level ',i1,' at(',i2,',',i2,') ',4(i3))
   27 format('Call peanom pos [1,1] level ',i1,' at(',i2,',',i2,') ',4(i3))
   28 format('Call peanom pos [1,0] level ',i1,' at(',i2,',',i2,') ',4(i3))
   29 format('Call peanom pos [2,0] level ',i1,' at(',i2,',',i2,') ',4(i3))
!
   end function peanom
!---------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   recursive function hilbert(l, type, ma, md, ja, jd) result(ierr)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: l, type, ma, md, ja, jd
   integer :: lma, lmd, lja, ljd, ltype
   integer :: ll
   integer :: ierr
   logical :: debug = .false.
!
   ll = l
   if (ll.gt.1) ltype = fact%factors(ll-1) ! set the next type of space curve
!--------------------------------------------------------------
!  Position [0,0]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,21) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,0] ',pos
   endif
!--------------------------------------------------------------
! Position [0,1]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = lma
   ljd = lmd
   if (ll.gt.1) then
     if (debug) write(iulog,22) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [0,1] ',pos
   endif
!--------------------------------------------------------------
! Position [1,1]
!--------------------------------------------------------------
   lma = ma
   lmd = md
   lja = mod(ma+1,maxdim)
   ljd =-md
   if (ll.gt.1) then
     if (debug) write(iulog,23) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [1,1] ',pos
   endif
!--------------------------------------------------------------
! Position [1,0]
!--------------------------------------------------------------
   lma = mod(ma+1,maxdim)
   lmd =-md
   lja = ja
   ljd = jd
   if (ll.gt.1) then
     if (debug) write(iulog,24) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
     ierr = gencurve(ll-1,ltype,lma,lmd,lja,ljd)
     if (debug) call printcurve(ordered)
   else
     ierr = incrementcurve(lja,ljd)
     if (debug) write(iulog,*) 'After position [1,0] ',pos
   endif
   21 format('Call hilbert pos [0,0] level ',i1,' at(',i2,',',i2,') ',4(i3))
   22 format('Call hilbert pos [0,1] level ',i1,' at(',i2,',',i2,') ',4(i3))
   23 format('Call hilbert pos [1,1] level ',i1,' at(',i2,',',i2,') ',4(i3))
   24 format('Call hilbert pos [1,0] level ',i1,' at(',i2,',',i2,') ',4(i3))
!
   end function hilbert
!---------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function incrementcurve(ja, jd) result(ierr)
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: ja, jd
   integer :: ierr
!
   ordered(pos(0)+1,pos(1)+1) = vcnt
   vcnt = vcnt+1
   pos(ja) = pos(ja)+jd
   ierr = 0
!
   end function incrementcurve
!---------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   recursive function hilbert_old(l, d, ma, md, ja, jd) result(ierr)
!
   integer :: l, d ! log base 2 of levels and dimensions left
   integer :: ma, md ! main axis and direction
   integer :: ja, jd ! joiner axis and direction
   integer :: ierr
   integer :: axis
   integer :: ll
!
   if (verbose) write(iulog,10) l,d,ma,md,ja,jd,pos(0),pos(1)
   ll = l ! copy this to a temporary variable
   if (d==0) then
     ll = ll-1
     if (ll==0) then
       return
     endif
     axis = ja
     if (dir(ll,axis)/= jd) then ! do not move away from joiner plane
       axis = mod(axis+1,maxdim) ! next axis
     endif
     if (verbose) write(iulog,*) 'hilbert_old:call hilbert_old(l,d) #1:'
     ierr = hilbert_old(ll,maxdim,axis,dir(ll,axis),ja,jd)
     dir(ll,ja) =-dir(ll,ja)
     return
   endif
   axis = mod(ma+1,maxdim)
   if (verbose) write(iulog,*) 'hilbert_old:before call hilbert_old(l,d) #2:'
   ierr = hilbert_old(ll,d-1,axis,dir(ll,axis),ma,md)
   if (verbose) write(iulog,*) 'hilbert_old:after call hilbert_old(l,d) #2:'
   if (verbose) write(iulog,30) l,d,ma,md,ja,jd,pos(0),pos(1)
   pos(ma) = pos(ma)+md
   dir(ll,ma) =-dir(ll,ma)
!----------------------------------
!  Mark this node as visited
!----------------------------------
   if (verbose) write(iulog,20) l,d,ma,md,ja,jd,pos(0),pos(1)
   vcnt = vcnt+1
   if (verbose) write(iulog,15) pos(0)+1,pos(1)+1,vcnt
   if (verbose) write(iulog,*) ' '
   if (verbose) write(iulog,*) ' '
   ordered(pos(0)+1,pos(1)+1) = vcnt
   if (verbose) write(iulog,*) 'hilbert_old:before call hilbert_old(l,d) #3:'
   ierr = hilbert_old(ll,d-1,axis,dir(ll,axis),ja,jd)
   if (verbose) write(iulog,*) 'hilbert_old:after call hilbert_old(l,d) #3:'
   10 format('hilbert_old:entering hilbert_old(l,d,ma,md,ja,jd) are:',&
                2(i4),' [',2(i3),'][',2(i3),']',2(i3))
   15 format('hilbert_old:mark element {x,y,ordered}:',3(i4))
   20 format('hilbert_old:before visit code(l,d,ma,md,ja,jd) are:',&
                2(i4),' [',2(i3),'][',2(i3),']',2(i3))
   30 format('hilbert_old:after call hilbert_old(l,d) #2 :(l,d,ma,md,ja,jd are:',&
                2(i4),' [',2(i3),'][',2(i3),']',2(i3))
!
   end function hilbert_old
!---------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function log2(n)
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: n
   integer :: log2, tmp
!
!  Find the log2 of input value
!
!
   log2 = 1
   tmp = n
   do while(tmp/2.ne.1)
     tmp = tmp/2
     log2 = log2+1
   enddo
!
   end function log2
!---------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function isloadbalanced(nelem, npart)
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: nelem, npart
   logical :: isloadbalanced
   integer :: tmp1
!
   tmp1 = nelem/npart
   if (npart*tmp1==nelem) then
     isloadbalanced = .true.
   else
     isloadbalanced = .false.
   endif
!
   end function isloadbalanced
!---------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function gencurve(l, type, ma, md, ja, jd) result(ierr)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: l, type, ma, md, ja, jd
   integer :: ierr
!
   if (type==2) then
     ierr = hilbert(l,type,ma,md,ja,jd)
     elseif (type==3) then
     ierr = peanom(l,type,ma,md,ja,jd)
     elseif (type==5) then
     ierr = cinco(l,type,ma,md,ja,jd)
   endif
!
   end function gencurve
!---------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function factor(num) result(res)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: num
   type(factor_t) :: res
   integer :: tmp, tmp2, tmp3, tmp5
   integer :: i, n
   logical :: found
! --------------------------------------
! Allocate for max # of factors
! --------------------------------------
!
   tmp = num
   tmp2 = log2(num)
   allocate(res%factors(tmp2))
   n = 0
!-----------------------
!  Look for factors of 2
!-----------------------
   found = .true.
   do while(found)
     found = .false.
     tmp2 = tmp/2
     if (tmp2*2==tmp) then
       n = n+1
       res%factors(n) = 2
       found = .true.
       tmp = tmp2
     endif
   enddo
!-----------------------
!  Look for factors of 3
!-----------------------
   found = .true.
   do while(found)
     found = .false.
     tmp3 = tmp/3
     if (tmp3*3==tmp) then
       n = n+1
       res%factors(n) = 3
       found = .true.
       tmp = tmp3
     endif
   enddo
!-----------------------
!  Look for factors of 5
!-----------------------
   found = .true.
   do while(found)
     found = .false.
     tmp5 = tmp/5
     if (tmp5*5==tmp) then
       n = n+1
       res%factors(n) = 5
       found = .true.
       tmp = tmp5
     endif
   enddo
   tmp = 1
   do i = 1,n
     tmp = tmp*res%factors(i)
   enddo
   if (tmp==num) then
     res%numfact = n
   else
     res%numfact =-1
   endif
!
   end function factor
!---------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function isfactorable(n)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: n
   type(factor_t) :: fact
   logical :: isfactorable
!
   fact = factor(n)
   if (fact%numfact.ne.-1) then
     isfactorable = .true.
   else
     isfactorable = .false.
   endif
!
   end function isfactorable
!------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine map(l)
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: l, d
   integer :: type, ierr
!
   d = size(pos)
   pos = 0
   maxdim=d
   vcnt = 0
   type = fact%factors(l)
!
   ierr = gencurve(l,type,0,1,0,1)
!
   end subroutine map
!---------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine genspacecurve(mesh)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, target, intent(inout) :: mesh(:,:)
   integer :: level, dim
   integer :: gridsize
!  Setup the size of the grid to traverse
!
   dim=2
   gridsize = size(mesh,dim=1)
   fact = factor(gridsize)
   level = fact%numfact
   if (verbose) write(iulog,*) 'GenSpacecurve:level is ',level
   allocate(ordered(gridsize,gridsize))
! Setup the working arrays for the traversal
   allocate(pos(0:dim-1))
!  The array ordered will contain the visitation order
   ordered(:,:) = 0
   call map(level)
   mesh(:,:) = ordered(:,:)
!
   end subroutine genspacecurve
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine printcurve(mesh)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, target :: mesh(:,:)
   integer :: gridsize, i
!
   gridsize = size(mesh,dim=1)
   if (gridsize==2) then
     write(iulog,*) "a level 1 hilbert curve:"
     write(iulog,*) "------------------------"
     do i = 1,gridsize
       write(iulog,2) mesh(1,i),mesh(2,i)
     enddo
   elseif (gridsize==3) then
     write(iulog,*) "a level 1 peano meandering curve:"
     write(iulog,*) "---------------------------------"
     do i = 1,gridsize
       write(iulog,3) mesh(1,i),mesh(2,i),mesh(3,i)
     enddo
   elseif (gridsize==4) then
     write(iulog,*) "a level 2 hilbert curve:"
     write(iulog,*) "------------------------"
     do i = 1,gridsize
       write(iulog,4) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i)
     enddo
   elseif (gridsize==5) then
     write(iulog,*) "a level 1 cinco curve:"
     write(iulog,*) "------------------------"
     do i = 1,gridsize
       write(iulog,5) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),mesh(5,i)
     enddo
   elseif (gridsize==6) then
     write(iulog,*) "a level 1 hilbert and level 1 peano curve:"
     write(iulog,*) "------------------------------------------"
     do i = 1,gridsize
       write(iulog,6) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),mesh(5,i),mesh(6,i)
     enddo
   elseif (gridsize==8) then
     write(iulog,*) "a level 3 hilbert curve:"
     write(iulog,*) "------------------------"
     do i = 1,gridsize
       write(iulog,8) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),&
              mesh(5,i),mesh(6,i),mesh(7,i),mesh(8,i)
     enddo
   elseif (gridsize==9) then
     write(iulog,*) "a level 2 peano meandering curve:"
     write(iulog,*) "---------------------------------"
     do i = 1,gridsize
       write(iulog,9) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),&
            mesh(5,i),mesh(6,i),mesh(7,i),mesh(8,i),&
                                                                 mesh(9,i)
     enddo
   elseif (gridsize==10) then
     write(iulog,*) "a level 1 hilbert and level 1 cinco curve:"
     write(iulog,*) "---------------------------------"
     do i = 1,gridsize
       write(iulog,10) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),&
            mesh(5,i),mesh(6,i),mesh(7,i),mesh(8,i),&
                                               mesh(9,i),mesh(10,i)
     enddo
   elseif (gridsize==12) then
     write(iulog,*) "a level 2 hilbert and level 1 peano curve:"
     write(iulog,*) "------------------------------------------"
     do i = 1,gridsize
       write(iulog,12) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),&
            mesh(5,i),mesh(6,i),mesh(7,i),mesh(8,i),&
           mesh(9,i),mesh(10,i),mesh(11,i),mesh(12,i)
     enddo
   elseif (gridsize==15) then
     write(iulog,*) "a level 1 peano and level 1 cinco curve:"
     write(iulog,*) "------------------------"
     do i = 1,gridsize
       write(iulog,15) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),&
            mesh(5,i),mesh(6,i),mesh(7,i),mesh(8,i),&
         mesh(9,i),mesh(10,i),mesh(11,i),mesh(12,i),&
                            mesh(13,i),mesh(14,i),mesh(15,i)
     enddo
   elseif (gridsize==16) then
     write(iulog,*) "a level 4 hilbert curve:"
     write(iulog,*) "------------------------"
     do i = 1,gridsize
       write(iulog,16) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),&
            mesh(5,i),mesh(6,i),mesh(7,i),mesh(8,i),&
         mesh(9,i),mesh(10,i),mesh(11,i),mesh(12,i),&
          mesh(13,i),mesh(14,i),mesh(15,i),mesh(16,i)
     enddo
   elseif (gridsize==18) then
     write(iulog,*) "a level 1 hilbert and level 2 peano curve:"
     write(iulog,*) "------------------------------------------"
     do i = 1,gridsize
       write(iulog,18) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),&
            mesh(5,i),mesh(6,i),mesh(7,i),mesh(8,i),&
         mesh(9,i),mesh(10,i),mesh(11,i),mesh(12,i),&
        mesh(13,i),mesh(14,i),mesh(15,i),mesh(16,i),&
                                              mesh(17,i),mesh(18,i)
     enddo
   elseif (gridsize==20) then
     write(iulog,*) "a level 2 hilbert and level 1 cinco curve:"
     write(iulog,*) "------------------------------------------"
     do i = 1,gridsize
       write(iulog,20) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),&
            mesh(5,i),mesh(6,i),mesh(7,i),mesh(8,i),&
         mesh(9,i),mesh(10,i),mesh(11,i),mesh(12,i),&
        mesh(13,i),mesh(14,i),mesh(15,i),mesh(16,i),&
          mesh(17,i),mesh(18,i),mesh(19,i),mesh(20,i)
     enddo
   elseif (gridsize==24) then
     write(iulog,*) "a level 3 hilbert and level 1 peano curve:"
     write(iulog,*) "------------------------------------------"
     do i = 1,gridsize
       write(iulog,24) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),&
            mesh(5,i),mesh(6,i),mesh(7,i),mesh(8,i),&
         mesh(9,i),mesh(10,i),mesh(11,i),mesh(12,i),&
        mesh(13,i),mesh(14,i),mesh(15,i),mesh(16,i),&
        mesh(17,i),mesh(18,i),mesh(19,i),mesh(20,i),&
          mesh(21,i),mesh(22,i),mesh(23,i),mesh(24,i)
     enddo
   elseif (gridsize==25) then
     write(iulog,*) "a level 2 cinco curve:"
     write(iulog,*) "------------------------------------------"
     do i = 1,gridsize
       write(iulog,25) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),&
            mesh(5,i),mesh(6,i),mesh(7,i),mesh(8,i),&
         mesh(9,i),mesh(10,i),mesh(11,i),mesh(12,i),&
        mesh(13,i),mesh(14,i),mesh(15,i),mesh(16,i),&
        mesh(17,i),mesh(18,i),mesh(19,i),mesh(20,i),&
        mesh(21,i),mesh(22,i),mesh(23,i),mesh(24,i),&
                                                                mesh(25,i)
     enddo
   elseif (gridsize==27) then
     write(iulog,*) "a level 3 peano meandering curve:"
     write(iulog,*) "---------------------------------"
     do i = 1,gridsize
       write(iulog,27) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),&
            mesh(5,i),mesh(6,i),mesh(7,i),mesh(8,i),&
         mesh(9,i),mesh(10,i),mesh(11,i),mesh(12,i),&
        mesh(13,i),mesh(14,i),mesh(15,i),mesh(16,i),&
        mesh(17,i),mesh(18,i),mesh(19,i),mesh(20,i),&
        mesh(21,i),mesh(22,i),mesh(23,i),mesh(24,i),&
                            mesh(25,i),mesh(26,i),mesh(27,i)
     enddo
   elseif (gridsize==32) then
     write(iulog,*) "a level 5 hilbert curve:"
     write(iulog,*) "------------------------"
     do i = 1,gridsize
       write(iulog,32) mesh(1,i),mesh(2,i),mesh(3,i),mesh(4,i),&
            mesh(5,i),mesh(6,i),mesh(7,i),mesh(8,i),&
         mesh(9,i),mesh(10,i),mesh(11,i),mesh(12,i),&
        mesh(13,i),mesh(14,i),mesh(15,i),mesh(16,i),&
        mesh(17,i),mesh(18,i),mesh(19,i),mesh(20,i),&
        mesh(21,i),mesh(22,i),mesh(23,i),mesh(24,i),&
        mesh(25,i),mesh(26,i),mesh(27,i),mesh(28,i),&
          mesh(29,i),mesh(30,i),mesh(31,i),mesh(32,i)
     enddo
   endif
   2 format('|',2(i2,'|'))
   3 format('|',3(i2,'|'))
   4 format('|',4(i2,'|'))
   5 format('|',5(i2,'|'))
   6 format('|',6(i2,'|'))
   8 format('|',8(i2,'|'))
   9 format('|',9(i2,'|'))
   10 format('|',10(i2,'|'))
   12 format('|',12(i3,'|'))
   15 format('|',15(i3,'|'))
   16 format('|',16(i3,'|'))
   18 format('|',18(i3,'|'))
   20 format('|',20(i3,'|'))
   24 format('|',24(i3,'|'))
   25 format('|',25(i3,'|'))
   27 format('|',27(i3,'|'))
   32 format('|',32(i4,'|'))
!
   end subroutine printcurve
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine genspacepart(gridedge, gridvertex)
!-------------------------------------------------------------------------------
   use dimensions, only: npart
   use gridgraph, only: gridedge_t, gridvertex_t
   implicit none
!
   type(gridvertex_t), intent(inout) :: gridvertex(:)
   type(gridedge_t), intent(inout) :: gridedge(:)
   integer :: nelem, nelem_edge, nelemd
   integer :: head_part, tail_part
   integer :: j, k, tmp1, id, s1, extra
!
   nelem = size(gridvertex(:))
   nelem_edge = size(gridedge(:))
   nelemd = nelem/npart
! every cpu gets nelemd elements,but the first 'extra' get nelemd+1
   extra = mod(nelem,npart)
   s1 = extra*(nelemd+1)
! split curve into two curves:
! 1 ... s1  s2 ... nelem
!
!  s1 = extra*(nelemd+1)         (count be 0)
!  s2 = s1+1
!
! First region gets nelemd+1 elements per Processor
! Second region gets nelemd elements per Processor
! ===========================================
!  Add the partitioning information into the
!    Grid Vertex and Grid Edge structures
! ===========================================
   do k = 1,nelem
     id = gridvertex(k)%spacecurve
     if (id<=s1) then
       tmp1 = id/(nelemd+1)
       gridvertex(k)%processor_number = tmp1+1
     else
       id = id-s1
       tmp1 = id/nelemd
       gridvertex(k)%processor_number = extra+tmp1+1
     endif
   enddo
#if 0
   write(iulog,*) 'Space-filling curve parititioning:'
   do k = 1,nelem
     write(iulog,*) k,gridvertex(k)%processor_number
   enddo
   stop 'halting:at the end of genspacepart'
#endif
!
   end subroutine genspacepart
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module spacecurve
!-------------------------------------------------------------------------------
