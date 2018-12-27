!-------------------------------------------------------------------------------
   program main_program
!-------------------------------------------------------------------------------
!
!  abstract :  (main program) NetCDF file compressor utility
!
!  history log :
!    2018-12-12  junghan kim    initial setup
!
!  structure :
!
!-------------------------------------------------------------------------------
   use kinds,               only: i4, r4, r8, l4
   use lagrangian_particle, only: lagrangian_t, initialize, run, finalize
!
   type(lagrangian_t) :: lag
   integer(i4) :: npart, ivis, uptype, nstep, bndry(2)
   real(r8)    :: dt, xmin, xmax, ymin, ymax
!
   call set_values(npart,ivis,uptype,dt,nstep,xmin,xmax,ymin,ymax,bndry)
   !call initialize(lag,npart,nstep,ivis,uptype,bndry,dt,xmin,xmax,ymin,ymax)
   call initialize(lag,npart,ivis,uptype,dt,nstep,xmin,xmax,ymin,ymax,bndry)
   call run(lag)
   call finalize(lag)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_values(npart,ivis,uptype,dt,nstep,xmin,xmax,ymin,ymax,bndry)
   implicit none
   integer(i4), intent(inout) :: npart, nstep, ivis, uptype, bndry(2)
   real(r8)   , intent(inout) :: xmin, xmax, ymin, ymax, dt
! local variables
   namelist/inputs/npart,ivis,uptype,dt,nstep,xmin,xmax,ymin,ymax,bndry
!
   npart  = 50
   ivis   = 0
   uptype = 0
   dt     = 0.001_r8
   nstep  = 10000
   xmin   = -50.0_r8
   xmax   =  50.0_r8
   ymin   = -50.0_r8
   ymax   =  50.0_r8
   bndry  = (/1,2/)

!
   open(31,file='inputs.nl',status='old')
   read(31,inputs)
   close(31)
!
   return
   end subroutine set_values
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program main_program
!-------------------------------------------------------------------------------
