!-------------------------------------------------------------------------------
   program test
!
   use kinds,        only: i4, l4, r4, r8
   use cubed_sphere, only: cubed_sphere_t, cs_initialize, cs_finalize, cs_run
   use cubed_sphere, only: cs_print_coordinates, cs_print_directions,         &
                           cs_write_cube_infos, cs_write_lonslats, cs_write_scrip, &
                           cs_get_distance, cs_find_index_point
   implicit none
!
   type(cubed_sphere_t) :: cs
!
   integer(i4) :: np, ne, nprocs, ie, je
   logical(l4) :: isrotated
!   real(r8), dimension(:,:,:), allocatable :: dist
!
   np = 4
   ne = 30
   nprocs = 6
!
   call readnl()
!
   call cs_initialize(cs,np,ne,nprocs,isrotated)
!
   call cs_run(cs)
  !call cs_print_coordinates(cs)
  !call cs_print_directions(cs)
!
   call cs_write_cube_infos(cs,'cs.nc')
   call cs_write_lonslats(cs,'lonslats.nc')
   call cs_write_scrip(cs)
!
   call cs_get_distance(cs)
                              ! lon    lat
   call cs_find_index_point(cs,- 75.0d0,- 30.0d0)
   call cs_find_index_point(cs,+160.0d0,- 30.0d0)
   call cs_find_index_point(cs,-120.0d0,+ 30.0d0)
   call cs_find_index_point(cs,+180.0d0,+ 00.0d0)
   call cs_find_index_point(cs,-119.0d0,+ 30.0d0)
   call cs_find_index_point(cs,- 90.0d0,- 30.0d0)
   call cs_find_index_point(cs,-150.0d0,  30.0d0)
   call cs_find_index_point(cs,-120.0d0,  30.0d0)
   call cs_find_index_point(cs,-125.0d0,- 12.0d0)
!
   call cs_finalize(cs)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine readnl()
   implicit none
   namelist/inputs/np, ne, nprocs, isrotated
!
   open(31,file='inputs.nl',status='old')
   read(31,inputs)
   close(31)
!
   print*,'* np        = ',np
   print*,'* ne        = ',ne
   print*,'* nprocs    = ',nprocs
   print*,'* isrotated = ',isrotated
!
   end subroutine readnl
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program test
