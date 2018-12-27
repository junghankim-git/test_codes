!-------------------------------------------------------------------------------
   program testmain
!
   use kinds,           only: i4, l4, r8
   use spline_remap,    only: spline_t, &
                              spline_initialize, spline_finalize,              &
                              spline_set_x,      spline_set,                   &
                              spline_get_fun, spline_get_fun_integral,         &
                              spline_get, spline_interface, unittest_main
   use spline_remap,    only: remap_psm, remap_psm_cyclic
!-------------------------------------------------------------------------------
   implicit none
   include 'mpif.h'
!
   logical(l4) :: reverse = .false.
!
   call unittest_main()
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program
!-------------------------------------------------------------------------------
