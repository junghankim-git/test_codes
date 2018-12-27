!-------------------------------------------------------------------------------
   program test
!
   use kinds,       only: i4, l4, r4, r8
   use mesh_contrl, only: axis_x, axis_y, axis_xy, axis_mxy, axis_0
   use mesh_contrl, only: mesh_convert, mesh_dir_convert, mesh_print
   use space_filling_curve, only: sfc_t, sfc_initializei, sfc_finalize,        &
                                  sfc_gen_curve, sfc_get_sfc
!
   type(sfc_t) :: sfc
   integer(i4), parameter :: n = 2
   integer(i4), dimension(n, n) :: mesh, cmesh
   integer(i4) :: ne = 4
   integer(i4), dimension(:), allocatable :: factors
   logical(l4) :: isfac
   integer(i4), dimension(:,:), allocatable :: sfc_curve
!
   mesh(1,1) = 1
   mesh(2,1) = 2
   mesh(1,2) = 3
   mesh(2,2) = 4
!
   print*,'original'
   call mesh_print(mesh)
!
   print*,'x-axis'
   call mesh_convert(axis_x,mesh,cmesh)
   call mesh_print(cmesh)
!
   print*,'y-axis'
   call mesh_convert(axis_y,mesh,cmesh)
   call mesh_print(cmesh)
!
   print*,'xy-axis'
   call mesh_convert(axis_xy,mesh,cmesh)
   call mesh_print(cmesh)
!
   print*,'mxy-axis'
   call mesh_convert(axis_mxy,mesh,cmesh)
   call mesh_print(cmesh)
!
   print*,'0-axis'
   call mesh_convert(axis_0,mesh,cmesh)
   call mesh_print(cmesh)
!
   print*,' '
   print*,' '
!
   print*,'original'
   call mesh_print(mesh,.true.)
!
   print*,'x-axis'
   call mesh_dir_convert(axis_x,mesh,cmesh)
   call mesh_print(cmesh,.true.)
!
   print*,'y-axis'
   call mesh_dir_convert(axis_y,mesh,cmesh)
   call mesh_print(cmesh,.true.)
!
   print*,'xy-axis'
   call mesh_dir_convert(axis_xy,mesh,cmesh)
   call mesh_print(cmesh,.true.)
!
   print*,'mxy-axis'
   call mesh_dir_convert(axis_mxy,mesh,cmesh)
   call mesh_print(cmesh,.true.)
!
   print*,'0-axis'
   call mesh_dir_convert(axis_0,mesh,cmesh)
   call mesh_print(cmesh,.true.)
!===============================
   print*,' '
   print*,' '
!
   print*,'Space filling curve'
   call sfc_initialize(sfc,ne)
   print*,sfc%factors
!
   print*,' '
   call sfc_gen_curve(sfc)
   call mesh_print(sfc%curve%joiner,.true.)
   print*,'A'
   call sfc_get_sfc(sfc,sfc_curve)
   print*,'B'
   call mesh_print(sfc_curve)
   print*,'C'
!
   call sfc_finalize(sfc)
   deallocate(sfc_curve)
!
   end program test
!-------------------------------------------------------------------------------
