!-------------------------------------------------------------------------------
   program genscrip_input
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
   use kinds      , only: i4, r4, r8
   use scrip_input, only: scrip_t, scrip_initialize_read, scrip_finalize, scrip_check
   use scrip_khkim, only: khkim_t, khkim_initialize, khkim_convert, khkim_finalize
   use scrip_hycom, only: hycom_t, hycom_initialize, hycom_convert, hycom_finalize
!
   implicit none
!
   type(khkim_t) :: k_in
   type(hycom_t) :: h_in
   type(scrip_t) :: s_in
!
   integer(i4)        :: method   ! 2 (khkim), 3(hycom)
   character(len=512) :: filename
!
   method   = 2
   filename = './khkim/cs_grid_voronoi_ne030np4_rotated.nc'
   call read_nl()
!
   if     (method.eq.2) then
     call khkim_initialize(k_in,trim(filename))
     call khkim_convert(k_in)
     call khkim_finalize(k_in)
   elseif (method.eq.3) then
     call hycom_initialize(h_in,trim(filename))
     call hycom_convert(h_in)
     call hycom_finalize(h_in)
   elseif (method.eq.4) then
     call scrip_initialize_read(s_in,trim(filename))
     call scrip_check(s_in)
     call scrip_finalize(s_in)
   else
     print*,'check method = ', method
   endif
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_nl()
!-------------------------------------------------------------------------------
   implicit none
! local variables
   namelist /inputs/ method, filename
!
   read(*,nml=inputs)
!
   print*,'# options #'
   print*,' - method   = ',method
   print*,' - filename = ',trim(filename)
   print*,' '
!
   return
   end subroutine read_nl
!-------------------------------------------------------------------------------
   end program genscrip_input
!-------------------------------------------------------------------------------
