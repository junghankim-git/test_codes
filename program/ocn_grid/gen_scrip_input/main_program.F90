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
   use scrip_hycom, only: hycom_t, hycom_initialize, hycom_convert, hycom_finalize
!
   implicit none
!
   type(hycom_t) :: h_in
!
   character(len=512) :: filename
!
   filename = './grid_240x193.nc'
!
   call hycom_initialize(h_in,trim(filename))
   call hycom_convert(h_in)
   call hycom_finalize(h_in)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program genscrip_input
!-------------------------------------------------------------------------------
