!-------------------------------------------------------------------------------
   program test
   use kinds,    only: i4, r8
   use w3odatmd, only: iaproc, naproc, mod1_ini=>initialize, mod1_fin=>finalize
   use w3gdatmd, only: nx, ny, nsea, nseal, mapsf, mapsta,                     &
                       mod2_ini=>initialize, mod2_fin=>finalize
   use wav_comp_driver, only: make_model_map, convert_wav2cpl, convert_cpl2wav
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    201x-0x-xx  junghan kim    initial setup
!
!  structure:
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4) :: gsize, lsize, i
   integer(i4), dimension(:), allocatable :: map, map_mask
!
   call mod1_ini()
   call mod2_ini()
!
   if (iaproc.eq.1) then
     print *, ' nx, ny      = ', nx, ny
     print *, ' nsea, nseal = ', nsea, nseal
     print *, ' mapsta = '
   endif
   do i = 1,ny
     print *, mapsta(i,:)
   enddo
   print *, ' '
!
   call make_model_map(gsize,lsize,map,map_mask)
!
   if (iaproc.eq.1) then
     print *, ' gsize, lsize, iproc = ', gsize, lsize, iaproc
     print *, ' '
     print *, ' map = '
     print *, map(:)
     print *, ' '
     print *, ' map_mask = '
     print *, map_mask(:)
     print *, ' '
   endif
!
   call mod2_fin()
   call mod1_fin()
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine test_sub(n)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: n
! local variables
   integer(i4) :: i
!
!
   return
   end subroutine test_sub
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program test
!-------------------------------------------------------------------------------
