   program main
   implicit none
   integer, parameter :: n = 10
   real(8), dimension(n+1) :: a_i, b_i
   real(8), dimension(n)   :: a_m, b_m

   a_i(:) = 1.0
   b_i(:) = 1.0
   a_m(:) = 1.0
   b_m(:) = 1.0

   call write_hybrid_coord(n,a_i,a_m,b_i,b_m)

   contains

!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_hybrid_coord(nlevs, a_i, a_m, b_i, b_m)
!-------------------------------------------------------------------------------
   implicit none
!
   integer,                     intent(in   ) :: nlevs
   real(8), dimension(nlevs+1), intent(in   ) :: a_i, b_i
   real(8), dimension(nlevs),   intent(in   ) :: a_m, b_m
!
   integer :: k
   character(len=128) :: filename_i, filename_m
!
   filename_i = './filename_half.dat'
   filename_m = './filename_full.dat'
!
   write(6,*) '* Output file(half level): ',trim(filename_i)
   write(6,*) '* Output file(full level): ',trim(filename_m)
!
   open(unit=71,file=filename_i,status='unknown',access='sequential',form='unformatted')
   open(unit=72,file=filename_m,status='unknown',access='sequential',form='unformatted')
!
   ! half levels      ! top to bottome
   do k = 1,nlevs+1
     write(71) a_i(k)
     write(71) b_i(k)
   enddo
   close(71)
!
   ! full levels      ! top to bottome
   do k = 1,nlevs
     write(72) a_m(k)
     write(72) b_m(k)
   enddo
   close(72)
!
   return
   end subroutine write_hybrid_coord
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------




   end program main
