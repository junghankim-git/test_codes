!-------------------------------------------------------------------------------
   module io_ascii
!-------------------------------------------------------------------------------
!
!  abstract :  bi-linear interpolation module
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-15  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   interface write_ascii
     module procedure write_ascii_l4_s
     module procedure write_ascii_i4_s
     module procedure write_ascii_r8_s
     module procedure write_ascii_l4_1d
     module procedure write_ascii_i4_1d
     module procedure write_ascii_r8_1d
     module procedure write_ascii_string
   end interface write_ascii
!
   interface write_ascii_lines
     module procedure write_ascii_r8_1d_lines
     module procedure write_ascii_2_r8_1d_lines
   end interface write_ascii_lines
!
   public :: open_file, close_file, write_ascii, write_ascii_lines
!   public :: write_ascii_l4_s, write_ascii_l4_1d
!   public :: write_ascii_i4_s, write_ascii_i4_1d
!   public :: write_ascii_r8_s, write_ascii_r8_1d
!   public :: write_ascii_string
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine open_file(filename, unitnum)
!-------------------------------------------------------------------------------
   implicit none
!
   character(*),intent(in   ) :: filename
   integer(i4), intent(in   ) :: unitnum
!
   open(unitnum,file=filename,status='UNKNOWN',form='FORMATTED')
!
   return
   end subroutine open_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine close_file(unitnum)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: unitnum
!
   close(unitnum)
!
   return
   end subroutine close_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ascii_i4_s(unitnum, value)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: unitnum
   integer(i4), intent(in   ) :: value
!
   write(unitnum,*) value
!
   return
   end subroutine write_ascii_i4_s
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ascii_i4_1d(unitnum, n, values)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: unitnum
   integer(i4),               intent(in   ) :: n
   integer(i4), dimension(n), intent(in   ) :: values
! local variables
   integer(i4) :: i
!
   write(unitnum,*) values(:)
!
   return
   end subroutine write_ascii_i4_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ascii_l4_s(unitnum, value)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: unitnum
   logical(l4), intent(in   ) :: value
!
   if (value) then
   write(unitnum,*) 1
   else
   write(unitnum,*) 0
   endif
!
   return
   end subroutine write_ascii_l4_s
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ascii_l4_1d(unitnum, n, values)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: unitnum
   integer(i4),               intent(in   ) :: n
   logical(l4), dimension(n), intent(in   ) :: values
! local variables
   integer(i4) :: i
   integer(i4), dimension(n) :: output
!
   do i = 1, n
     if (values(i)) then
       output(i) = 1
     else
       output(i) = 0
     endif
   enddo
   write(unitnum,*) output
!
   return
   end subroutine write_ascii_l4_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ascii_r8_s(unitnum, value)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: unitnum
   real(r8),    intent(in   ) :: value
!
   write(unitnum,*) value
!
   return
   end subroutine write_ascii_r8_s
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ascii_r8_1d(unitnum, n, values)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),            intent(in   ) :: unitnum
   integer(i4),            intent(in   ) :: n
   real(i8), dimension(n), intent(in   ) :: values
! local variables
   integer(i4) :: i
!
   write(unitnum,*) values(:)
!
   return
   end subroutine write_ascii_r8_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ascii_r8_1d_lines(unitnum, n, values)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),            intent(in   ) :: unitnum
   integer(i4),            intent(in   ) :: n
   real(i8), dimension(n), intent(in   ) :: values
! local variables
   integer(i4) :: i
!
   do i = 1,n
     write(unitnum,*) values(i)
   enddo
!
   return
   end subroutine write_ascii_r8_1d_lines
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ascii_2_r8_1d_lines(unitnum, n, values1, values2)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),            intent(in   ) :: unitnum
   integer(i4),            intent(in   ) :: n
   real(i8), dimension(n), intent(in   ) :: values1, values2
! local variables
   integer(i4) :: i
!
   do i = 1,n
     write(unitnum,*) values1(i),values2(i)
   enddo
!
   return
   end subroutine write_ascii_2_r8_1d_lines
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ascii_string(unitnum, string)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),  intent(in   ) :: unitnum
   character(*), intent(in   ) :: string
!
   write(unitnum,*) string
!
   return
   end subroutine write_ascii_string
#if 0
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_metadata(filename, testcase, n, zmin, zmax)
!-------------------------------------------------------------------------------
   implicit none
!
   character(*), intent(in   ) :: filename
   integer(i4),  intent(in   ) :: testcase, n
   real(r8),     intent(in   ) :: zmin, zmax
! local variables
   integer(i4) :: unitnum
!
   unitnum = 21
!
   open(unitnum,file=filename,status='UNKNOWN',form='FORMATTED')
!
   write(unitnum,*) testcase
   write(unitnum,*) n
   write(unitnum,*) zmin
   write(unitnum,*) zmax
!
   close(unitnum)
!
   return
   end subroutine write_metadata
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ascii_real(unitnum, n, values)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),            intent(in   ) :: unitnum
   integer(i4),            intent(in   ) :: n
   real(i8), dimension(n), intent(in   ) :: values
! local variables
   integer(i4) :: i
!
   do i = 1,n
     write(unitnum,*) values(i),0.0d0
   enddo
!
   return
   end subroutine write_ascii_real
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ascii_cmplx(unitnum, n, values)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: unitnum
   integer(i4),                  intent(in   ) :: n
   double complex, dimension(n), intent(in   ) :: values
! local variables
   integer(i4) :: i
!
   do i = 1,n
     write(unitnum,*) dble(values(i)),dimag(values(i))
   enddo
!
   return
   end subroutine write_ascii_cmplx
#endif
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module io_ascii
!-------------------------------------------------------------------------------
