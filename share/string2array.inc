   character(len=*),                       intent(in   ) :: string
   character(len=1), optional,             intent(in   ) :: divisor
! local variables
   character(len=1) :: div
   integer(i4)      :: len_str, ndivs
   integer(i4)      :: i, is, ie
   integer(i4), dimension(100) :: idivs
!
   if (present(divisor)) then
     div = divisor
   else
     div = ','
   endif
!
   len_str = len(string)
   ndivs = 0
   do i = 1,len_str
     if (string(i:i).eq.div) then
       ndivs = ndivs+1
       idivs(ndivs) = i
     endif
   enddo
!
   allocate(array(ndivs+1))
   if (ndivs.eq.0) then
     read(string,*) array(1)
   else
     do i = 1, ndivs+1
       if (i.eq.1) then
         is =            1; ie = idivs(i)-1
       elseif (i.eq.ndivs+1) then
         is = idivs(i-1)+1; ie = len_str
       else
         is = idivs(i-1)+1; ie = idivs(i)-1
       endif
       !
       read(string(is:ie),*) array(i)
     enddo
   endif
