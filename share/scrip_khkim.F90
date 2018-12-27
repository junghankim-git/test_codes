!-------------------------------------------------------------------------------
   module scrip_khkim
!-------------------------------------------------------------------------------
!
!  abstract : 
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds, only : i4, l4, r4, r8
   use scrip_input, only : scrip_t
   use scrip_input, only : scrip_initialize, scrip_write, scrip_finalize, scrip_check
   use netcdf
!
   private
!
   type khkim_t
     ! filename
     character(len=512) :: filename
     ! attributes
     integer(i4) :: ne
     logical(l4) :: rotated
     ! dimensions
     integer(i4) :: up_size
     integer(i4) :: corner_max_size
     ! variables
     real(r8), dimension(:,:), allocatable :: center_latlons ! 2, up_size
     real(r8), dimension(:,:, :), allocatable :: corner_latlons ! 2, corner_max_size, up_size
     ! scrip input
     type(scrip_t) :: scrip
   end type khkim_t
!
   real(r8), parameter :: pi = 3.141592653589793238462643383279_r8
!
   public :: khkim_t, khkim_initialize, khkim_convert, khkim_finalize
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine khkim_initialize(input, infilename)
!-------------------------------------------------------------------------------
   implicit none
!
   type(khkim_t), intent(inout) :: input
   character(len=*),  intent(in   ) :: infilename!, oufilename
! local variables
   integer(i4) :: fid, nup_did, ncor_did, vid, ierr
   character(len=16) :: rotated
   character(len=512) :: oufilename, description
!
   input%filename = trim(infilename)
! open
   ierr = nf90_open(trim(infilename),nf90_nowrite,fid)
   call nc_check(ierr)
! read dimensions
   ierr = nf90_inq_dimid(fid,'up_size',nup_did)
   call nc_check(ierr)
   ierr = nf90_inquire_dimension(fid,nup_did,len=input%up_size)
   call nc_check(ierr)
   ierr = nf90_inq_dimid(fid,'corner_max_size',ncor_did)
   call nc_check(ierr)
   ierr = nf90_inquire_dimension(fid,ncor_did,len=input%corner_max_size)
   call nc_check(ierr)
! allocate
   allocate(input%center_latlons(2,input%up_size))
   allocate(input%corner_latlons(2,input%corner_max_size,input%up_size))
! read variables
   ierr = nf90_inq_varid(fid,'center_latlons',vid)
   call nc_check(ierr)
   ierr = nf90_get_var(fid,vid,input%center_latlons)
   call nc_check(ierr)
   ierr = nf90_inq_varid(fid,'corner_latlons',vid)
   call nc_check(ierr)
   ierr = nf90_get_var(fid,vid,input%corner_latlons)
   call nc_check(ierr)
! read global attributes
   ierr = nf90_get_att(fid,nf90_global,'ne',input%ne)
   call nc_check(ierr)
   ierr = nf90_get_att(fid,nf90_global,'rotated',rotated)
   call nc_check(ierr)
!
   if (trim(rotated).eq.'true') then
     input%rotated = .true.
   else
     input%rotated = .false.
   endif
! close
   ierr = nf90_close(fid)
   call nc_check(ierr)
! initialize SCRIP_Input
   if (input%rotated) then
     write(oufilename,'(a,i3.3,a) ') '../CS/khkim_ne',input%ne,'_rotated.nc'
     write(description,'(a,i3.3,a) ') 'Cubed-sphere on ne',input%ne,' and np4(rotated) '
   else
     write(oufilename,'(a,i3.3,a) ') '../CS/khkim_ne',input%ne,'_unrotated.nc'
     write(description,'(a,i3.3,a) ') 'Cubed-sphere on ne',input%ne,' and np4(unrotated) '
   endif
   call scrip_initialize(input%scrip,trim(oufilename),1,input%up_size,input%corner_max_size,trim(description))
!
   return
   end subroutine khkim_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine khkim_finalize(input)
!-------------------------------------------------------------------------------
   implicit none
!
   type(khkim_t), intent(inout) :: input
!
   call scrip_finalize(input%scrip)
   deallocate(input%center_latlons)
   deallocate(input%corner_latlons)
!
   return
   end subroutine khkim_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine khkim_convert(input)
!-------------------------------------------------------------------------------
   implicit none
!
   type(khkim_t), intent(inout) :: input
! local variables
   integer(i4) :: iup, ico, ncorners
!
! dimension
   input%scrip%grid_dims(:) = 1
!
! centers
   do iup = 1,input%up_size
   input%scrip%grid_imask(iup) = 1
!radian
   input%scrip%grid_center_lat(iup) = input%center_latlons(1,iup)
   input%scrip%grid_center_lon(iup) = input%center_latlons(2,iup)
   enddo
! corners
   ncorners = input%corner_max_size
   do iup = 1,input%up_size
     do ico = 1,input%corner_max_size
#if 0
       if (input%corner_latlons(1,ico,iup).lt.-pi/2.0_r8.or.                   &
                              input%corner_latlons(1,ico,iup).gt.pi/2.0_r8) then
         print*,'found lat',input%corner_latlons(1,ico,iup)
       endif
       if (input%corner_latlons(2,ico,iup).lt.0.0_r8.or.                       &
                              input%corner_latlons(2,ico,iup).gt.2.0_r8*pi) then
         print*,'found lon',input%corner_latlons(1,ico,iup)
       endif
#endif
       ! counter-clockwise
       input%scrip%grid_corner_lat(ico,iup) = input%corner_latlons(1,ico,iup)
       input%scrip%grid_corner_lon(ico,iup) = input%corner_latlons(2,ico,iup)
       ! remove zero
       if (input%scrip%grid_corner_lat(ico,iup).eq.0.0.and.input%scrip%grid_corner_lon(ico,iup).eq.0.0) then
         if (ico.eq.1) then
           input%scrip%grid_corner_lon(ico,iup) = input%scrip%grid_corner_lon(2,iup)
           input%scrip%grid_corner_lat(ico,iup) = input%scrip%grid_corner_lat(2,iup)
         else
           input%scrip%grid_corner_lon(ico,iup) = input%scrip%grid_corner_lon(ico-1,iup)
           input%scrip%grid_corner_lat(ico,iup) = input%scrip%grid_corner_lat(ico-1,iup)
         endif
       endif
     enddo
   enddo
!
   call scrip_check(input%scrip)
! write
   call scrip_write(input%scrip)
!
   return
!-------------------------------------------------------------------------------
   end subroutine khkim_convert
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_check(ncstat_in)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: ncstat_in
!
   if (ncstat_in.ne.nf90_noerr) then
     print*,trim(nf90_strerror(ncstat_in))
     stop
   endif
!
   return
   end subroutine nc_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module scrip_khkim
