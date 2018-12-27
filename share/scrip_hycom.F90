!------------------------------------------------------------------------------
   module scrip_hycom
!-------------------------------------------------------------------------------
!
!  abstract : 
!
!  history log :
!    2017-07-20  junghan kim    initial setup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds      , only: i4, l4, r4, r8
   use scrip_input, only: scrip_t
   use scrip_input, only: scrip_initialize, scrip_write, scrip_finalize, scrip_check
   use netcdf
!
   logical(l4) :: remove_pole = .false.
!
   type hycom_t
     ! filename
     character(len=512) :: filename
     ! attributes
     real(r8)    :: reflon, grdlon, half
     ! dimensions
     integer(i4) :: idm, jdm
     !integer(i4) :: corner_max_size
     ! variables
     real(r8), dimension(:,:), allocatable :: plon, plat, qlon, qlat
     ! scrip input
     type(scrip_t) :: scrip
   end type hycom_t
!
   !real(r8), parameter :: pi = 3.141592653589793238462643383279_r8
   real(r8), parameter :: pi = acos(-1.0_16)
!
   public :: hycom_t, hycom_initialize, hycom_convert, hycom_finalize
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine hycom_initialize(input, infilename)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hycom_t)   , intent(inout) :: input
   character(len=*), intent(in   ) :: infilename
! local variables
   integer(i4) :: fid, idid, jdid, vid, ierr
   integer(i4) :: idm, jdm, im, jm, ncorners, ndims
   character(len= 16) :: rotated
   character(len=512) :: oufilename, description
!
   input%filename = trim(infilename)
! open
   ierr = nf90_open(trim(infilename),nf90_nowrite,fid)
   call nc_check(ierr)
! read dimensions
   ierr = nf90_inq_dimid(fid,'idm',idid)
   call nc_check(ierr)
   ierr = nf90_inquire_dimension(fid,idid,len=input%idm)
   call nc_check(ierr)
   ierr = nf90_inq_dimid(fid,'jdm',jdid)
   call nc_check(ierr)
   ierr = nf90_inquire_dimension(fid,jdid,len=input%jdm)
   call nc_check(ierr)
! allocate
   idm = input%idm
   jdm = input%jdm
   allocate(input%plon(idm,jdm),input%plat(idm,jdm))
   allocate(input%qlon(idm,jdm),input%qlat(idm,jdm))
! read variables
   ierr = nf90_inq_varid(fid,'plon',vid)
   call nc_check(ierr)
   ierr = nf90_get_var(fid,vid,input%plon)
   call nc_check(ierr)
   ierr = nf90_inq_varid(fid,'plat',vid)
   call nc_check(ierr)
   ierr = nf90_get_var(fid,vid,input%plat)
   call nc_check(ierr)
   ierr = nf90_inq_varid(fid,'qlon',vid)
   call nc_check(ierr)
   ierr = nf90_get_var(fid,vid,input%qlon)
   call nc_check(ierr)
   ierr = nf90_inq_varid(fid,'qlat',vid)
   call nc_check(ierr)
   ierr = nf90_get_var(fid,vid,input%qlat)
   call nc_check(ierr)
! read global attributes
   ierr = nf90_get_att(fid,nf90_global,'reflon',input%reflon)
   call nc_check(ierr)
   ierr = nf90_get_att(fid,nf90_global,'grdlon',input%grdlon)
   call nc_check(ierr)
! close
   ierr = nf90_close(fid)
   call nc_check(ierr)
!
! convert deg. to rad.
   do jm = 1,input%jdm
     do im = 1,input%idm
       if (input%plon(im,jm).gt.360.0_r8) then
         input%plon(im,jm) = (input%plon(im,jm)-360.0_r8)*pi/180.0_r8
       else
         input%plon(im,jm) = input%plon(im,jm)*pi/180.0_r8
       endif
       input%plat(im,jm)   = input%plat(im,jm)*pi/180.0_r8
       if (input%qlon(im,jm).gt.360.0_r8) then
         input%qlon(im,jm) = (input%qlon(im,jm)-360.0_r8)*pi/180.0_r8
       else
         input%qlon(im,jm) = input%qlon(im,jm)*pi/180.0_r8
       endif
       input%qlat(im,jm)   = input%qlat(im,jm)*pi/180.0_r8
     enddo
   enddo
!
   ndims    = 2
   ncorners = 4
!
! initialize scrip_input
   if (remove_pole) then
     write(oufilename,'(a,i4.4,a,i4.4,a)')'../HYCOM/hycom_',idm,'x',jdm,'_wo_pole.nc'
     write(description,'(a,i4.4,a,i4.4,a)')'hycom grid: ',idm,'x',jdm,' (w/o pole)'
     call scrip_initialize(input%scrip,trim(oufilename),ndims,idm*(jdm-1),ncorners,trim(description))
   else
     write(oufilename,'(a,i4.4,a,i4.4,a)')'../HYCOM/hycom_',idm,'x',jdm,'.nc'
     write(description,'(a,i4.4,a,i4.4,a)')'hycom grid: ',idm,'x',jdm,' (w/ pole)'
     call scrip_initialize(input%scrip,trim(oufilename),ndims,idm*jdm,ncorners,trim(description))
   endif
!
   return
   end subroutine hycom_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine hycom_finalize(input)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hycom_t), intent(inout) :: input
!
   call scrip_finalize(input%scrip)
   deallocate(input%plon,input%plat)
   deallocate(input%qlon,input%qlat)
!
   return
   end subroutine hycom_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine hycom_convert(input)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hycom_t), intent(inout) :: input
! local variables
   integer(i4) :: idm, jdm, im, jm, nn, i, xx, yy
   real(r8)    :: xs(4), ys(4), tmp
!
   idm = input%idm
   jdm = input%jdm
   if (remove_pole) jdm = jdm-1
!
! dimension
   input%scrip%grid_dims(:) = (/idm,jdm/)
   !input%scrip%grid_dims(:) = (/idm*jdm/)
!
! centers
   do jm = 1,jdm
     do im = 1,idm
       nn = im+idm*(jm-1)
       input%scrip%grid_imask(nn) = 1
       input%scrip%grid_center_lon(nn) = input%plon(im,jm)
       input%scrip%grid_center_lat(nn) = input%plat(im,jm)
     enddo
   enddo
!
! corners
!   4     3
!      *
!   1     2
   do jm = 1,jdm
     do im = 1,idm
       nn = im+idm*(jm-1)
       if (remove_pole) then
         if (im.ne.idm) then
           xs(1) = input%qlon(im  ,jm  )
           ys(1) = input%qlat(im  ,jm  )
           xs(2) = input%qlon(im+1,jm  )
           ys(2) = input%qlat(im+1,jm  )
           xs(3) = input%qlon(im+1,jm+1)
           ys(3) = input%qlat(im+1,jm+1)
           xs(4) = input%qlon(im  ,jm+1)
           ys(4) = input%qlat(im  ,jm+1)
         elseif (im.eq.idm) then
           xs(1) = input%qlon(im  ,jm  )
           ys(1) = input%qlat(im  ,jm  )
           xs(2) = input%qlon(   1,jm  )
           ys(2) = input%qlat(   1,jm  )
           xs(3) = input%qlon(   1,jm+1)
           ys(3) = input%qlat(   1,jm+1)
           xs(4) = input%qlon(im  ,jm+1)
           ys(4) = input%qlat(im  ,jm+1)
         else
           print *, 'some error: what!'
           stop
         endif
       else
         if (im.ne.idm.and.jm.ne.jdm) then
           xs(1) = input%qlon(im  ,jm  )
           ys(1) = input%qlat(im  ,jm  )
           xs(2) = input%qlon(im+1,jm  )
           ys(2) = input%qlat(im+1,jm  )
           xs(3) = input%qlon(im+1,jm+1)
           ys(3) = input%qlat(im+1,jm+1)
           xs(4) = input%qlon(im  ,jm+1)
           ys(4) = input%qlat(im  ,jm+1)
         elseif (im.eq.idm.and.jm.ne.jdm) then
           xs(1) = input%qlon(im  ,jm  )
           ys(1) = input%qlat(im  ,jm  )
           xs(2) = input%qlon(   1,jm  )
           ys(2) = input%qlat(   1,jm  )
           xs(3) = input%qlon(   1,jm+1)
           ys(3) = input%qlat(   1,jm+1)
           xs(4) = input%qlon(im  ,jm+1)
           ys(4) = input%qlat(im  ,jm+1)
         elseif (im.ne.idm.and.jm.eq.jdm) then
           xs(1) = input%qlon(im  ,jm  )
           ys(1) = input%qlat(im  ,jm  )
           xs(2) = input%qlon(im+1,jdm  )
           ys(2) = input%qlat(im+1,jdm  )
           xs(3) = input%qlon(idm-im+1,jdm-1)
           ys(3) = input%qlat(idm-im+1,jdm-1)
           if (im.eq.1) then
           xs(4) = input%qlon(1,jdm-1)
           ys(4) = input%qlat(1,jdm-1)
           else
           xs(4) = input%qlon(idm-im+2,jdm-1)
           ys(4) = input%qlat(idm-im+2,jdm-1)
           endif
         elseif (im.eq.idm.and.jm.eq.jdm) then
           xs(1) = input%qlon(im  ,jm  )
           ys(1) = input%qlat(im  ,jm  )
           xs(2) = input%qlon(   1,jdm-1)
           ys(2) = input%qlat(   1,jdm-1)
           xs(3) = input%qlon(   1,jdm-1)
           ys(3) = input%qlat(   1,jdm-1)
           xs(4) = input%qlon(   2,jdm-1)
           ys(4) = input%qlat(   2,jdm-1)
         else
           print *, 'some error: what!'
           stop
         endif
       endif ! not remove pole
!
       do i = 1,4
         input%scrip%grid_corner_lon(i,nn) = xs(i)
         input%scrip%grid_corner_lat(i,nn) = ys(i)
       enddo
#if 0
       if (jm.eq.jdm) then
       !print *, im, jm, idm, jdm, nn, nn-im-im+1
       do i = 1,4
         input%scrip%grid_corner_lon(i,nn) = input%scrip%grid_corner_lon(i,nn-im-im+1)
         input%scrip%grid_corner_lat(i,nn) = input%scrip%grid_corner_lat(i,nn-im-im+1)
       enddo
       endif
#endif
!
     enddo
   enddo
!
   call scrip_check(input%scrip)
!
! write
   call scrip_write(input%scrip)
!
   return
   end subroutine hycom_convert
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine hycom_convert_no_dipole(input)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hycom_t), intent(inout) :: input
! local variables
   integer(i4) :: idm, jdm, im, jm, ii, i, xx, yy
   real(r8)    :: xs(4), ys(4), tmp
!
   idm = input%idm
   jdm = input%jdm
   if (remove_pole) jdm = jdm-1
!
! dimension
   input%scrip%grid_dims(:) = (/idm,jdm/)
!
! centers
   do jm = 1,jdm
     do im = 1,idm
       ii = im+idm*(jm-1)
       input%scrip%grid_imask(ii) = 1
       input%scrip%grid_center_lon(ii) = input%plon(im,jm)
       input%scrip%grid_center_lat(ii) = input%plat(im,jm)
     enddo
   enddo
!
! corners
!   4     3
!      *
!   1     2
   do jm = 1,jdm
     do im = 1,idm
       ii = im+idm*(jm-1)
       xs(1) = input%qlon(im,jm)
       xs(4) = input%qlon(im,jm)
       ys(1) = input%qlat(im,jm)
       ys(2) = input%qlat(im,jm)
       if (im.eq.idm) then
         tmp   = input%qlon(1,jm)
         xs(2) = tmp
         xs(3) = tmp
       else
         tmp   = input%qlon(im+1,jm)
         xs(2) = tmp
         xs(3) = tmp
       endif
       if (jm.eq.input%jdm) then
         tmp   = input%qlat(im,jm)+2.*(input%plat(im,jm)-input%qlat(im,jm))
         ys(3) = tmp
         ys(4) = tmp
       else
         tmp   = input%qlat(im,jm+1)
         ys(3) = tmp
         ys(4) = tmp
       endif
! careful
       if (ys(3).ge.pi/2.0_r8) then
         ys(3) = pi-ys(3)
         ys(4) = pi-ys(4)
         xs(3) = xs(3)+pi
         xs(4) = xs(4)+pi
         if (xs(3).ge.2.0_r8*pi) xs(3) = xs(3)-2.0_r8*pi
         if (xs(4).ge.2.0_r8*pi) xs(4) = xs(4)-2.0_r8*pi
         tmp   = xs(3)
         xs(3) = xs(4)
         xs(4) = tmp
       endif
! careful
       do i=1,4
         input%scrip%grid_corner_lon(i,ii) = xs(i)
         input%scrip%grid_corner_lat(i,ii) = ys(i)
       enddo
       ! some strange!! (counter-clockwise -> clockwise in pole)
       if (jm.eq.input%jdm) then
         !print *, 'DEBUG before: ', input%scrip%grid_corner_lon(:,ii)
         input%scrip%grid_corner_lon(2,ii) = xs(4)
         input%scrip%grid_corner_lat(2,ii) = ys(4)
         input%scrip%grid_corner_lon(4,ii) = xs(2)
         input%scrip%grid_corner_lat(4,ii) = ys(2)
         !print *, 'DEBUG after : ', input%scrip%grid_corner_lon(:,ii)
         !print *, ' '
       endif
     enddo
   enddo
!
   call scrip_check(input%scrip)
!
! write
   call scrip_write(input%scrip)
!
   return
   end subroutine hycom_convert_no_dipole
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
   end module scrip_hycom
!-------------------------------------------------------------------------------
