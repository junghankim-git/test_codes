!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
   use kinds, only : i4, r4, r8
   use kim_std_vars, only : nc_check
   use netcdf
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4) :: err, nhoriz, ihoriz, si, ei
   integer(i4) :: fileid, dimid, varid
   character(len=512) :: ofilename, nfilename
   real(r8), dimension(:), allocatable :: var, lats, lons
!
   ofilename = './test/org_5.nc'
   nfilename = './test/mod_5.nc'
!
   ! original file
   err = nf90_open(trim(ofilename),nf90_nowrite,fileid)
   call nc_check(err)
   err = nf90_inq_dimid(fileid,'ncol',dimid)
   call nc_check(err)
   err = nf90_inquire_dimension(fileid,dimid,len=nhoriz)
   call nc_check(err)
   allocate(var(nhoriz))
   allocate(lons(nhoriz))
   allocate(lats(nhoriz))
   ! lon & lat
   err = nf90_inq_varid(fileid,'lon',varid)
   call nc_check(err)
   err = nf90_get_var(fileid,varid,lons)
   call nc_check(err)
   err = nf90_inq_varid(fileid,'lat',varid)
   call nc_check(err)
   err = nf90_get_var(fileid,varid,lats)
   call nc_check(err)
   ! var
   err = nf90_inq_varid(fileid,'tsfc',varid)
   call nc_check(err)
   err = nf90_get_var(fileid,varid,var)
   call nc_check(err)
   err = nf90_close(fileid)
   call nc_check(err)
!
   ! new file
   err = nf90_open(trim(nfilename),nf90_write,fileid)
   call nc_check(err)
!
! modification here
! method 1
!  var(101) = 127.0_r8
!  var(102) = 127.0_r8
!  var(201) = 127.0_r8
!  var(321) = 627.0_r8

! method 2
!  si = 500
!  ei = 520
!  DO ihoriz = si,ei
!    var(ihoriz) = var(ihoriz)*0.9_r8
!  END DO

! method 3
!  DO ihoriz = 1,nhoriz
!    IF (lats(ihoriz) .LT. -80.0_r8) THEN
!      var(ihoriz) = var(ihoriz)*0.9_r8
!    END IF
!  END DO

! method 4
   do ihoriz = 1,nhoriz
     if (lats(ihoriz).gt.-1.0_r8.and.lats(ihoriz).lt.1.0_r8) then
       var(ihoriz) = var(ihoriz)-50.0_r8
     endif
   enddo
!
   err = nf90_inq_varid(fileid,'tsfc',varid)
   call nc_check(err)
   err = nf90_put_var(fileid,varid,var)
   call nc_check(err)
!
   err = nf90_close(fileid)
   call nc_check(err)
!
   deallocate(var)
   deallocate(lons)
   deallocate(lats)
!
   end
!-------------------------------------------------------------------------------
