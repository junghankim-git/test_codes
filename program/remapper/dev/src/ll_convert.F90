!-------------------------------------------------------------------------------
   program ll_convert
!-------------------------------------------------------------------------------
!
!  abstract : 
!
!  history log :
!    2016-03-18  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!    2017-06-02  junghan kim    apply the remapping matrix of kim
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use netcdf
   use kinds,        only: i4, l4, r8
!-------------------------------------------------------------------------------
   implicit none
!
! variables
   character(len=512) :: filename
   integer(i4) :: err
   integer(i4) :: fid, did1, did2, did3, vid
   integer(i4) :: ncols, nsfcs, nlats, nlons
   real, dimension(:,:)  , allocatable :: inbuf
   real, dimension(:,:,:), allocatable :: oubuf
!
   filename = '2011072519'
!
   err = nf90_open(trim(filename)//'.nc',nf90_nowrite,fid); call check_nc_error(err)
   err = nf90_inq_dimid(fid,'ncol',did1); call check_nc_error(err)
   err = nf90_inq_dimid(fid,'nsfcs',did2); call check_nc_error(err)
   err = nf90_inquire_dimension(fid,did1,len=ncols); call check_nc_error(err)
   err = nf90_inquire_dimension(fid,did2,len=nsfcs); call check_nc_error(err)
   allocate(inbuf(ncols,nsfcs))
   err = nf90_inq_varid(fid,'sfcfcs',vid); call check_nc_error(err)
   err = nf90_get_var(fid,vid,inbuf); call check_nc_error(err)
   err = nf90_close(fid); call check_nc_error(err)
!
   nlons = 360
   nlats = 180
   if (ncols.ne.nlons*nlats) then
     print *, 'check dimension ... '
     stop
   endif
   allocate(oubuf(nlons,nlats,nsfcs))
   oubuf = reshape(inbuf,shape=(/nlons,nlats,nsfcs/))
!
   err = nf90_create(trim(filename)//'_ll.nc',ior(nf90_64bit_offset,nf90_clobber),fid); call check_nc_error(err)
   err = nf90_def_dim(fid,'nlons',nlons,did1); call check_nc_error(err)
   err = nf90_def_dim(fid,'nlats',nlats,did2); call check_nc_error(err)
   err = nf90_def_dim(fid,'nsfcs',nsfcs,did3); call check_nc_error(err)
   err = nf90_def_var(fid,'sfcfcs',nf90_real,(/did1,did2,did3/),vid); call check_nc_error(err)
   err = nf90_enddef(fid); call check_nc_error(err)
   err = nf90_put_var(fid,vid,oubuf)
   err = nf90_close(fid)
!
   deallocate(inbuf)
   deallocate(oubuf)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine check_nc_error(err)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: err
!
   if (err.ne.nf90_noerr) then
     print*,'error = ',nf90_strerror(err)
     stop
   endif
!
   return
   end subroutine check_nc_error
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program ll_convert
