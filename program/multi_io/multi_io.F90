!-------------------------------------------------------------------------------
   module multi_io
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    2016-05-18  junghan kim    initial setup
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kinds, only: i4, l4, r4, r8
   use netcdf
!-------------------------------------------------------------------------------
   implicit none
   include 'mpif.h'
!
   private
!
   integer(i4)        :: nfiles
   integer(i4)        :: ncols, nlevs, nvars
   character(len=256) :: filepath
   character(len= 16) :: prefix
!
   real(r8), dimension(:,:), allocatable :: var
   real(r8) :: wtime
!
   public :: io_initialize, io_finalize, write_file, read_file, copy_file, get_io_wtime
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine io_initialize(counts, nv, nx, ny, filedir)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),      intent(in   ) :: counts, nv, nx, ny
   character(len=*), intent(in   ) :: filedir
!
   nfiles   = counts
   nvars    = nv
   ncols    = nx
   nlevs    = ny
   filepath = trim(filedir)
   prefix   = 'file-'
   allocate(var(ncols,nlevs))
   call random_number(var) !var = 0.0_r8
!   wtime = 0.0_r8
!
   return
   end subroutine io_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine io_finalize()
!-------------------------------------------------------------------------------
   implicit none
!
   deallocate(var)
!
   return
   end subroutine io_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine gen_filename(id, filename)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),      intent(in   ) :: id
   character(len=*), intent(inout) :: filename
!
   write(filename,'(a,a,a,i4.4,a)') trim(filepath),'/',trim(prefix),id,'.nc'
!
   return
   end subroutine gen_filename
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_io_wtime()
!-------------------------------------------------------------------------------
   implicit none
   real(r8) :: get_io_wtime
!
   get_io_wtime = wtime
!
   end function get_io_wtime
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_file(id)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4), intent(in   ) :: id
! local variables
   integer(i4)              :: i, j
   integer(i4)              :: err, fid, h_did, v_did
   integer(i4), allocatable :: vids(:)
   character(len= 16)       :: varname
   character(len=256)       :: filename
   real(r8)                 :: st, et
!
   allocate(vids(nvars))
   !write(filename,'(a,a,a,i3.3)') trim(filepath),'/',trim(prefix),id,'.nc'
   call gen_filename(id, filename)
   st = mpi_wtime(err)
   err = nf90_create(trim(filename),ior(nf90_mpiio,ior(nf90_clobber,nf90_netcdf4)),fid)
   err = nf90_def_dim(fid,'ncols',ncols,h_did)
   err = nf90_def_dim(fid,'nlevs',nlevs,v_did)
   do i = 1,nvars
     write(varname,'(a,i3.3)') 'var',i
     err = nf90_def_var(fid,trim(varname),nf90_double,(/h_did,v_did/),vids(i))
   enddo
   err = nf90_enddef(fid)
   do i = 1,nvars
     err = nf90_put_var(fid,vids(i),var)
   enddo
   err = nf90_close(fid)
   et = mpi_wtime(err)
   deallocate(vids)
   wtime = wtime+(et-st)
!
   return
   end subroutine write_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_file(id)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4), intent(in   ) :: id
! local variables
   integer(i4)              :: i, j
   integer(i4)              :: err, fid, h_did, v_did
   integer(i4), allocatable :: vids(:)
   character(len= 16)       :: varname
   character(len=256)       :: filename
   real(r8)                 :: st, et
!
   allocate(vids(nvars))
   call gen_filename(id, filename)
   st = mpi_wtime(err)
   err = nf90_open(trim(filename),nf90_nowrite,fid)
   if (err.ne.nf90_noerr) then
     print *,'no file...'
     stop
   endif
   do i = 1,nvars
     write(varname,'(a,i3.3)') 'var',i
     err = nf90_inq_varid(fid,trim(varname),vids(i))
   enddo
   do i = 1,nvars
     err = nf90_get_var(fid,vids(i),var)
   enddo
   err = nf90_close(fid)
   et = mpi_wtime(err)
   deallocate(vids)
   wtime = wtime+(et-st)
!
   return
   end subroutine read_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine copy_file(id)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4), intent(in   ) :: id
! local variables
   integer(i4)              :: err
   character(len=256)       :: filename
   real(r8)                 :: st, et
!
   call gen_filename(id, filename)
   st = mpi_wtime(err)
   call system('cp '//trim(filename)//' /scratch/jhkim/data/multi_io_cp')
   et = mpi_wtime(err)
   wtime = wtime+(et-st)
!
   return
   end subroutine copy_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module multi_io
!-------------------------------------------------------------------------------
