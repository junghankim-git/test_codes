!-------------------------------------------------------------------------------
   module mod_netcdf
!-------------------------------------------------------------------------------
!
!  abstract :  NetCDF
!
!  history log :
!    2017-08-08  junghan kim    initial setup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds,    only: i4, r4, r8, l4
   use netcdf
!-------------------------------------------------------------------------------
   private
   include 'mpif.h'
!
   public :: nc_initialize, nc_finalize
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_initialize(nc)
   implicit none
   integer(i4),                intent(inout) :: nc
!   character(len=*),           intent(in   ) :: infilename
!   character(len=*),           intent(in   ) :: oufilename
! local variables
   integer(i4) :: comm, nprocs, rank
   integer(i4) :: err, ndims, nvars, natts
!
   call mpi_init(err)
   comm = mpi_comm_world
   call mpi_comm_size(comm,nprocs,err)
   call mpi_comm_rank(comm,  rank,err)
!
   err = nf90_create('./test.nc',ior(nf90_mpiio,ior(nf90_netcdf4,nf90_clobber)),ncid=nc,comm=comm,info=mpi_info_null)
   print *, err
   call nc_check(err,'nc_initialize')
   !print *, 'infile : ', trim(infilename)
   !print *, ' - ndims ', ndims
   !print *, ' - nvars ', nvars
   !print *, ' - natts ', natts
!
   err = nf90_close(nc)
!
   end subroutine nc_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_finalize()
   implicit none
! local variables
   integer(i4) :: err
!
   call mpi_finalize(err)
!
   end subroutine nc_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_check(err,message)
   implicit none
   integer(i4),                intent(in   ) :: err
   character(len=*), optional, intent(in   ) :: message
! local variables
!
   if (err.ne.nf90_noerr) then
     if (present(message)) print *, trim(message)
     stop
   endif
!
   end subroutine nc_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module mod_netcdf
!-------------------------------------------------------------------------------
