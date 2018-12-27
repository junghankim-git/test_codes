!-------------------------------------------------------------------------------
   program main_program
!-------------------------------------------------------------------------------
!
!  abstract : 
!
!  history log :
!    2018-09-21  junghan kim    initial setup
!
!  structure :
!
!-------------------------------------------------------------------------------
   use kinds,   only: i4, i8, r4, r8, l4
   use pnetcdf
   use pnetcdf, only: nf90mpi_clobber=>nf90_clobber,nf90mpi_noerr=>nf90_noerr, &
                      nf90mpi_64bit_offset=>nf90_64bit_offset,                 &
                      nf90mpi_nowrite=>nf90_nowrite,                           &
                      nf90mpi_global=>nf90_global,nf90mpi_int=>nf90_int,       &
                      nf90mpi_real=>nf90_real,nf90mpi_real8=>nf90_real8,       &
                      nf90mpi_char=>nf90_char
   include 'mpif.h'
   integer(i8) :: nx, ny
   integer(i4) :: comm, nprocs, rank, err
!
   call mpi_init(err)
   comm = mpi_comm_world
   call mpi_comm_size(comm,nprocs,err)
   call mpi_comm_rank(comm,  rank,err)
!
   nx = 14155778
   !ny = 3
   ny = 12*16
!
   call write_pnc(nx,ny)
!
   call mpi_finalize(err)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_pnc(n, m)
   implicit none
   integer(i8), intent(in   ) :: n, m
! local variables
   integer, dimension(:)  , allocatable :: var1
   real   , dimension(:)  , allocatable :: var2
   real   , dimension(:,:), allocatable :: var3
   integer :: i, j
   integer :: fid, did1, did2, vid1, vid2, vid3
!
   write(*,'(a,i14)')    'n x m         = ',n*m
   write(*,'(a,f5.2,a)') 'variable size = ',4.*n*m/1000./1000./1000.,' GB'
!
   allocate(var1(n))
   allocate(var2(n))
   allocate(var3(n,m))
!
   err = nf90mpi_create(comm,'pnc.nc',ior(nf90mpi_clobber,nf90_64bit_data), &
                        mpi_info_null,fid)  
   call nc_check(err,'create')
!
   !err = nf90mpi_def_dim(fid,'n',int(n,8),did1)
   !err = nf90mpi_def_dim(fid,'m',int(m,8),did2)
   err = nf90mpi_def_dim(fid,'n',n,did1)
   err = nf90mpi_def_dim(fid,'m',m,did2)
   err = nf90mpi_def_var(fid,'var1',nf90mpi_int,(/did1/),vid1)
   err = nf90mpi_def_var(fid,'var2',nf90mpi_real,(/did1/),vid2)
   err = nf90mpi_def_var(fid,'var3',nf90mpi_real,(/did1,did2/),vid3)
   err = nf90mpi_put_att(fid,nf90mpi_global,'version','v1.2')
   err = nf90mpi_put_att(fid,nf90mpi_global,'i4_type',9)
   err = nf90mpi_put_att(fid,nf90mpi_global,'r4_type',9.8)
   err = nf90mpi_put_att(fid,vid1,'unit','m/s')
   err = nf90mpi_put_att(fid,vid1,'i4_type',2)
   err = nf90mpi_put_att(fid,vid1,'r4_type',1.99)
   err = nf90mpi_put_att(fid,vid2,'unit','m/s')
   err = nf90mpi_put_att(fid,vid2,'i4_type',2)
   err = nf90mpi_put_att(fid,vid2,'r4_type',1.99)
   err = nf90mpi_put_att(fid,vid3,'unit','m/s')
   err = nf90mpi_put_att(fid,vid3,'i4_type',2)
   err = nf90mpi_put_att(fid,vid3,'r4_type',1.99)
   err = nf90mpi_enddef(fid)
!
   err = nfmpi_begin_indep_data(fid)
   call nc_check(err,'start: indep')
   var1 = 1
   var2 = 0.2
   var3 = 0.3
   err = nfmpi_put_var(fid,vid1,var1,int(n,8),mpi_integer)
   err = nfmpi_put_var(fid,vid2,var2,int(n,8),mpi_real)
!   err = nfmpi_put_var(fid,vid3,var3,int(n*m,8),mpi_real)
   do j = 1,m
     write(*,'(a,i3,a)')'write ',j,'th level'
     err = nfmpi_put_vara(fid,vid3,(/1_i8,int(j,8)/),(/int(n,8),1_i8/),        &
                          var3(:,j),int(n,8),mpi_real)
     call nc_check(err,'put_var_sub')
   enddo
   call nc_check(err,'put_var')
   err = nfmpi_end_indep_data(fid)
   call nc_check(err,'end  : indep')
!
   err = nf90mpi_close(fid)
   deallocate(var1)
   deallocate(var2)
   deallocate(var3)
!
   end subroutine write_pnc
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
   if (err.ne.nf90mpi_noerr) then
     if (present(message)) print*,trim(message)//': '//trim(nf90mpi_strerror(err))
     stop
   endif
!
   end subroutine nc_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program main_program
!-------------------------------------------------------------------------------
