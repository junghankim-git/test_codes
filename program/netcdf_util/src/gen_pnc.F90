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
   implicit none
   include 'mpif.h'
   integer(i4) :: j2, j3, j4, n1, n2, n3, n4
   integer(i4) :: comm, nprocs, rank, err
   integer(i4) :: fid, d1, d2, d3, d4
   integer(i4) :: i41, i42, i43, i44, r41, r42, r43, r44, r81, r82, r83, r84
   integer(i4) :: n1d, n2d, n3d, n4d
   !integer(i4) :: vids(12)
   integer(i4), dimension(:)      , allocatable :: i4_1d
   integer(i4), dimension(:,:)    , allocatable :: i4_2d
   integer(i4), dimension(:,:,:)  , allocatable :: i4_3d
   integer(i4), dimension(:,:,:,:), allocatable :: i4_4d
   real(r4)   , dimension(:)      , allocatable :: r4_1d
   real(r4)   , dimension(:,:)    , allocatable :: r4_2d
   real(r4)   , dimension(:,:,:)  , allocatable :: r4_3d
   real(r4)   , dimension(:,:,:,:), allocatable :: r4_4d
   real(r8)   , dimension(:)      , allocatable :: r8_1d
   real(r8)   , dimension(:,:)    , allocatable :: r8_2d
   real(r8)   , dimension(:,:,:)  , allocatable :: r8_3d
   real(r8)   , dimension(:,:,:,:), allocatable :: r8_4d
!
   call mpi_init(err)
   comm = mpi_comm_world
   call mpi_comm_size(comm,nprocs,err)
   call mpi_comm_rank(comm,  rank,err)
!
#if 0
   n1 = 20
   n2 = 2
   n3 = 2
   n4 = 2
#else
   !n1 = 125000000
   n1 = 130000000
   n2 = 2
   n3 = 2
   n4 = 1
#endif
   n1d = int(n1,8)
   n2d = int(n1,8)*int(n2,8)
   n3d = int(n1,8)*int(n2,8)*int(n3,8)
   n4d = int(n1,8)*int(n2,8)*int(n3,8)*int(n4,8)
   print *,'n1 = ',n1
   print *,'n2 = ',n2
   print *,'n3 = ',n3
   print *,'n4 = ',n4
   print *,'huge(i4) = ',huge(n1)
   print *,'total n  = ',n4d
   print *,'1d size(4byte) = ',4_i8*n1d/1024./1024./1024.,'GiB'
   print *,'2d size(4byte) = ',4_i8*n2d/1024./1024./1024.,'GiB'
   print *,'3d size(4byte) = ',4_i8*n3d/1024./1024./1024.,'GiB'
   print *,'4d size(4byte) = ',4_i8*n4d/1024./1024./1024.,'GiB'
!
   allocate(i4_1d(n1))
   allocate(i4_2d(n1,n2))
   allocate(i4_3d(n1,n2,n3))
   allocate(i4_4d(n1,n2,n3,n4))
   allocate(r4_1d(n1))
   allocate(r4_2d(n1,n2))
   allocate(r4_3d(n1,n2,n3))
   allocate(r4_4d(n1,n2,n3,n4))
   allocate(r8_1d(n1))
   allocate(r8_2d(n1,n2))
   allocate(r8_3d(n1,n2,n3))
   allocate(r8_4d(n1,n2,n3,n4))
   i4_1d(:)       = 1000
   i4_2d(:,:)     = 2000
   i4_3d(:,:,:)   = 3000
   i4_4d(:,:,:,:) = 4000
   r4_1d(:)       = 1000._r4
   r4_2d(:,:)     = 2000._r4
   r4_3d(:,:,:)   = 3000._r4
   r4_4d(:,:,:,:) = 4000._r4
   r8_1d(:)       = 1000._r8
   r8_2d(:,:)     = 2000._r8
   r8_3d(:,:,:)   = 3000._r8
   r8_4d(:,:,:,:) = 4000._r8
   do j2 = 1,n2
     i4_2d(:,j2) = i4_2d(:,j2)+j2
     r4_2d(:,j2) = r4_2d(:,j2)+j2
     r8_2d(:,j2) = r8_2d(:,j2)+j2
   enddo
   do j3 = 1,n3
   do j2 = 1,n2
     i4_3d(:,j2,j3) = i4_3d(:,j2,j3)+j2+10*j3
     r4_3d(:,j2,j3) = r4_3d(:,j2,j3)+j2+10*j3
     r8_3d(:,j2,j3) = r8_3d(:,j2,j3)+j2+10*j3
   enddo
   enddo
   do j4 = 1,n4
   do j3 = 1,n3
   do j2 = 1,n2
     i4_4d(:,j2,j3,j4) = i4_4d(:,j2,j3,j4)+j2+10*j3+100*j4
     r4_4d(:,j2,j3,j4) = r4_4d(:,j2,j3,j4)+j2+10*j3+100*j4
     r8_4d(:,j2,j3,j4) = r8_4d(:,j2,j3,j4)+j2+10*j3+100*j4
   enddo
   enddo
   enddo
!
   err = nf90mpi_create(comm,'output/pnc.nc',ior(nf90mpi_clobber,nf90_64bit_data), &
                                                                 mpi_info_null,fid)  
   call nc_check(err,'create')
!
!   integer(i4) :: i41, i42, i43, i44, r41, r42, r43, r44, r81, r82, r83, r84
   err = nf90mpi_def_dim(fid,'n1',int(n1,8),d1)
   err = nf90mpi_def_dim(fid,'n2',int(n2,8),d2)
   err = nf90mpi_def_dim(fid,'n3',int(n3,8),d3)
   err = nf90mpi_def_dim(fid,'n4',int(n4,8),d4)
   err = nf90mpi_put_att(fid,nf90mpi_global,'version','v1.2')
   err = nf90mpi_put_att(fid,nf90mpi_global,'i4_att',9)
   err = nf90mpi_put_att(fid,nf90mpi_global,'i8_att',9_i8)
   err = nf90mpi_put_att(fid,nf90mpi_global,'r4_att',9.8)
   err = nf90mpi_put_att(fid,nf90mpi_global,'r8_att',9.8_r8)
   err = nf90mpi_def_var(fid,'i4_1d',nf90mpi_int,(/d1/),i41)
   err = nf90mpi_def_var(fid,'i4_2d',nf90mpi_int,(/d1,d2/),i42)
   err = nf90mpi_def_var(fid,'i4_3d',nf90mpi_int,(/d1,d2,d3/),i43)
   err = nf90mpi_def_var(fid,'i4_4d',nf90mpi_int,(/d1,d2,d3,d4/),i44)
   err = nf90mpi_def_var(fid,'r4_1d',nf90mpi_real,(/d1/),r41)
   err = nf90mpi_def_var(fid,'r4_2d',nf90mpi_real,(/d1,d2/),r42)
   err = nf90mpi_def_var(fid,'r4_3d',nf90mpi_real,(/d1,d2,d3/),r43)
   err = nf90mpi_def_var(fid,'r4_4d',nf90mpi_real,(/d1,d2,d3,d4/),r44)
   err = nf90mpi_def_var(fid,'r8_1d',nf90mpi_real8,(/d1/),r81)
   err = nf90mpi_def_var(fid,'r8_2d',nf90mpi_real8,(/d1,d2/),r82)
   err = nf90mpi_def_var(fid,'r8_3d',nf90mpi_real8,(/d1,d2,d3/),r83)
   err = nf90mpi_def_var(fid,'r8_4d',nf90mpi_real8,(/d1,d2,d3,d4/),r84)
!
   err = nf90mpi_put_att(fid,i41,'unit','m/s')
   err = nf90mpi_put_att(fid,i41,'i4_type',2)
   err = nf90mpi_put_att(fid,i41,'r4_type',1.99)
   err = nf90mpi_put_att(fid,r42,'unit','m/s')
   err = nf90mpi_put_att(fid,r42,'i4_type',2)
   err = nf90mpi_put_att(fid,r42,'r4_type',1.99)
   err = nf90mpi_put_att(fid,r84,'unit','m/s')
   err = nf90mpi_put_att(fid,r84,'i4_type',2)
   err = nf90mpi_put_att(fid,r84,'r4_type',1.99)
   err = nf90mpi_enddef(fid)
!
   err = nfmpi_begin_indep_data(fid)
   call nc_check(err,'start: indep')
!
   print *,'i4_1d'
   err = nfmpi_put_var(fid,i41,i4_1d,int(n1,8),mpi_integer)
   call nc_check(err,'put_var')
   print *,'i4_2d'
   err = nfmpi_put_var(fid,i42,i4_2d,int(n1*n2,8),mpi_integer)
   call nc_check(err,'put_var')
   print *,'i4_3d'
   err = nfmpi_put_var(fid,i43,i4_3d,int(n1*n2*n3,8),mpi_integer)
   call nc_check(err,'put_var')
   print *,'i4_4d'
   err = nfmpi_put_var(fid,i44,i4_4d,int(n1*n2*n3*n4,8),mpi_integer)
   call nc_check(err,'put_var')
   print *,'r4_1d'
   err = nfmpi_put_var(fid,r41,r4_1d,int(n1,8),mpi_real)
   call nc_check(err,'put_var')
   print *,'r4_2d'
   err = nfmpi_put_var(fid,r42,r4_2d,int(n1*n2,8),mpi_real)
   call nc_check(err,'put_var')
   print *,'r4_3d'
   err = nfmpi_put_var(fid,r43,r4_3d,int(n1*n2*n3,8),mpi_real)
   call nc_check(err,'put_var')
   print *,'r4_4d'
   err = nfmpi_put_var(fid,r44,r4_4d,int(n1*n2*n3*n4,8),mpi_real)
   call nc_check(err,'put_var')
   print *,'r8_1d'
   err = nfmpi_put_var(fid,r81,r8_1d,int(n1,8),mpi_real8)
   call nc_check(err,'put_var')
   print *,'r8_2d'
   err = nfmpi_put_var(fid,r82,r8_2d,int(n1*n2,8),mpi_real8)
   call nc_check(err,'put_var')
   print *,'r8_3d'
   err = nfmpi_put_var(fid,r83,r8_3d,int(n1*n2*n3,8),mpi_real8)
   call nc_check(err,'put_var')
   print *,'r8_4d'
   err = nfmpi_put_var(fid,r84,r8_4d,int(n1*n2*n3*n4,8),mpi_real8)
   call nc_check(err,'put_var')
!
   err = nfmpi_end_indep_data(fid)
   call nc_check(err,'end  : indep')
!
   err = nf90mpi_close(fid)
!
   deallocate(i4_1d)
   deallocate(i4_2d)
   deallocate(i4_3d)
   deallocate(i4_4d)
   deallocate(r4_1d)
   deallocate(r4_2d)
   deallocate(r4_3d)
   deallocate(r4_4d)
   deallocate(r8_1d)
   deallocate(r8_2d)
   deallocate(r8_3d)
   deallocate(r8_4d)
   call mpi_finalize(err)
!
   contains
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
