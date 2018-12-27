!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    201x-0x-xx  junghan kim    initial setup
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kinds
   use statistics    , only: get_lerror
   use matrix_control, only: print_matrix, matrix_inverse, matrix_multiplicity,&
                             matrix_chk_identity, get_diagonal,                &
                             diagonal_inverse_interface, diagonal_solve_interface
   use matrix_control, only: fast_pentadiagonal_dummy
   use netcdf
   implicit none
   include 'mpif.h'
!
   integer(i4), parameter :: max_n  = 20
   integer(i4)            :: nn(max_n)            ! user input
   logical(l4)            :: do_solve = .false.   ! user input
   logical(l4)            :: tri      = .false.   ! user input
   logical(l4)            :: do_prt   = .false.   ! user input
   integer(i4)            :: nexp, nmet, ntry
   integer(i4)            :: ie, n, im, it, i, j, err
   integer(i4), dimension(:)  , allocatable :: ns
   real(r8)   , dimension(:,:), allocatable :: mat, imat, ret
   real(r8)   , dimension(:)  , allocatable :: a, b, c, d, e, ret0
   character(len=16), dimension(:), allocatable :: string
! performance
   integer(i4) :: perf_method = 0   ! 0(min), 1(average)
   real(r8)    :: lerr, stime, etime
   real(r8)   , dimension(:,:)  , allocatable :: l2err, elaps
   real(r8)   , dimension(:,:,:), allocatable :: elaps_t
!
   call read_namelist('inputs.nl')
   if (do_solve) then
     nmet = 3
   else
     nmet = 2
   endif
   allocate(l2err(nexp,nmet),elaps(nexp,nmet))
   allocate(elaps_t(ntry,nexp,nmet))
   allocate(string(nmet))
!
   if     (nmet==1) then
     string(1) = 'Fast algorithm'
   elseif (nmet==2) then
     string(1) = 'Fast algorithm'
     string(2) = 'LAPACK (general)'
   elseif (nmet==3) then
     string(1) = 'Fast algorithm'
     string(2) = 'LAPACK (general)'
     string(3) = 'LAPACK (band)'
   endif
   do ie = 1, nexp
     n = ns(ie)
     print '(a,i6.6,a)','------------------------------- ',n,' ------------------------------'
     call flush(6)
     !
     allocate(a(n),b(n),c(n))
     if (do_solve) then
       allocate(ret(n,nmet),ret0(n))
       do i=1,n
         ret0(i) = real(i,r8)
       enddo
     else
       allocate(mat(n,n),imat(n,n))
     endif
     if (tri) then
       call random_number(a)
       call random_number(b)
       call random_number(c)
       if (.not.do_solve) call get_diagonal(n,a,b,c,mat)
     else
       allocate(d(n),e(n))
       call random_number(a)
       call random_number(b)
       call random_number(c)
       call random_number(d)
       call random_number(e)
       if (.not.do_solve) call get_diagonal(n,a,b,c,mat,d,e)
     endif
     if (do_prt) call print_matrix('mat = ',n,n,mat)
     !
     do im = 1,nmet
       print *,'(info) ',trim(string(im)),' (start)'
       call flush(6)
       do it = 1,ntry
         print *,'  try: ',it
         call flush(6)
         if (do_solve) ret(:,im) = ret0(:)
         stime = mpi_wtime(err)
         if (do_solve) then
           if (tri) then
             call diagonal_solve_interface(im,3,n,a,b,c,ret(:,im))
           else
             call diagonal_solve_interface(im,5,n,a,b,c,ret(:,im),d,e)
           endif
         else
           if (tri) then
             call diagonal_inverse_interface(im,3,n,a,b,c,imat,mat)
           else
             call diagonal_inverse_interface(im,5,n,a,b,c,imat,mat,d,e)
           endif
         endif
         etime = mpi_wtime(err)
         if (do_solve) then
           l2err(ie,im) = 0.0_r8
         else
           if (do_prt) call print_matrix('imat = ',n,n,imat)
           if (it==1) then
             lerr  = matrix_chk_identity(n,mat,imat)
             l2err(ie,im) = lerr
           endif
         endif
         elaps_t(it,ie,im) = etime-stime
         if (perf_method==1) then
           elaps(ie,im) = sum(elaps_t(:,ie,im))/real(ntry,r8)
         else
           elaps(ie,im) = minval(elaps_t(:,ie,im))
         endif
         call print_performance(lerr,etime-stime)
         call flush(6)
       enddo ! do it = 1,ntry
       print *,'(info) ',trim(string(im)),' (done)'
       print *,''
       call flush(6)
     enddo ! do im = 1,nmet
     if (do_solve) then
       l2err(ie,1) = get_lerror(n,ret(:,2),ret(:,1),'L2')
       l2err(ie,2) = 0.0_r8
       l2err(ie,3) = get_lerror(n,ret(:,2),ret(:,3),'L2')
     endif
     !
     if (do_solve) then
       deallocate(ret,ret0)
     else
       deallocate(mat,imat)
     endif
     deallocate(a,b,c)
     if (.not.tri) then
       deallocate(d,e)
     endif
     !
   enddo ! do ie = 1,nexp
!
   call write_performance(tri,nexp,ns,nmet,string,l2err,elaps)
!
   deallocate(ns)
   deallocate(l2err,elaps,elaps_t)
   deallocate(string)
!
#if 0
     iden = matrix_multiplicity(n,d,imat)
     do j = 1,n
     do i = 1,n
       if (i.ne.j.and.iden(i,j).gt.1e-8) then
       print *,i,j,iden(i,j)
       endif
     enddo
     enddo
!
     print *,'(info) tridiagonal matrix: inverse algorithm (start)'
     stime = mpi_wtime(err)
     call tridiagonal_inverse(n,a,b,c,imat)
     etime = mpi_wtime(err)
     lerr  = matrix_chk_identity(n,mat,imat)
     if (do_prt) call print_matrix('imat = ',n,n,imat)
     call print_performance(lerr,etime-stime)
     print *,'(info) tridiagonal matrix: inverse algorithm (done)'
     print *,''
#endif
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_namelist(nlfilename)
!-------------------------------------------------------------------------------
   implicit none
   character(len=*), optional, intent(in   ) :: nlfilename
! local
   integer(i4) :: i
   namelist/inputs/do_solve, tri, nn, ntry, do_prt
!
   nn = -1
   if (present(nlfilename)) then
     open(21,file=trim(nlfilename),status='old')
     read(21,nml=inputs)
     close(21)
   else
     read(*,nml = inputs)
   endif
!
   do i = 1,max_n
     if (nn(i)<0) then
       nexp = i-1
       exit
     endif
   enddo
   allocate(ns(nexp))
   ns(:) = nn(1:nexp)
!
   print*,'--------------------------------'
   if (do_solve) then
     print*,' linear equation'
   else
     print*,' inverse of matrix'
   endif
   if (tri) then
     print*,' use tridiagonal matrix'
   else
     print*,' use pentadiagonal matrix'
   endif
   print*,' nexp   = ',nexp
   print*,' n      = ',ns
   print*,' ntry   = ',ntry
   print*,' do_prt = ',do_prt
   print*,'--------------------------------'
   print*,' '
!
   return
   end subroutine read_namelist
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_performance(err,elap)
!-------------------------------------------------------------------------------
   implicit none
   real(r8), intent(in   ) :: err, elap
! local
   integer(i4) :: i
!
   !print '(a,f10.4,a,f10.4)','   - L2 error = ',err,',  elapse = ',elap
!   print *,'   - L2 error = ',err,',  elapse = ',elap
!
   return
   end subroutine print_performance
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_performance(tr,ne,nd,nm,mnames,l2e,ela)
!-------------------------------------------------------------------------------
   implicit none
   logical(l4),                     intent(in   ) :: tr
   integer(i4),                     intent(in   ) :: ne
   integer(i4), dimension(ne),      intent(in   ) :: nd
   integer(i4),                     intent(in   ) :: nm
   character(len=*), dimension(nm), intent(in   ) :: mnames
   real(r8)   , dimension(ne,nm),   intent(in   ) :: l2e, ela
! local
   integer(i4) :: ie
   integer(i4) :: fid, neid, nmid, v_nid, v_eid, v_tid, err
!
   print *,'========== L2 error ==========='
   if     (nm==1) then
     print *,'                 ',trim(mnames(1))
     do ie = 1,ne
       print '(3x,i6.6,a,e15.5)', nd(ie),':',l2e(ie,1)
     enddo
     print *,'========== Elapse time ==========='
     print *,'                 ',trim(mnames(1))
     do ie = 1,ne
       print '(3x,i6.6,a,e15.5)', nd(ie),':',ela(ie,1)
     enddo
   elseif (nm==2) then
     print *,'                 ',trim(mnames(1)),'    ',trim(mnames(2))
     do ie = 1,ne
       print '(3x,i6.6,a,e15.5,x,e15.5)', nd(ie),':',l2e(ie,1),l2e(ie,2)
     enddo
     print *,'========== Elapse time ==========='
     print *,'                 ',trim(mnames(1)),'    ',trim(mnames(2))
     do ie = 1,ne
       print '(3x,i6.6,a,e15.5,x,e15.5)', nd(ie),':',ela(ie,1),ela(ie,2)
     enddo
   elseif (nm==3) then
     print *,'                 ',trim(mnames(1)),'    ',trim(mnames(2)),'    ',trim(mnames(3))
     do ie = 1,ne
       print '(3x,i6.6,a,e15.5,x,e15.5,x,e15.5)', nd(ie),':',l2e(ie,1),l2e(ie,2),l2e(ie,3)
     enddo
     print *,'========== Elapse time ==========='
     print *,'                 ',trim(mnames(1)),'    ',trim(mnames(2)),'    ',trim(mnames(3))
     do ie = 1,ne
       print '(3x,i6.6,a,e15.5,x,e15.5,x,e15.5)', nd(ie),':',ela(ie,1),ela(ie,2),ela(ie,3)
     enddo
   endif
!
   err = nf90_create('result.nc',ior(nf90_clobber,nf90_netcdf4),fid)
   err = nf90_def_dim(fid,'nexp',ne,neid)
   err = nf90_def_dim(fid,'nmethod',nm,nmid)
   err = nf90_def_var(fid,'nsize',nf90_int,neid,v_nid)
   err = nf90_def_var(fid,'l2err',nf90_double,(/neid,nmid/),v_eid)
   err = nf90_def_var(fid,'elapse',nf90_double,(/neid,nmid/),v_tid)
   if (tr) then
     err = nf90_put_att(fid,nf90_global,'tridiagonal',1)
   else
     err = nf90_put_att(fid,nf90_global,'tridiagonal',0)
   endif
   err = nf90_enddef(fid)
   err = nf90_put_var(fid,v_nid,nd)
   err = nf90_put_var(fid,v_eid,l2e)
   err = nf90_put_var(fid,v_tid,ela)
   err = nf90_close(fid)
!
   return
   end subroutine write_performance
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program test
!-------------------------------------------------------------------------------
