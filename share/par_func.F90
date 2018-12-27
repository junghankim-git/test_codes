!-------------------------------------------------------------------------------
   module par_func
!-------------------------------------------------------------------------------
!
!  abstract :
!
!  history log :
!    2015-10-02  junghan kim   First written
!    2016-01-05  junghan kim   Add the GlobalSUM, Deompose1D
!    2016-02-22  junghan kim   MPI_OP_FREE, 3D, 4D options
!
!-------------------------------------------------------------------------------
   use parallel, only : i4 => i4_par, r4 => r4_par, r8 => r8_par, kim_par,     &
                        par_min,   par_max
   use par_mpi , only : par_allreduce
!
   implicit none
!
   include 'mpif.h'
!
   private
!
   public :: globalmin, globalmax, globalsum, decompose1d
!
   interface ddpdd
     module procedure ddpdd_r4
     module procedure ddpdd
   end interface ddpdd
   interface globalmin
     module procedure globalmin_r4_s
     module procedure globalmin_r4_1d
     module procedure globalmin_r4_2d
     module procedure globalmin_r4_3d
     module procedure globalmin_r4_4d
     module procedure globalmin_r8_s
     module procedure globalmin_r8_1d
     module procedure globalmin_r8_2d
     module procedure globalmin_r8_3d
     module procedure globalmin_r8_4d
   end interface globalmin

   interface globalmax
     module procedure globalmax_r4_s
     module procedure globalmax_r4_1d
     module procedure globalmax_r4_2d
     module procedure globalmax_r4_3d
     module procedure globalmax_r4_4d
     module procedure globalmax_r8_s
     module procedure globalmax_r8_1d
     module procedure globalmax_r8_2d
     module procedure globalmax_r8_3d
     module procedure globalmax_r8_4d
   end interface globalmax

   interface globalsum
     module procedure globalsum_r4_s
     module procedure globalsum_r4_1d
     module procedure globalsum_r4_2d
     module procedure globalsum_r4_3d
     module procedure globalsum_r4_4d
     module procedure globalsum_r8_s
     module procedure globalsum_r8_1d
     module procedure globalsum_r8_2d
     module procedure globalsum_r8_3d
     module procedure globalsum_r8_4d
   end interface globalsum
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmin_r4_s(var) result(gmin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), intent(in) :: var
   real(r4)             :: gmin
   integer(i4)          :: ierr
!-------------------------------------------------------------------------------
   call par_allreduce(var, gmin, 1, par_op = par_min)

   end function globalmin_r4_s
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmin_r4_1d(var) result(gmin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:), intent(in) :: var
   real(r4)                           :: gmin
   integer(i4)                        :: ierr
   real(r4)                           :: lmin
!-------------------------------------------------------------------------------
   lmin = minval(var)

   call par_allreduce(lmin, gmin, 1, par_op = par_min)

   end function globalmin_r4_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmin_r4_2d(var) result(gmin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:, :), intent(in) :: var
   real(r4)                              :: gmin
   integer(i4)                           :: ierr
   real(r4)                              :: lmin
!-------------------------------------------------------------------------------
   lmin = minval(var)

   call par_allreduce(lmin, gmin, 1, par_op = par_min)

   end function globalmin_r4_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmin_r4_3d(var) result(gmin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:, :, :), intent(in) :: var
   real(r4)                                 :: gmin
   integer(i4)                              :: ierr
   real(r4)                                 :: lmin
!-------------------------------------------------------------------------------
   lmin = minval(var)

   call par_allreduce(lmin, gmin, 1, par_op = par_min)

   end function globalmin_r4_3d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmin_r4_4d(var) result(gmin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:, :, :, :), intent(in) :: var
   real(r4)                                    :: gmin
   integer(i4)                                 :: ierr
   real(r4)                                    :: lmin
!-------------------------------------------------------------------------------
   lmin = minval(var)

   call par_allreduce(lmin, gmin, 1, par_op = par_min)

   end function globalmin_r4_4d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmin_r8_s(var) result(gmin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in) :: var
   real(r8)             :: gmin
   integer(i4)          :: ierr
!-------------------------------------------------------------------------------
   call par_allreduce(var, gmin, 1, par_op = par_min)

   end function globalmin_r8_s
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmin_r8_1d(var) result(gmin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:), intent(in) :: var
   real(r8)                           :: gmin
   integer(i4)                        :: ierr
   real(r8)                           :: lmin
!-------------------------------------------------------------------------------
   lmin = minval(var)

   call par_allreduce(lmin, gmin, 1, par_op = par_min)

   end function globalmin_r8_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmin_r8_2d(var) result(gmin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:, :), intent(in) :: var
   real(r8)                              :: gmin
   integer(i4)                           :: ierr
   real(r8)                              :: lmin
!-------------------------------------------------------------------------------
   lmin = minval(var)

   call par_allreduce(lmin, gmin, 1, par_op = par_min)

   end function globalmin_r8_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmin_r8_3d(var) result(gmin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:, :, :), intent(in) :: var
   real(r8)                                 :: gmin
   integer(i4)                              :: ierr
   real(r8)                                 :: lmin
!-------------------------------------------------------------------------------
   lmin = minval(var)

   call par_allreduce(lmin, gmin, 1, par_op = par_min)

   end function globalmin_r8_3d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmin_r8_4d(var) result(gmin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:, :, :, :), intent(in) :: var
   real(r8)                                    :: gmin
   integer(i4)                                 :: ierr
   real(r8)                                    :: lmin
!-------------------------------------------------------------------------------
   lmin = minval(var)

   call par_allreduce(lmin, gmin, 1, par_op = par_min)

   end function globalmin_r8_4d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmax_r4_s(var) result(gmax)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), intent(in) :: var
   real(r4)             :: gmax
   integer(i4)          :: ierr
!-------------------------------------------------------------------------------
   call par_allreduce(var, gmax, 1, par_op = par_max)

   end function globalmax_r4_s
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmax_r4_1d(var) result(gmax)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:), intent(in) :: var
   real(r4)                           :: gmax
   integer(i4)                        :: ierr
   real(r4)                           :: lmax
!-------------------------------------------------------------------------------
   lmax = maxval(var)

   call par_allreduce(lmax, gmax, 1, par_op = par_max)

   end function globalmax_r4_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmax_r4_2d(var) result(gmax)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:, :), intent(in) :: var
   real(r4)                              :: gmax
   integer(i4)                           :: ierr
   real(r4)                              :: lmax
!-------------------------------------------------------------------------------
   lmax = maxval(var)

   call par_allreduce(lmax, gmax, 1, par_op = par_max)

   end function globalmax_r4_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmax_r4_3d(var) result(gmax)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:, :, :), intent(in) :: var
   real(r4)                                 :: gmax
   integer(i4)                              :: ierr
   real(r4)                                 :: lmax
!-------------------------------------------------------------------------------
   lmax = maxval(var)

   call par_allreduce(lmax, gmax, 1, par_op = par_max)

   end function globalmax_r4_3d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmax_r4_4d(var) result(gmax)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:, :, :, :), intent(in) :: var
   real(r4)                                    :: gmax
   integer(i4)                                 :: ierr
   real(r4)                                    :: lmax
!-------------------------------------------------------------------------------
   lmax = maxval(var)

   call par_allreduce(lmax, gmax, 1, par_op = par_max)

   end function globalmax_r4_4d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmax_r8_s(var) result(gmax)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in) :: var
   real(r8)             :: gmax
   integer(i4)          :: ierr
!-------------------------------------------------------------------------------
   call par_allreduce(var, gmax, 1, par_op = par_max)

   end function globalmax_r8_s
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmax_r8_1d(var) result(gmax)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:), intent(in) :: var
   real(r8)                           :: gmax
   integer(i4)                        :: ierr
   real(r8)                           :: lmax
!-------------------------------------------------------------------------------
   lmax = maxval(var)

   call par_allreduce(lmax, gmax, 1, par_op = par_max)

   end function globalmax_r8_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmax_r8_2d(var) result(gmax)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:, :), intent(in) :: var
   real(r8)                              :: gmax
   integer(i4)                           :: ierr
   real(r8)                              :: lmax
!-------------------------------------------------------------------------------
   lmax = maxval(var)

   call par_allreduce(lmax, gmax, 1, par_op = par_max)

   end function globalmax_r8_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmax_r8_3d(var) result(gmax)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:, :, :), intent(in) :: var
   real(r8)                                 :: gmax
   integer(i4)                              :: ierr
   real(r8)                                 :: lmax
!-------------------------------------------------------------------------------
   lmax = maxval(var)

   call par_allreduce(lmax, gmax, 1, par_op = par_max)

   end function globalmax_r8_3d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalmax_r8_4d(var) result(gmax)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:, :, :, :), intent(in) :: var
   real(r8)                                    :: gmax
   integer(i4)                                 :: ierr
   real(r8)                                    :: lmax
!-------------------------------------------------------------------------------
   lmax = maxval(var)

   call par_allreduce(lmax, gmax, 1, par_op = par_max)

   end function globalmax_r8_4d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalsum_r4_s(var) result(gsum)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), intent(in) :: var
   real(r4)             :: gsum
   complex(r4)          :: lsum_c, gsum_c
   integer(i4)          :: mpi_sumop, ierr
!-------------------------------------------------------------------------------
   lsum_c = cmplx(var, 0.0_r4, r4)

   call mpi_op_create(ddpdd, .true., mpi_sumop, ierr)
   call mpi_allreduce(lsum_c, gsum_c, 1, mpi_complex16, mpi_sumop,             &
                      kim_par%comm, ierr)

   gsum = real(gsum_c)
   call mpi_op_free(mpi_sumop, ierr)

   end function globalsum_r4_s
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalsum_r4_1d(var) result(gsum)
!-------------------------------------------------------------------------------
   implicit none
   real(r4), dimension(:), intent(in) :: var
   real(r4)                           :: gsum
   complex(r4)                        :: tmp, lsum_c, gsum_c
   integer(i4)                        :: mpi_sumop, ierr
   integer(i4)                        :: i, n
!-------------------------------------------------------------------------------
   n      = size(var)
   lsum_c = cmplx(0.0_r4, 0.0_r4, r4)
!
   do i = 1,n
     tmp = cmplx(var(i), 0.0_r4, r4)
     call ddpdd(tmp, lsum_c)
   enddo
!
   call mpi_op_create(ddpdd, .true., mpi_sumop, ierr)
   call mpi_allreduce(lsum_c, gsum_c, 1, mpi_complex16, mpi_sumop,             &
                      kim_par%comm, ierr)

   gsum = real(gsum_c)
   call mpi_op_free(mpi_sumop, ierr)

   end function globalsum_r4_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalsum_r4_2d(var) result(gsum)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:, :), intent(in) :: var
   real(r4)                              :: gsum
   complex(r4)                           :: tmp, lsum_c, gsum_c
   integer(i4)                           :: mpi_sumop, ierr
   integer(i4)                           :: i, j, n(2)
!-------------------------------------------------------------------------------
   n      = shape(var)
   lsum_c = cmplx(0.0_r4, 0.0_r4, r4)
!
   do j = 1,n(2)
     do i = 1,n(1)
       tmp = cmplx(var(i, j), 0.0_r4, r4)
       call ddpdd(tmp, lsum_c)
     enddo
   enddo
!
   call mpi_op_create(ddpdd, .true., mpi_sumop, ierr)
   call mpi_allreduce(lsum_c, gsum_c, 1, mpi_complex16, mpi_sumop,             &
                      kim_par%comm, ierr)

   gsum = real(gsum_c)
   call mpi_op_free(mpi_sumop, ierr)

   end function globalsum_r4_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalsum_r4_3d(var) result(gsum)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:, :, :), intent(in) :: var
   real(r4)                                 :: gsum
   complex(r4)                              :: tmp, lsum_c, gsum_c
   integer(i4)                              :: mpi_sumop, ierr
   integer(i4)                              :: i, j, k, n(3)
!-------------------------------------------------------------------------------
   n      = shape(var)
   lsum_c = cmplx(0.0_r4, 0.0_r4, r4)
!
   do k = 1,n(3)
     do j = 1,n(2)
       do i = 1,n(1)
         tmp = cmplx(var(i, j, k), 0.0_r4, r4)
         call ddpdd(tmp, lsum_c)
       enddo
     enddo
   enddo
!
   call mpi_op_create(ddpdd, .true., mpi_sumop, ierr)
   call mpi_allreduce(lsum_c, gsum_c, 1, mpi_complex16, mpi_sumop,             &
                      kim_par%comm, ierr)

   gsum = real(gsum_c)
   call mpi_op_free(mpi_sumop, ierr)

   end function globalsum_r4_3d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalsum_r4_4d(var) result(gsum)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:, :, :, :), intent(in) :: var
   real(r4)                                    :: gsum
   complex(r4)                                 :: tmp, lsum_c, gsum_c
   integer(i4)                                 :: mpi_sumop, ierr
   integer(i4)                                 :: i, j, k, l, n(4)
!-------------------------------------------------------------------------------
   n      = shape(var)
   lsum_c = cmplx(0.0_r4, 0.0_r4, r4)
!
   do l = 1,n(4)
     do k = 1,n(3)
       do j = 1,n(2)
         do i = 1,n(1)
           tmp = cmplx(var(i, j, k, l), 0.0_r4, r4)
           call ddpdd(tmp, lsum_c)
         enddo
       enddo
     enddo
   enddo
!
   call mpi_op_create(ddpdd, .true., mpi_sumop, ierr)
   call mpi_allreduce(lsum_c, gsum_c, 1, mpi_complex16, mpi_sumop,             &
                      kim_par%comm, ierr)

   gsum = real(gsum_c)
   call mpi_op_free(mpi_sumop, ierr)

   end function globalsum_r4_4d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalsum_r8_s(var) result(gsum)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in) :: var
   real(r8)             :: gsum
   complex(r8)          :: lsum_c, gsum_c
   integer(i4)          :: mpi_sumop, ierr
!-------------------------------------------------------------------------------
   lsum_c = cmplx(var, 0.0_r8, r8)

   call mpi_op_create(ddpdd, .true., mpi_sumop, ierr)
   call mpi_allreduce(lsum_c, gsum_c, 1, mpi_complex16, mpi_sumop,             &
                      kim_par%comm, ierr)

   gsum = real(gsum_c)
   call mpi_op_free(mpi_sumop, ierr)

   end function globalsum_r8_s
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalsum_r8_1d(var) result(gsum)
!-------------------------------------------------------------------------------
   implicit none
   real(r8), dimension(:), intent(in) :: var
   real(r8)                           :: gsum
   complex(r8)                        :: tmp, lsum_c, gsum_c
   integer(i4)                        :: mpi_sumop, ierr
   integer(i4)                        :: i, n
!-------------------------------------------------------------------------------
   n      = size(var)
   lsum_c = cmplx(0.0_r8, 0.0_r8, r8)
!
   do i = 1,n
     tmp = cmplx(var(i), 0.0_r8, r8)
     call ddpdd(tmp, lsum_c)
   enddo
!
   call mpi_op_create(ddpdd, .true., mpi_sumop, ierr)
   call mpi_allreduce(lsum_c, gsum_c, 1, mpi_complex16, mpi_sumop,             &
                      kim_par%comm, ierr)

   gsum = real(gsum_c)
   call mpi_op_free(mpi_sumop, ierr)

   end function globalsum_r8_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalsum_r8_2d(var) result(gsum)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:, :), intent(in) :: var
   real(r8)                              :: gsum
   complex(r8)                           :: tmp, lsum_c, gsum_c
   integer(i4)                           :: mpi_sumop, ierr
   integer(i4)                           :: i, j, n(2)
!-------------------------------------------------------------------------------
   n      = shape(var)
   lsum_c = cmplx(0.0_r8, 0.0_r8, r8)
!
   do j = 1,n(2)
     do i = 1,n(1)
       tmp = cmplx(var(i, j), 0.0_r8, r8)
       call ddpdd(tmp, lsum_c)
     enddo
   enddo
!
   call mpi_op_create(ddpdd, .true., mpi_sumop, ierr)
   call mpi_allreduce(lsum_c, gsum_c, 1, mpi_complex16, mpi_sumop,             &
                      kim_par%comm, ierr)

   gsum = real(gsum_c)
   call mpi_op_free(mpi_sumop, ierr)

   end function globalsum_r8_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalsum_r8_3d(var) result(gsum)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:, :, :), intent(in) :: var
   real(r8)                                 :: gsum
   complex(r8)                              :: tmp, lsum_c, gsum_c
   integer(i4)                              :: mpi_sumop, ierr
   integer(i4)                              :: i, j, k, n(3)
!-------------------------------------------------------------------------------
   n      = shape(var)
   lsum_c = cmplx(0.0_r8, 0.0_r8, r8)
!
   do k = 1,n(3)
     do j = 1,n(2)
       do i = 1,n(1)
         tmp = cmplx(var(i, j, k), 0.0_r8, r8)
         call ddpdd(tmp, lsum_c)
       enddo
     enddo
   enddo
!
   call mpi_op_create(ddpdd, .true., mpi_sumop, ierr)
   call mpi_allreduce(lsum_c, gsum_c, 1, mpi_complex16, mpi_sumop,             &
                      kim_par%comm, ierr)

   gsum = real(gsum_c)
   call mpi_op_free(mpi_sumop, ierr)

   end function globalsum_r8_3d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function globalsum_r8_4d(var) result(gsum)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:, :, :, :), intent(in) :: var
   real(r8)                                    :: gsum
   complex(r8)                                 :: tmp, lsum_c, gsum_c
   integer(i4)                                 :: mpi_sumop, ierr
   integer(i4)                                 :: i, j, k, l, n(4)
!-------------------------------------------------------------------------------
   n      = shape(var)
   lsum_c = cmplx(0.0_r8, 0.0_r8, r8)
!
   do l = 1,n(4)
     do k = 1,n(3)
       do j = 1,n(2)
         do i = 1,n(1)
           tmp = cmplx(var(i, j, k, l), 0.0_r8, r8)
           call ddpdd(tmp, lsum_c)
         enddo
       enddo
     enddo
   enddo
!
   call mpi_op_create(ddpdd, .true., mpi_sumop, ierr)
   call mpi_allreduce(lsum_c, gsum_c, 1, mpi_complex16, mpi_sumop,             &
                      kim_par%comm, ierr)

   gsum = real(gsum_c)
   call mpi_op_free(mpi_sumop, ierr)

   end function globalsum_r8_4d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ddpdd_r4(dda, ddb)
!-------------------------------------------------------------------------------
   implicit none
!
   complex(r4), intent(in   ) :: dda
   complex(r4), intent(inout) :: ddb
   real(r4)                   :: e, t1, t2
!-------------------------------------------------------------------------------
   t1 = real(dda) + real(ddb)
   e = t1 - real(dda)
   t2 =((real(ddb) - e) +(real(dda) -(t1 - e))) + imag(dda) + imag(ddb)
   ddb = cmplx(t1 + t2, t2 -((t1 + t2) - t1), r4)

   return

   end subroutine ddpdd_r4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ddpdd(dda, ddb)
!-------------------------------------------------------------------------------
   implicit none
!
   complex(r8), intent(in   ) :: dda
   complex(r8), intent(inout) :: ddb
   real(r8)                   :: e, t1, t2
!-------------------------------------------------------------------------------
   t1 = real(dda) + real(ddb)
   e = t1 - real(dda)
   t2 =((real(ddb) - e) +(real(dda) -(t1 - e))) + imag(dda) + imag(ddb)
   ddb = cmplx(t1 + t2, t2 -((t1 + t2) - t1), r8)

   return

   end subroutine ddpdd
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function decompose1d(n1, n2, nprocs, rank, ista, iend) result(nd)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: n1, n2
   integer(i4),           intent(in   ) :: nprocs, rank
   integer(i4), optional, intent(  out) :: ista, iend
   integer(i4) :: nd
   integer(i4) :: l_sta, l_end, domain, extra
!-------------------------------------------------------------------------------
   domain = n2 - n1 + 1
   nd = domain / nprocs
   extra = mod(domain, nprocs)
   l_sta = nd * rank + n1 + min(rank, extra)
   l_end = l_sta + nd - 1
!
   if (rank.le.extra-1) then
     l_end = l_end + 1
   endif
!
   nd = l_end - l_sta + 1
!
   if (present(ista)) then
     ista = l_sta
   endif
!
   if (present(iend)) then
     iend = l_end
   endif
!
   return

   end function decompose1d
!-------------------------------------------------------------------------------
   end module par_func
!-------------------------------------------------------------------------------
