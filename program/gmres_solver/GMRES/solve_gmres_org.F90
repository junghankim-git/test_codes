!----------------------------------------------------------------------!
!This subroutine solves a linear problem using GMRES.
!Written by Paul Fischer and
!Modified by Francis X. Giraldo on 7/08
!           Department of Applied Mathematics
!           Naval Postgraduate School
!           Monterey, CA 93943-5216
!----------------------------------------------------------------------!
subroutine solve_gmres_1d(A,b,x)

  implicit none
 
  nz
  solver_tol = 1e-12 
  ndof = nz * nvar 
  maxiter = 500

  real krylov(ndof,maxiter)
  real krylow(ndof)

  !Solver arrays
  real x(ndof), b(ndof)

  !GMRES arrays
  real h(maxiter,maxiter), s(maxiter), c(maxiter), gamma0(maxiter)
  real h_loc, h_glob, twonorm, glsc2
  logical ifdone
  integer i, j, k, n, maxiter
  real one, eps, t, etol, etol2, den

  !>sh.kang @14.1.28
  integer preconType


  !set machine tolerances
  n=ndof
  one = 1.
  eps = 1.e-20
  if (one+eps == one) eps = 1.e-14
  if (one+eps == one) eps = 1.e-7
  if (one+eps == one) eps = 1.e-16
  !      eps=1.0e-7
  eps=solver_tol
  
  !Initialize
  x=0

  ! sh.kang @14.1.28
  preconType = 0
  if (preconType .eq. 0) then
    krylow(:) = b(:)

  else
    !Solve M*w=b
    call precon(krylow,b,adiag,n)
  
  end if
  
  !Compute L2 Norm
  gamma0(1) = twonorm(krylow,n)
  
!  if (gamma0(1) == 0) then 
!       if (lprint_diagnostics .and. ix == ncol) then
!          kiter0=nint( real(kiter-kiter0)/ncol )
!          write(*,'("     istage k Tol     = ",i2,i5,e16.8)')istage,kiter0,gamma0(k+1)/gamma1
!          kiter0=kiter
!       end if
!       return
!  end if

  t = 1./gamma0(1)
  call cmult2(krylov(:,1),krylow,t,n)
  
!  if (gamma1 < 0) then
     etol  = gamma0(1)*eps
     etol2 = gamma0(1)*eps*eps
!  else    
!     etol  = gamma1*eps     
!     etol2 = gamma1*eps*eps
!  end if
  
  !Begin Arnoldi w/ modified gram-schmidt
  do k=1,maxiter
     !call create_lhs_gmres_no_schur_1d(krylov(:,k+1),krylov(:,k),lambda,ndof,n_r,ix)
     krylov(:,k+1) = A*krylov(:,k)

     !sh.kang @14.1.28
     if (preconType .eq. 0) then
       krylow(:) = krylov(:,k+1)
     
     else 
       !Solve M*w=v(:,k+1)
       call precon(krylow,krylov(:,k+1),adiag,n)

     end if
     
     do j=1,k
        h(j,k) = glsc2(krylow,krylov(:,j),n)
        t      = -h(j,k)
        call add2s2(krylow,krylov(:,j),t,n)
     end do !j
     
     !Compute L2 Norm
     h(k+1,k) = twonorm(krylow,n)
     
     if (abs(h(k+1,k)) <  etol2) then
        ifdone=.true.
     else
        ifdone=.false.
        t = 1./h(k+1,k)
        call cmult2(krylov(:,k+1),krylow,t,n)
     end if
     
     !apply Given's rotations to new column of H
     do i=1,k-1
        t = h(i,k)
        h(i  ,k) =  c(i)*t + s(i)*h(i+1,k)
        h(i+1,k) = -s(i)*t + c(i)*h(i+1,k)
     end do !i
     den        =  sqrt(  h(k,k)*h(k,k) + h(k+1,k)*h(k+1,k)  )
     c(k)       =  h(k  ,k) / den
     s(k)       =  h(k+1,k) / den
     h(k,k)     =  c(k)*h(k,k)+s(k)*h(k+1,k)
     gamma0(k+1) = -s(k)*gamma0(k)
     gamma0(k  ) =  c(k)*gamma0(k)
     if (ifdone .or. abs(gamma0(k+1)) < etol) exit
  end do !k
  
  !if we're here, we exceeded the max iteration, reduce k by 1
  if (k > maxiter) then
     print*,' MAXITER exceeded'
     k = k-1
  end if
  
  !Compute solution via back substitution
 ! kiter=kiter + k
 ! if (lprint_diagnostics .and. ix == ncol) then
 !    kiter0=nint( real(kiter-kiter0)/ncol )
 !    write(*,'("     istage k Tol     = ",i2,i5,e16.8)')istage,kiter0,gamma0(k+1)/gamma1
 !    kiter0=kiter
 ! end if

  do i=k,1,-1
     t = gamma0(i)
     do j=k,i+1,-1
        t = t - h(i,j)*c(j)
     end do !j
     c(i) = t/h(i,i)
  end do !i
  
  !Sum up Arnoldi vectors
  do i=1,k
     call add2s2(x,krylov(:,i),c(i),n)
  end do !i
  
end subroutine solve_gmres_1d




!-----------------------------------------------------------------------
subroutine add2s2(a,b,c1,n)
  implicit none
  real a(n), b(n), c1
  integer i, n

  do i=1,n
     a(i)=a(i)+c1*b(i)
  end do
  
end subroutine add2s2
!-----------------------------------------------------------------------
function twonorm(x,n)
  implicit none
  real x(n), twonorm, tscal
  integer i, n

  tscal = x(1)*x(1)
  do i=2,n
     tscal = tscal+x(i)*x(i)
  end do
  twonorm=tscal
  if (twonorm.gt.0) twonorm = sqrt(twonorm)
  
end function twonorm
!-----------------------------------------------------------------------
function glsc2(x,y,n)
  implicit none
  real x(n), y(n), glsc2, tscal
  integer i, n
  
  tscal = x(1)*y(1)
  do i=2,n
     tscal = tscal+x(i)*y(i)
  end do
  glsc2=tscal
  
end function glsc2
!-----------------------------------------------------------------------
subroutine cmult2(a,b,c,n)
  implicit none
  real a(n), b(n), c
  integer i, n

  do i = 1, n
     a(i) = b(i)*c
  end do
  
end subroutine cmult2
!-----------------------------------------------------------------------
subroutine precon(y,x,a,n)
  implicit none
  real y(n), x(n), a(n)
  integer i, n

  do i=1,n
     y(i)=x(i)/a(i)
  end do !i
  
end subroutine precon
!-----------------------------------------------------------------------
function vlsc2(x,y,n)
  implicit none
  real x(n), y(n), vlsc2, tscal
  integer i, n
  
  tscal = x(1)*y(1)
  do i=2,n
     tscal = tscal+x(i)*y(i)
  end do
  vlsc2 = tscal
  
end function vlsc2
!-----------------------------------------------------------------------
