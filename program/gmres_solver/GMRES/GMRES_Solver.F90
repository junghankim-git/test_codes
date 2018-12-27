MODULE GMRES_Solver

  PRIVATE

  INTEGER, PARAMETER  :: i4 = 4, r8 = 8, l4 = 4

  INTEGER(i4), PARAMETER  :: maxiter = 1000
  REAL(r8), PARAMETER     :: solver_tol = 1.0D-280



  PUBLIC :: SolveGMRES!, twonorm, glsc2




CONTAINS



  SUBROUTINE SolveGMRES(n, A, b, x)
  
    IMPLICIT NONE
   
    INTEGER(i4), INTENT(IN)                  :: n
    REAL(r8), DIMENSION(n,n), INTENT(INOUT)  :: A
    REAL(r8), DIMENSION(n), INTENT(INOUT)    :: b
    REAL(r8), DIMENSION(n), INTENT(INOUT)    :: x


    ! local
    REAL(r8)    :: krylov(n,maxiter)
    REAL(r8)    :: krylow(n)
     !GMRES arrays
    REAL(r8)    :: h(maxiter,maxiter), s(maxiter), c(maxiter), gamma0(maxiter)
    REAL(r8)    :: h_loc, h_glob!, twonorm, glsc2
    LOGICAL(l4) :: ifdone
    REAL(r8)    :: eps, t, etol, etol2, den
    INTEGER(i4) :: i, j, k
    !>sh.kang @14.1.28
    INTEGER(i4) preconType
  
  
    !set machine tolerances
    eps=solver_tol
    
    !Initialize
    x(:) = 0.0_r8
  
    !Precondition
    krylow(:) = b(:)
    

    !Compute L2 Norm
    gamma0(1) = twonorm(krylow,n)
    
  
    t = 1./gamma0(1)
    CALL cmult2(krylov(:,1),krylow,t,n)
    
    etol  = gamma0(1)*eps
    etol2 = gamma0(1)*eps*eps
    
    !Begin Arnoldi w/ modified gram-schmidt
    DO k=1,maxiter
       !CALL create_lhs_gmres_no_schur_1d(krylov(:,k+1),krylov(:,k),lambda,n,n_r,ix)

       krylov(:,k+1) = 0.0_r8
       DO j = 1, n
        DO i = 1, n
          krylov(j,k+1) = krylov(j,k+1) + A(j,i)*krylov(i,k)
        ENDDO
       ENDDO
  
       krylow(:) = krylov(:,k+1)

       
       DO j=1,k
          h(j,k) = glsc2(krylow,krylov(:,j),n)
          t      = -h(j,k)
          CALL add2s2(krylow,krylov(:,j),t,n)
       END DO !j
       
       !Compute L2 Norm
       h(k+1,k) = twonorm(krylow,n)
       
       IF (abs(h(k+1,k)) <  etol2) then
          ifdone=.true.
       else
          ifdone=.false.
          t = 1./h(k+1,k)
          CALL cmult2(krylov(:,k+1),krylow,t,n)
       END if
       
       !apply Given's rotations to new column of H
       DO i=1,k-1
          t = h(i,k)
          h(i  ,k) =  c(i)*t + s(i)*h(i+1,k)
          h(i+1,k) = -s(i)*t + c(i)*h(i+1,k)
       END DO !i
       den        =  sqrt(  h(k,k)*h(k,k) + h(k+1,k)*h(k+1,k)  )
       c(k)       =  h(k  ,k) / den
       s(k)       =  h(k+1,k) / den
       h(k,k)     =  c(k)*h(k,k)+s(k)*h(k+1,k)
       gamma0(k+1) = -s(k)*gamma0(k)
       gamma0(k  ) =  c(k)*gamma0(k)
 
       IF (ifdone .or. (abs(gamma0(k+1)) < etol)) exit
    END DO !k
    
    !IF we're here, we exceeded the max iteration, reduce k by 1
    IF (k > maxiter) then
       print*,' MAXITER exceeded'
       k = k-1
    END if
    PRINT *, ' GMRES: # of interation = ', k
    
    !Compute solution via back substitution
   ! kiter=kiter + k
   ! IF (lprint_diagnostics .and. ix == ncol) then
   !    kiter0=nint( REAL(r8)(kiter-kiter0)/ncol )
   !    write(*,'("     istage k Tol     = ",i2,i5,e16.8)')istage,kiter0,gamma0(k+1)/gamma1
   !    kiter0=kiter
   ! END if
  
    DO i=k,1,-1
       t = gamma0(i)
       DO j=k,i+1,-1
          t = t - h(i,j)*c(j)
       END DO !j
       c(i) = t/h(i,i)
    END DO !i
    
    !Sum up Arnoldi vectors
    DO i=1,k
       CALL add2s2(x,krylov(:,i),c(i),n)
    END DO !i
    
  END SUBROUTINE SolveGMRES
  
  
  
  
  !-----------------------------------------------------------------------
  SUBROUTINE add2s2(a,b,c1,n)
    IMPLICIT NONE
    REAL(r8) a(n), b(n), c1
    INTEGER(i4) i, n
  
    DO i=1,n
       a(i)=a(i)+c1*b(i)
    END DO
    
  END SUBROUTINE add2s2




  !-----------------------------------------------------------------------
  FUNCTION twonorm(x,n)
    IMPLICIT NONE
    REAL(r8) x(n), twonorm, tscal
    INTEGER(i4) i, n
  
    tscal = x(1)*x(1)
    DO i=2,n
       tscal = tscal+x(i)*x(i)
    END DO
    twonorm=tscal
    IF (twonorm.gt.0) twonorm = sqrt(twonorm)
    
  END FUNCTION twonorm





  !-----------------------------------------------------------------------
  FUNCTION glsc2(x,y,n)
    IMPLICIT NONE
    REAL(r8) x(n), y(n), glsc2, tscal
    INTEGER(i4) i, n
    
    tscal = x(1)*y(1)
    DO i=2,n
       tscal = tscal+x(i)*y(i)
    END DO
    glsc2=tscal
    
  END FUNCTION glsc2
  !-----------------------------------------------------------------------
  SUBROUTINE cmult2(a,b,c,n)
    IMPLICIT NONE
    REAL(r8) a(n), b(n), c
    INTEGER(i4) i, n
  
    DO i = 1, n
       a(i) = b(i)*c
    END DO
    
  END SUBROUTINE cmult2
  !-----------------------------------------------------------------------
  SUBROUTINE precon(y,x,a,n)
    IMPLICIT NONE
    REAL(r8) y(n), x(n), a(n)
    INTEGER(i4) i, n
  
    DO i=1,n
       y(i)=x(i)/a(i)
    END DO !i
    
  END SUBROUTINE precon
  !-----------------------------------------------------------------------
  FUNCTION vlsc2(x,y,n)
    IMPLICIT NONE
    REAL(r8) x(n), y(n), vlsc2, tscal
    INTEGER(i4) i, n
    
    tscal = x(1)*y(1)
    DO i=2,n
       tscal = tscal+x(i)*y(i)
    END DO
    vlsc2 = tscal
    
  END FUNCTION vlsc2
  !-----------------------------------------------------------------------


END MODULE GMRES_Solver
