Module  misc  !

use iso_fortran_env, only: real64, int8
 
Implicit None

  contains
      
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    function inv(A) result(Ainv)
      integer, parameter:: dp = real64
      real(dp), dimension(:,:), intent(in):: A
      real(dp), dimension(size(A, 1), size(A, 2)):: Ainv
    
      real(dp), dimension(size(A, 1)):: work  ! work array for LAPACK
      integer, dimension(size(A, 1)):: ipiv   ! pivot indices
      integer:: n, info
    
      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI
    
      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A, 1)
    
      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)
    
      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if
    
      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)
    
      if (info /= 0) then
         stop 'Matrix inversion failed!'
      end if
    end function inv

! *-----------------------------------------------------------------------
! * gd-Gradient descent
! * G must be set correctly at the initial point X.
! *
! *___Name______Type________In/Out___Description_________________________
! *   N         Integer     In       Number of Variables.
! *   X(N)      Double      Both     Variables
! *   G(N)      Double      Both     Gradient
! *   TOL       Double      In       Relative convergence tolerance
! *   IFLAG     Integer     Out      Reverse Communication Flag
! *                                    on output:  0  done
! *                                                1  compute G and call again
! *-----------------------------------------------------------------------
      SUBROUTINE GD (N, chi, eta, G, TOL, IFLAG)
       IMPLICIT NONE
       INTEGER N, IFLAG
       DOUBLE PRECISION chi(N), eta(N),  G(2*N), TOL
       DOUBLE PRECISION LR
       PARAMETER (LR = 0.3D0)     ! Learning rate
       INTEGER I
       DOUBLE PRECISION GNORM      ! norm of gradient

       GNORM = 0.D0                ! convergence test
       DO I = 1, 2*N
         GNORM = GNORM+G(I)**2
       END DO
       GNORM = SQRT(GNORM)
       IF (GNORM < TOL) THEN
         IFLAG = 0
         RETURN                     ! success
       END IF

       DO I = 1, N                 ! take step
         chi(I) = chi(I) - LR*G(I)
         eta(I) = eta(I) - LR*G(N+I)
       END DO
       IFLAG = 1
       RETURN                      ! let main program evaluate G
      END  ! of gd

End Module  misc 
