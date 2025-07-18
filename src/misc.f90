!> @brief A module containing miscellaneous utility subroutines.
!> @details This module provides general-purpose functions used elsewhere in the application,
!>          such as matrix inversion and a simple gradient descent optimizer.
Module  misc  !

use iso_fortran_env, only: real64, int8
use ieee_arithmetic, only: ieee_is_finite

Implicit None

  contains

    !> @brief Calculates the inverse of a square matrix.
    !> @details This function is a wrapper around the LAPACK routines DGETRF (LU decomposition)
    !>          and DGETRI (matrix inversion from LU decomposition). It is a general-purpose
    !>          matrix inverter.
    !> @param[in] A The square matrix to be inverted.
    !> @return Ainv The inverse of matrix A.
    !> @note This function will stop the program if the matrix is singular or if the
    !>       inversion fails for any reason.
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
! * @brief Simple Gradient Descent (GD) optimizer.
! * @details This subroutine performs a single step of gradient descent.
! *          It is a basic optimizer provided for testing or simple cases.
! *          The main optimizer used in the application is the conjugate gradient
! *          method in `optimizer_module`.
! *
! * @param[in] N Number of variables (per parameter type, e.g., number of atom types).
! * @param[in] chi_in Input chi parameters.
! * @param[in] eta_in Input eta parameters.
! * @param[in] G The gradient vector (concatenated for chi and eta).
! * @param[in] TOL Relative convergence tolerance.
! * @param[out] IFLAG Reverse communication flag:
! *                   0 = done (convergence reached)
! *                   1 = compute G and call again
! *                  -1 = error (NaN/Inf encountered)
! * @param[out] chi_out Updated chi parameters.
! * @param[out] eta_out Updated eta parameters.
! *-----------------------------------------------------------------------
      SUBROUTINE GD (N, chi_in, eta_in, G, TOL, IFLAG, chi_out, eta_out)
       IMPLICIT NONE
       INTEGER N, IFLAG
       DOUBLE PRECISION chi_in(N), eta_in(N),  G(2*N), TOL
       DOUBLE PRECISION chi_out(N), eta_out(N)
       DOUBLE PRECISION LR
       PARAMETER (LR = 1.0D-12)     ! Learning rate
       INTEGER I
       DOUBLE PRECISION GNORM      ! norm of gradient

       print *, "GD: Before update"
       print *, "GD: chi = ", chi_in
       print *, "GD: eta = ", eta_in
       print *, "GD: G = ", G

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
         chi_out(I) = chi_in(I) - LR*G(I)
         eta_out(I) = eta_in(I) - LR*G(N+I)

         ! Check for NaN/Inf after update
         if (.not. ieee_is_finite(chi_out(I)) .or. .not. ieee_is_finite(eta_out(I))) then
             print *, "GD: Parameter became NaN/Inf at index ", I
             IFLAG = -1 ! Indicate error
             RETURN
         endif
       END DO

       print *, "GD: After update"
       print *, "GD: chi = ", chi_in
       print *, "GD: eta = ", eta_in

       IFLAG = 1
       RETURN                      ! let main program evaluate G
      END  ! of gd

End Module  misc