! This module contains the subroutine for obtaining the potential

MODULE potential

  ! We import the modules that are utilized by these subroutines/functions
  USE initialize
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  CONTAINS

  ! This calculates the potential using the Gauss-Seidel method
  SUBROUTINE get_potential(nx, ny, dx, dy, rho, phi)

    ! The tolerance level for the Gauss-Seidel solver is set here
    REAL(REAL64), PARAMETER :: tol=1e-5
    INTEGER(INT32), INTENT(IN) :: nx, ny
    INTEGER(INT32) :: i, j
    REAL(REAL64), INTENT(IN) :: dx, dy
    REAL(REAL64), DIMENSION(0:nx+1, 0:ny+1), INTENT(IN) :: rho

    ! phi is variable for the potential values
    ! phix, phiy, etot, drms, ratio will be used in its calculation
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: phi
    REAL(REAL64) :: phix, phiy, etot = 0.0_REAL64, drms = 0.0_REAL64
    REAL(REAL64) :: ratio = 1.0_REAL64

    ! Allocating phi to hold nx+2 by ny+2 values
    ALLOCATE(phi(0:nx+1, 0:ny+1))
    phi = 0.0_REAL64

    ! Iterations to find phi by Gauss-Seidel method
    DO WHILE(ratio>tol)

      ! This loop obtains an iteration of phi
      DO i=1,nx
        DO j=1,ny
          phi(i,j) = -(rho(i,j) - (phi(i+1,j)+phi(i-1,j))/dx**2 - &
          (phi(i,j+1)+phi(i,j-1))/dy**2)/(2.0_REAL64/dx**2 + 2.0_REAL64/dy**2)
        END DO
      END DO

      ! This loop and the code that follows calculates how far phi is from the
      ! true solution, the while loop continues so long as ratio > tol
      etot = 0.0_REAL64
      drms = 0.0_REAL64
      DO i=1,nx
        DO j=1,ny
          phix = (phi(i-1,j) - 2.0_REAL64*phi(i,j) + phi(i+1,j))/dx**2
          phiy = (phi(i,j-1) - 2.0_REAL64*phi(i,j) + phi(i,j+1))/dy**2
          etot = etot + ABS(phix + phiy - rho(i,j))
          drms = drms + (phix + phiy)**2
        END DO
      END DO

      drms = SQRT(drms/REAL(nx*ny, KIND=REAL64))

      ! Check if drms equals 0 and let ratio be etot if it does
      ! Otherwise ratio is the regular value
      IF (drms == 0) THEN
        ratio = etot
      ELSE
        ratio = etot/drms
      END IF

    END DO

  END SUBROUTINE get_potential

END MODULE potential
