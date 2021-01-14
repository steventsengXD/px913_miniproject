MODULE helpers

  ! We import the modules that are utilized by these subroutines/functions
  ! USE command_line
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  CONTAINS


  SUBROUTINE init_rho(x, y, nx, ny, dx, dy, rho)

    INTEGER(INT32) :: i=1, j=1
    INTEGER(INT32), INTENT(IN) :: nx, ny
    REAL(REAL64), INTENT(IN) :: dx, dy
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: x, y
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: rho

    ! ALLOCATE(rho(2*nx+1,2*ny+1))
    ALLOCATE(x(0:nx+2),y(0:ny+2))
    ALLOCATE(rho(0:nx+2,0:ny+2))

    x(1)=-1
    DO i=1,nx+1
      x(i+1) = x(i) + dx
    END DO

    y(1)=-1
    DO i=1,ny+1
      y(i+1) = y(i) + dy
    END DO

    DO i = 1,nx+1
      DO j = 1,ny+1
        rho(i,j) = EXP(-(x(i)/0.1_REAL64)**2 - (y(j)/0.1_REAL64)**2)
      END DO
    END DO

  DEALLOCATE(x, y)

  END SUBROUTINE init_rho


  SUBROUTINE gauss_seidel(nx, ny, dx, dy, rho, phi)

    REAL(REAL64) :: tol=1e-5
    INTEGER(INT32), INTENT(IN) :: nx, ny
    INTEGER(INT32) :: i = 1, j = 1
    REAL(REAL64), INTENT(IN) :: dx, dy
    REAL(REAL64), DIMENSION(0:nx+2,0:ny+2), INTENT(IN) :: rho
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: phi
    REAL(REAL64) :: phix, phiy, etot = 0.0_REAL64, drms = 0.0_REAL64
    REAL(REAL64) :: denom, ratio = 1.0_REAL64

    ALLOCATE(phi(0:nx+2,0:ny+2))

    ! Calculate denominator here first for neatness purposes
    denom = -2.0_REAL64*(1.0_REAL64/dx**2+1.0_REAL64/dy**2)

    DO WHILE(ratio>tol)

      DO i=1,nx+1
        DO j=1,ny+1
          phi(i,j) = (-(phi(i+1,j)+phi(i-1,j))/dx**2 - (phi(i,j+1)+phi(i,j-1))/dy**2 &
          + rho(i,j))/denom
        END DO
      END DO

      DO i=1,nx+1
        DO j=1,ny+1
          phi(i,j) = (-(phi(i+1,j)+phi(i-1,j))/dx**2 - (phi(i,j+1)+phi(i,j-1))/dy**2 &
          + rho(i,j))/denom
        END DO
      END DO

      etot = 0.0_REAL64
      DO i=1,nx+1
        DO j=1,ny+1
          phix = (phi(i-1,j)-2.0_REAL64*phi(i,j)+phi(i+1,j))/dx**2
          phiy = (phi(i,j-1)-2.0_REAL64*phi(i,j)+phi(i,j+1))/dy**2
          etot = etot + ABS(phix + phiy - rho(i,j))
        END DO
      END DO

      drms = 0.0_REAL64

      DO i=1,nx+1
        DO j=1,ny+1
          phix = (phi(i-1,j)-2.0_REAL64*phi(i,j)+phi(i+1,j))/dx**2
          phiy = (phi(i,j-1)-2.0_REAL64*phi(i,j)+phi(i,j+1))/dy**2
          drms = drms + phix + phiy
        END DO
      END DO

      ratio = etot  ! /REAL(drms/((nx+1)*(ny+1)), KIND=REAL64)

      ! Update tol value to this in case drms equals to zero
      tol = REAL(drms/((nx+1)*(ny+1)), KIND=REAL64) * 1e-5

      ! PRINT*, etot/REAL(drms/((nx+1)*(ny+1)), KIND=REAL64)

    END DO

  END SUBROUTINE gauss_seidel

  SUBROUTINE elecF(nx, ny, dx, dy, phi, Ex, Ey)

    INTEGER(INT32) :: i = 1, j = 1
    INTEGER(INT32), INTENT(IN) :: nx, ny
    REAL(REAL64), INTENT(IN) :: dx, dy
    REAL(REAL64), DIMENSION(0:nx+2,0:ny+2), INTENT(IN) :: phi
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: Ex, Ey

    ALLOCATE(Ex(1:nx+1,1:ny+1), Ey(1:nx+1,1:ny+1))

    DO i = 1,nx+1
      DO j = 1,ny+1
        Ex(i,j) = (phi(i+1,j) - phi(i-1,j))/(2.0_REAL64*dx)
        Ey(i,j) = (phi(i,j+1) - phi(i,j-1))/(2.0_REAL64*dy)
      END DO
    END DO

  END SUBROUTINE elecF


END MODULE helpers
