MODULE helpers

  ! We import the modules that are utilized by these subroutines/functions
  ! USE command_line
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE init_rho(x, y, nx, ny, rho)

    INTEGER(INT32) :: i=1, j=1, nx, ny
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: x, y
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: rho

    ! ALLOCATE(rho(2*nx+1,2*ny+1))

    DO i = 1,nx+1
      DO j = 1,ny+1
        rho(i,j) = EXP(-(x(i)/0.1_REAL64)**2 - (y(j)/0.1_REAL64)**2)
      END DO
    END DO

  END SUBROUTINE init_rho



END MODULE helpers



PROGRAM main

  ! We import the modules that are utilized by this program.
  ! USE create_axis
  ! USE write_netcdf
  USE helpers
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  REAL(REAL64), PARAMETER :: tol=1e-5
  INTEGER(INT32) :: nx, ny, i = 1, j = 1
  CHARACTER(LEN=100) :: filename
  REAL(REAL64) :: dx, dy, denom
  REAL(REAL64) :: phix, phiy, etot = 0.0_REAL64, drms = 0.0_REAL64, ratio
  REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: phi, rho
  REAL(REAL64), DIMENSION(:), ALLOCATABLE :: x, y

  ! ! We use the subroutine 'parse' to check our command line arguments to make
  ! ! sure they are suitable for this simulation
  ! CALL parse(N, T, init, beta, J)
  !

  ! ! The axes for x, y, t are created here
  ! CALL define_axis(xaxis, N, 1)
  ! CALL define_axis(yaxis, N, 1)
  ! CALL define_axis(taxis, T, 0)
  !
  ! ! The data/axes are written to file here
  ! CALL write_data(grid, grid_init, mag, xaxis, yaxis, taxis, init, beta, J, &
  ! filename, ierr)
  !


  filename = 'poisson.txt'

  ! 5 in each directrion plus zero
  nx = 10
  ny = 10

  ! Space step calculation
  dx = 2.0_REAL64/REAL(nx, KIND=real64)
  dy = 2.0_REAL64/REAL(ny, KIND=real64)

  ALLOCATE(phi(0:nx+2,0:ny+2), rho(0:nx+2,0:ny+2))
  ALLOCATE(x(0:nx+2),y(0:ny+2))

  x(1)=-1
  DO i=1,nx+1
    x(i+1) = x(i) + dx
  END DO

  y(1)=-1
  DO i=1,ny+1
    y(i+1) = y(i) + dy
  END DO

  CALL init_rho(x, y, nx, ny, rho)

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

    ratio = etot/REAL(drms/((nx+1)*(ny+1)), KIND=REAL64)

  END DO

  PRINT *, phi


  DEALLOCATE(phi, rho, x, y)

END PROGRAM main
