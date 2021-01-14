PROGRAM main

  ! We import the modules that are utilized by this program.
  ! USE create_axis
  ! USE write_netcdf
  USE helpers
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  ! REAL(REAL64), PARAMETER :: tol=1e-5
  INTEGER(INT32) :: nx, ny, i = 1, j = 1
  CHARACTER(LEN=100) :: filename
  REAL(REAL64) :: dx, dy
  ! REAL(REAL64) :: phix, phiy, etot = 0.0_REAL64, drms = 0.0_REAL64, ratio
  REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: phi, rho, Ex, Ey
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

  filename = 'poisson.txt'

  ! 5 in each directrion plus zero
  nx = 50
  ny = 50

  ! Space step calculation
  dx = 2.0_REAL64/REAL(nx, KIND=real64)
  dy = 2.0_REAL64/REAL(ny, KIND=real64)

  CALL init_rho(x, y, nx, ny, dx, dy, rho)
  CALL gauss_seidel(nx, ny, dx, dy, rho, phi)
  CALL elecF(nx, ny, dx, dy, phi, Ex, Ey)

  open(9, file="potential.txt",form='formatted')
    23 FORMAT(3(ES23.12E3))
    ! Choose format for text file ( real numbers )
    do i=1,nx+1
      do j = 1, ny+1
        write(9,23) phi(i,j), Ex(i,j), Ey(i,j)
      end do
    end do
  close(9)


  ! PRINT *, phi

  DEALLOCATE(phi, rho, Ex, Ey)

END PROGRAM main
