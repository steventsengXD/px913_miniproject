PROGRAM main

  ! We import the modules that are utilized by this program.
  ! USE create_axis
  ! USE write_netcdf
  USE helpers
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  ! Defining variable for electron using our defined type
  TYPE(particle) :: elle
  INTEGER(INT32) :: nx, ny, steps, i = 1, j = 1
  CHARACTER(LEN=30) :: filename
  CHARACTER(LEN=10) :: ic
  REAL(REAL64) :: dx = 0.0_REAL64, dy = 0.0_REAL64
  ! REAL(REAL64) :: phix, phiy, etot = 0.0_REAL64, drms = 0.0_REAL64, ratio
  REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: phi, rho, Ex, Ey
  ! REAL(REAL64), DIMENSION(:), ALLOCATABLE :: x, y

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
  ic = 'double'

  ! 5 in each directrion plus zero
  nx = 20_INT32
  ny = 20_INT32
  steps = 1000_INT32

  CALL init_conditions(nx, ny, dx, dy, steps, rho, elle, ic)
  Print *,dx,dy
  CALL gauss_seidel(nx, ny, dx, dy, rho, phi)
  CALL elecF(nx, ny, dx, dy, phi, Ex, Ey)
  CALL motion(dx, dy, steps, Ex, Ey, elle)

  open(9, file="potential.txt",form='formatted')
    23 FORMAT(3(ES23.12E3))
    ! Choose format for text file ( real numbers )
    do i=1,nx
      do j = 1, ny
        write(9,23) phi(i,j), Ex(i,j), Ey(i,j)
      end do
    end do
  close(9)


  ! PRINT *, phi

  DEALLOCATE(phi, rho, Ex, Ey)

END PROGRAM main
