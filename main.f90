! Main driver program that calls the necessary subroutines and simulates
! electron motion


PROGRAM main

  ! We import the modules that are utilized by this program.
  USE initialize
  USE potential
  USE motion
  USE write_netcdf
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  ! Defining variable for electron using our defined type
  TYPE(electron) :: elle
  ! Number timestep iterations electron will move
  INTEGER(INT32), PARAMETER :: steps = 1000

  ! Initial condition variables which are input at command line:
  ! The number of elements in the x and y directions and the problem type
  INTEGER(INT32) :: nx, ny
  CHARACTER(LEN=10) :: ic

  ! Data file name for netCDF output
  CHARACTER(LEN=30) :: filename = 'electron.nc4'

  ! Space step variables
  REAL(REAL64) :: dx, dy

  ! Declaring the arrays for the potential, rho and electric field x, y values
  REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: phi, rho, Ex, Ey


  ! Parse our command line inputs to obtain initial conditions
  CALL parse(nx, ny, ic)

  ! Initialize the x,y axes, rho for the problem type, electron variable elle
  CALL initialize_sim(nx, ny, dx, dy, steps, rho, elle, ic)

  ! Calculate the potential using the Gauss-Seidel method
  CALL get_potential(nx, ny, dx, dy, rho, phi)

  ! Obtain electric field x and y values
  CALL elecfield(nx, ny, dx, dy, phi, Ex, Ey)

  ! Move electron using velocity verlet algorithm
  CALL advance_elec(dx, dy, steps, Ex, Ey, elle)

  ! Write the results to netCDF file for python visualization
  CALL write_file(filename, rho, phi, Ex, Ey, elle, nx, ny, ic)

  ! Deallocating our variables
  DEALLOCATE(phi, rho, Ex, Ey, elle%r, elle%v, elle%a)

END PROGRAM main
