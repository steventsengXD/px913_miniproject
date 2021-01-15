PROGRAM main

  ! We import the modules that are utilized by this program.
  USE helpers
  USE write_netcdf
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  ! Defining variable for electron using our defined type
  TYPE(electron) :: elle
  INTEGER(INT32), PARAMETER :: steps = 1000
  INTEGER(INT32) :: nx, ny, i = 1, j = 1
  CHARACTER(LEN=30) :: filename
  CHARACTER(LEN=10) :: ic
  REAL(REAL64) :: dx=0, dy=0
  ! REAL(REAL64) :: phix, phiy, etot = 0.0_REAL64, drms = 0.0_REAL64, ratio
  REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: phi, rho, Ex, Ey


  CALL parse(nx, ny, ic)

  filename = 'electron.nc4'

  CALL initialize_sim(nx, ny, dx, dy, steps, rho, elle, ic)
  Print *,dx,dy
  CALL gauss_seidel(nx, ny, dx, dy, rho, phi)
  CALL elecF(nx, ny, dx, dy, phi, Ex, Ey)
  CALL motion(dx, dy, steps, Ex, Ey, elle)

  CALL write_file(filename, rho, phi, Ex, Ey, elle)

  DEALLOCATE(phi, rho, Ex, Ey)
   

END PROGRAM main
