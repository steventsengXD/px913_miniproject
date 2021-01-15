! This module contains the subroutines to parse command line arguments,
! intialize the simulation, get rho, the potential, electric field and then
! move the electron

MODULE initialize

  ! We import the modules that are utilized by these subroutines/functions
  USE command_line
  USE domain_tools
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  ! Defined type to create to act as our x and y axes
  TYPE xygrid
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: x, y
  END TYPE xygrid


  ! Defined type for electron where r is position, v is velocity, a is acceleration
  TYPE electron
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: r, v, a
  END TYPE electron

  CONTAINS

  ! This subroutine collects/parses through each of the command line arguments
  ! to ensure that they meet the minimal criteria for the program to run.
  SUBROUTINE parse(nx, ny, ic)

    INTEGER(INT32) :: nx, ny
    CHARACTER(LEN=10) :: ic
    LOGICAL :: success, exists

    ! We utilize the parse_args subroutine to collect and store all the
    ! command line arguments
    CALL parse_args

    ! get_arg is used on each argument individually to access them. If input for
    ! nx, ny aren't integers then the program terminates.
    success = get_arg("nx", nx, exists=exists)
    IF (.NOT. success) THEN
      STOP "An integer was not entered for nx! Simulation terminated!"
    END IF

    success = get_arg("ny", ny, exists=exists)
    IF (.NOT. success) THEN
      STOP "An integer was not entered for ny! Simulation terminated!"
    END IF

    ! The following checks to make sure one of three problems was stated
    ! If not, the program terminates
    success = get_arg("problem", ic)
    IF(success)THEN
      IF (ic == "null" .OR. ic == "single" .OR. ic == "double") THEN
        CONTINUE
      ELSE
        STOP "A proper problem (null, single or double) was not stated! Simulation terminated!"
      END IF
    ELSE
      STOP "A proper problem (null, single or double) was not stated! Simulation terminated!"
    END IF

  END SUBROUTINE parse


  ! This module initalizes the simulation based on the command line input and
  ! calculates the values for rho
  SUBROUTINE initialize_sim(nx, ny, dx, dy, steps, rho, elle, ic)

    ! Declare axes and electron variable
    TYPE(xygrid) :: xy
    TYPE(electron), INTENT(INOUT) :: elle
    ! Set number of ghost nodes to 1 (on each side)
    INTEGER(INT32), PARAMETER :: ng = 1
    INTEGER(INT32), INTENT(IN) :: nx, ny, steps
    INTEGER(INT32) :: i, j
    ! This holds the lower bound and upper bound for the axes
    REAL(REAL64), DIMENSION(2) :: xyrange
    REAL(REAL64), INTENT(INOUT) :: dx, dy
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: rho
    CHARACTER(LEN=10), INTENT(IN) :: ic

    ! Allocating the electron variable for positon and velocity
    ALLOCATE(elle%r(0:steps,2))
    ALLOCATE(elle%v(0:steps,2))
    ! Allocating the rho variable/grid
    ALLOCATE(rho(0:nx+1,0:ny+1))

    ! Setting our x, y axes to be between -1 and 1
    xyrange = (/ -1.0_REAL64, 1.0_REAL64 /)

    ! Creating our axes using the create_axis code provided
    CALL create_axis(xy%x, nx, xyrange, ng)
    CALL create_axis(xy%y, ny, xyrange, ng)

    ! Space step calculation
    dx = xy%x(2)-xy%x(1)
    dy = xy%y(2)-xy%y(1)

    ! This if-else loop sets the initial conditions (ICs) for the problem depending on
    ! if null, single or double is input at command line
    ! For all three, initial position is first set, then velocity and then rho
    ! The null problem is initialized here
    IF (ic == 'null') THEN
      elle%r(0,:) = (/0.0_REAL64, 0.0_REAL64/)
      elle%v(0,:) = (/0.1_REAL64, 0.1_REAL64/)
      rho = 0.0_REAL64

    ! Single problem is a Gaussian peak at the origin and is initialized here
    ELSE IF (ic == 'single') THEN
      elle%r(0,:) = (/0.1_REAL64, 0.0_REAL64/)
      elle%v(0,:) = (/0.0_REAL64, 0.0_REAL64/)
      DO i = 1,nx
        DO j = 1,ny
          rho(i,j) = EXP(-(xy%x(i)/0.1_REAL64)**2 - (xy%y(j)/0.1_REAL64)**2)
        END DO
      END DO

    ! Double problem is two Gaussian peaks of different widths and is initialized here
    ELSE IF (ic == 'double') THEN
      elle%r(0,:) = (/0.0_REAL64, 0.5_REAL64/)
      elle%v(0,:) = (/0.0_REAL64, 0.0_REAL64/)
      DO i = 1,nx
        DO j = 1,ny
          rho(i,j) = EXP(-((xy%x(i)+0.250_REAL64)/0.10_REAL64)**2 - &
          ((xy%y(j)+0.250_REAL64)/0.10_REAL64)**2) + &
          EXP(-((xy%x(i)-0.750_REAL64)/0.20_REAL64)**2 - &
          ((xy%y(j)-0.750_REAL64)/0.20_REAL64)**2)
        END DO
      END DO
    END IF

    PRINT *,"This simulation is for the problem: ", ic
    PRINT '(2(A,F5.2))', "Electron's initial position: ", elle%r(0,1),", ", elle%r(0,2)

    DEALLOCATE(xy%x, xy%y)

  END SUBROUTINE initialize_sim

END MODULE initialize
