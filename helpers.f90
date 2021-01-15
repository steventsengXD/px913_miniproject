MODULE helpers

  ! We import the modules that are utilized by these subroutines/functions
  USE command_line
  USE domain_tools
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE
  SAVE

  TYPE xygrid
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: x, y
  END TYPE xygrid


  ! Defined type for particle where r is position, v is velosity, a is acceleration
  TYPE electron
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: r, v, a
  END TYPE electron

  CONTAINS

  ! This subroutine collects/parses through each of the command line arguments
  ! to ensure that they meet the minimal criteria for the program to run.
  ! Program will run a default of the null problem with nx = ny = 50 if user
  ! input doesnt make sense
  SUBROUTINE parse(nx, ny, ic)

    INTEGER(INT32) :: nx, ny
    CHARACTER(LEN=10) :: ic
    LOGICAL :: success, exists

    ! We utilize the parse_args subroutine to collect and store all the
    ! command line arguments
    CALL parse_args

    ! get_arg is used on each argument individually to access them. If input is
    ! unreasonable, a default value is set for each of the parameters.
    success = get_arg("nx", nx, exists=exists)
    IF (.NOT. success) THEN
      nx = 50
      PRINT*,"An integer was not set for nx! nx = 50 default value will be used."
      PRINT*,""
    END IF

    success = get_arg("ny", ny, exists=exists)
    IF (.NOT. success) THEN
      ny = 50
      PRINT*,"An integer was not set for ny! ny = 50 default value will be used."
      PRINT*,""
    END IF

    ! The following checks to make sure one of three problems was stated
    ! If not, the null problem is used by default
    success = get_arg("problem", ic)
    IF(success)THEN
      IF (ic == "null" .OR. ic == "Null" .OR. ic == "single" .OR. ic == "Single" &
      .OR. ic == "double" .OR. ic == "Double") THEN
        CONTINUE
      END IF
    ELSE
      ic = "null"
      PRINT*,"A problem was not stated! Problem null default will be used."
      PRINT*,""
    END IF

  END SUBROUTINE parse


  SUBROUTINE initialize_sim(nx, ny, dx, dy, steps, rho, elle, ic)

    TYPE(xygrid) :: xy
    TYPE(electron), INTENT(INOUT) :: elle
    INTEGER(INT32), PARAMETER :: ghosts = 1
    INTEGER(INT32) :: i, j
    INTEGER(INT32), INTENT(IN) :: nx, ny, steps
    REAL(REAL64), DIMENSION(2) :: xyrange
    REAL(REAL64), INTENT(INOUT) :: dx, dy
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: rho
    CHARACTER(LEN=10), INTENT(IN) :: ic

    ALLOCATE(elle%r(0:steps,2))
    ALLOCATE(elle%v(0:steps,2))
    ALLOCATE(rho(0:nx+1,0:ny+1))

    xyrange = (/ -1.0_REAL64, 1.0_REAL64 /)

    CALL create_axis(xy%x, nx, xyrange, ghosts)
    CALL create_axis(xy%y, ny, xyrange, ghosts)
    print*, xy%x, xy%y

    ! Space step calculation
    dx = 2.0_REAL64/REAL(nx-1, KIND=real64)
    dy = 2.0_REAL64/REAL(ny-1, KIND=real64)

    IF (ic == 'null') THEN

      ! Setting initial position and velocity
      elle%r(0,:) = (/0.0_REAL64, 0.0_REAL64/)
      elle%v(0,:) = (/0.1_REAL64, 0.1_REAL64/)

      rho(i,j) = 0.0_REAL64

    ELSE IF (ic == 'single') THEN

      ! Setting initial position and velocity
      elle%r(0,:) = (/0.1_REAL64, 0.0_REAL64/)
      elle%v(0,:) = (/0.0_REAL64, 0.0_REAL64/)

      DO i = 1,nx
        DO j = 1,ny
          rho(i,j) = EXP(-(xy%x(i)/0.1_REAL64)**2 - (xy%y(j)/0.1_REAL64)**2)
        END DO
      END DO

    ELSE IF (ic == 'double') THEN

      ! Setting initial position and velocity
      elle%r(0,:) = (/0.1_REAL64, 0.5_REAL64/)
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

  DEALLOCATE(xy%x, xy%y)

END SUBROUTINE initialize_sim


  SUBROUTINE gauss_seidel(nx, ny, dx, dy, rho, phi)

    REAL(REAL64) :: tol=1e-5
    INTEGER(INT32), INTENT(IN) :: nx, ny
    INTEGER(INT32) :: i, j
    REAL(REAL64), INTENT(IN) :: dx, dy
    REAL(REAL64), DIMENSION(0:nx+1, 0:ny+1), INTENT(IN) :: rho
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: phi
    REAL(REAL64) :: phix, phiy, etot = 0.0_REAL64, drms = 0.0_REAL64
    REAL(REAL64) :: denom, ratio = 1.0_REAL64

    ALLOCATE(phi(0:nx+1, 0:ny+1))

    ! Calculate denominator here first for neatness purposes
    denom = -2.0_REAL64*(1.0_REAL64/dx**2+1.0_REAL64/dy**2)

    DO WHILE(ratio>tol)

      DO i=1,nx
        DO j=1,ny
          phi(i,j) = (-(phi(i+1,j)+phi(i-1,j))/dx**2 - (phi(i,j+1)+phi(i,j-1))/dy**2 &
          + rho(i,j))/denom
        END DO
      END DO

      DO i=1,nx
        DO j=1,ny
          phi(i,j) = (-(phi(i+1,j)+phi(i-1,j))/dx**2 - (phi(i,j+1)+phi(i,j-1))/dy**2 &
          + rho(i,j))/denom
        END DO
      END DO

      etot = 0.0_REAL64
      DO i=1,nx
        DO j=1,ny
          phix = (phi(i-1,j)-2.0_REAL64*phi(i,j)+phi(i+1,j))/dx**2
          phiy = (phi(i,j-1)-2.0_REAL64*phi(i,j)+phi(i,j+1))/dy**2
          etot = etot + ABS(phix + phiy - rho(i,j))
        END DO
      END DO

      drms = 0.0_REAL64

      DO i=1,nx
        DO j=1,ny
          phix = (phi(i-1,j)-2.0_REAL64*phi(i,j)+phi(i+1,j))/dx**2
          phiy = (phi(i,j-1)-2.0_REAL64*phi(i,j)+phi(i,j+1))/dy**2
          drms = drms + phix + phiy
        END DO
      END DO

      ratio = etot  ! /REAL(drms/((nx+1)*(ny+1)), KIND=REAL64)

      ! Update tol value to this in case drms equals to zero
      tol = REAL(drms/(nx*ny), KIND=REAL64) * 1e-5

      ! PRINT*, etot/REAL(drms/((nx+1)*(ny+1)), KIND=REAL64)

    END DO

  END SUBROUTINE gauss_seidel


  ! Subroutine that sets the electric field
  ! Can probably create a 3d array or defined type for this so that
  ! only one variable is needed
  SUBROUTINE elecF(nx, ny, dx, dy, phi, Ex, Ey)

    INTEGER(INT32) :: i = 1, j = 1
    INTEGER(INT32), INTENT(IN) :: nx, ny
    REAL(REAL64), INTENT(IN) :: dx, dy
    REAL(REAL64), DIMENSION(0:nx+1,0:ny+1), INTENT(IN) :: phi
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: Ex, Ey

    ALLOCATE(Ex(1:nx,1:ny), Ey(1:nx,1:ny))

    DO i = 1,nx
      DO j = 1,ny
        Ex(i,j) = (phi(i+1,j) - phi(i-1,j))/(2.0_REAL64*dx)
        Ey(i,j) = (phi(i,j+1) - phi(i,j-1))/(2.0_REAL64*dy)
      END DO
    END DO

  END SUBROUTINE elecF


  SUBROUTINE motion(dx, dy, steps, Ex, Ey, elle)

    TYPE(electron) :: elle
    INTEGER(INT32) :: i = 1
    INTEGER(INT32) :: cell_x, cell_y
    INTEGER(INT32), INTENT(IN) :: steps
    REAL(REAL64), INTENT(IN) :: dx, dy
    REAL(REAL64), PARAMETER :: dt = 0.01_REAL64
    ! REAL(REAL64), DIMENSION(0:nx+1,0:ny+1), INTENT(IN) :: phi
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: Ex, Ey


    ! Allocate acceleration array of particle
    ALLOCATE(elle%a(0:steps,2))

    ! Calculates which cell particle is starting in
    cell_x = FLOOR((elle%r(0,1) - 1.0_REAL64)/dx) + 1
    cell_y = FLOOR((elle%r(0,2) - 1.0_REAL64)/dy) + 1

    ! Calculate initial acceleration
    elle%a(0,1) = Ex(cell_x, cell_y)
    elle%a(0,2) = Ey(cell_x, cell_y)

    ! Particle motion using velocity verlet algorithm
    DO i=1, steps

      elle%r(i,1) = elle%r(i-1,1) + elle%v(i-1,1)*dt + (elle%a(i-1,1)*dt)**2/2.0_REAL64
      elle%r(i,2) = elle%r(i-1,2) + elle%v(i-1,2)*dt + (elle%a(i-1,2)*dt)**2/2.0_REAL64

      ! Locate cell indices for electron
      cell_x = FLOOR((elle%r(i,1) - 1.0_REAL64)/dx) + 1
      cell_y = FLOOR((elle%r(i,2) - 1.0_REAL64)/dy) + 1

      ! Calculate next acceleration
      elle%a(i,1) = -Ex(cell_x, cell_y)
      elle%a(i,2) = -Ey(cell_x, cell_y)

      ! Calculate next velocity
      elle%v(i,1) = elle%v(i-1,1) + dt*(elle%a(i,1) + elle%a(i-1,1))/2.0_REAL64
      elle%v(i,2) = elle%v(i-1,2) + dt*(elle%a(i,2) + elle%a(i-1,2))/2.0_REAL64

      ! PRINT*, elle%r(i,1), elle%r(i,2)

    END DO



  END SUBROUTINE motion

END MODULE helpers
