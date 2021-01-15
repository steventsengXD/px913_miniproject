! This module contains the subroutines for simulating the electron's movement

MODULE motion

  USE potential
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  CONTAINS

  ! Subroutine that obtains the electric field
  ! Can probably create a 3d array or defined type for this so that
  ! only one variable is needed
  SUBROUTINE elecfield(nx, ny, dx, dy, phi, Ex, Ey)

    INTEGER(INT32) :: i = 1, j = 1
    INTEGER(INT32), INTENT(IN) :: nx, ny
    REAL(REAL64), INTENT(IN) :: dx, dy
    REAL(REAL64), DIMENSION(0:nx+1,0:ny+1), INTENT(IN) :: phi

    ! The x,y values of the electric field are stored in these two 2D arrays
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Ex, Ey

    ! The arrays are allocated here without ghost nodes as they arent necessary
    ALLOCATE(Ex(1:nx,1:ny))
    ALLOCATE(Ey(1:nx,1:ny))

    ! Iterate over the phi array to obtain the x and y electric field values
    DO i = 1,nx
      DO j = 1,ny
        Ex(i,j) = (phi(i+1,j) - phi(i-1,j))/(2.0_REAL64*dx)
        Ey(i,j) = (phi(i,j+1) - phi(i,j-1))/(2.0_REAL64*dy)
      END DO
    END DO

  END SUBROUTINE elecfield


  ! This subroutine moves the electron using the velocity verlet algorithm
  SUBROUTINE advance_elec(dx, dy, steps, Ex, Ey, elle)

    TYPE(electron), INTENT(INOUT) :: elle
    INTEGER(INT32) :: i = 1
    INTEGER(INT32) :: cell_x, cell_y
    INTEGER(INT32), INTENT(IN) :: steps
    REAL(REAL64), INTENT(IN) :: dx, dy
    ! Set the timestep size
    REAL(REAL64), PARAMETER :: dt = 0.01_REAL64
    REAL(REAL64), DIMENSION(1:,1:), INTENT(IN) :: Ex, Ey

    ! Allocate acceleration array of electron particle
    ALLOCATE(elle%a(0:steps,2))

    ! Calculates which cell particle is starting in
    ! The code was changed from the instructions, the minus sign was replaced by
    ! a plus sign as a minus gave negative grid cell indices which isn't correct
    cell_x = FLOOR((elle%r(0,1) + 1.0_REAL64)/dx) + 1
    cell_y = FLOOR((elle%r(0,2) + 1.0_REAL64)/dy) + 1

    ! PRINT *, cell_x,cell_y

    ! Calculate initial acceleration
    elle%a(0,1) = -Ex(cell_x, cell_y)
    elle%a(0,2) = -Ey(cell_x, cell_y)

    ! Iterate over electron's previous position, velocity and acceleration to
    ! obtain motion using velocity verlet algorithm
    DO i=1, steps

      ! Obtaining the next position
      elle%r(i,1) = elle%r(i-1,1) + elle%v(i-1,1)*dt + (elle%a(i-1,1)*dt)**2/2.0_REAL64
      elle%r(i,2) = elle%r(i-1,2) + elle%v(i-1,2)*dt + (elle%a(i-1,2)*dt)**2/2.0_REAL64

      ! Locate cell indices for electron
      cell_x = FLOOR((elle%r(i,1) + 1.0_REAL64)/dx) + 1
      cell_y = FLOOR((elle%r(i,2) + 1.0_REAL64)/dy) + 1

      ! Calculate next acceleration, -1 due to negative charge of electron
      elle%a(i,1) = -Ex(cell_x, cell_y)
      elle%a(i,2) = -Ey(cell_x, cell_y)

      ! Calculate next velocity
      elle%v(i,1) = elle%v(i-1,1) + dt*(elle%a(i,1) + elle%a(i-1,1))/2.0_REAL64
      elle%v(i,2) = elle%v(i-1,2) + dt*(elle%a(i,2) + elle%a(i-1,2))/2.0_REAL64

    END DO

    PRINT '(2(A,F5.2))', "Electron's final position: ", elle%r(steps,1), &
    ", ", elle%r(steps,2)

  END SUBROUTINE advance_elec

END MODULE motion
