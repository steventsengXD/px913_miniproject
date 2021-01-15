!This code was written by Charlotte Rogerson
!For PX913 Mini-Project Assessment

!Writes a NetCDF file for the following variables
!Charge density rho, potential phi, Electric field (both x and y),
!Particle position, velocity and acceleration


MODULE write_netcdf

  USE ISO_FORTRAN_ENV
  USE netcdf
  USE helpers
  
  IMPLICIT NONE
  
  CONTAINS
  !Not sure if need to include ierr?
  SUBROUTINE write_file(filename, rho, phi, Ex, Ey, ele)
  
  !!!!!!!!!!!!!!!!! Inputs into the subroutine !!!!!!!!!!!!!!!!!
  !From the main part of the code
  CHARACTER(LEN=30) :: filename !changed to match the one in main.f90
  REAL(REAL64), DIMENSION(:,:), INTENT(IN) :: phi, rho, Ex, Ey
  !Types - conatined in the helpers.f90 code
  !This gives the position (r), velocity (v) and acceleration (a)
  TYPE(electron) :: ele !electron
  INTEGER(INT32) :: nx, ny
  
  !Specific for the file
  INTEGER(INT32) :: ierr
 
    
  !!!!!!!!!!!!!!!!! Dimension names in the file !!!!!!!!!!!!!!!!!
  CHARACTER(LEN=100), DIMENSION(2) :: dim_rho = (/"rho_x", "rho_y"/)
  CHARACTER(LEN=100), DIMENSION(2) :: dim_phi = (/"phi_x", "phi_y"/)
  CHARACTER(LEN=100), DIMENSION(2) :: dim_Ex = (/"Ex_x", "Ex_y"/)
  CHARACTER(LEN=100), DIMENSION(2) :: dim_Ey = (/"Ey_x", "Ey_y"/)
  !g = spatial grid,t = time
  CHARACTER(LEN=100), DIMENSION(2) :: dim_v = (/"v_g", "v_t"/)
  CHARACTER(LEN=100), DIMENSION(2) :: dim_a = (/"a_g", "a_t"/)
  CHARACTER(LEN=100), DIMENSION(2) :: dim_r = (/"r_g", "r_t"/) 
  
  !!!!!!!!!!!!!!!!! Creating the variable id's !!!!!!!!!!!!!!!!!
  INTEGER(INT32) :: var_id_rho, var_id_phi, var_id_ex, var_id_ey
  INTEGER(INT32) :: var_id_v, var_id_r, var_id_a
  
  !Dimension id's
  INTEGER(INT32), DIMENSION(2) :: dim_id_rho, dim_id_phi, dim_id_ex, dim_id_ey
  INTEGER(INT32), DIMENSION(2) :: dim_id_a, dim_id_v, dim_id_r
  INTEGER(INT32) :: file_id
  
  !Loop variable
  INTEGER(INT32) :: i
  
  !Implicit types for the Sizes
  INTEGER(INT32), DIMENSION(2) :: sizes_rho, sizes_phi, sizes_ex, sizes_ey
  INTEGER(INT32), DIMENSION(2) :: sizes_r, sizes_v, sizes_a  
   
  !!!!!!!!!!!!!!!!! Creating the file !!!!!!!!!!!!!!!!!
  ! - caution: this will overwrite the file if it already exists in the directory
  ierr = nf90_create(filename, NF90_CLOBBER, file_id)
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF
  
  !!!!!!!!!!!!!!!!!! Variable sizes/shape !!!!!!!!!!!!!!!!!
  sizes_rho = SHAPE(rho)
  sizes_phi = SHAPE(phi)
  sizes_ex = SHAPE(Ex)
  sizes_ey = SHAPE(Ey)
  sizes_a = SHAPE(ele%a)
  sizes_v = SHAPE(ele%v)
  sizes_r = SHAPE(ele%r)
  
  !Error messages if code did not run correctly
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF
  
  !!!!!!!!!!!!!!!!! Global attributes from the main code !!!!!!!!!!!!!!!!!
  !Problem
  !ierr = nf90_put_att(file_id, NF90_GLOBAL, "problem", 
  !IF (ierr /= nf90_noerr) THEN
   ! PRINT *, TRIM(nf90_strerror(ierr))
    !RETURN
  !END IF
  !nx
  ierr = nf90_put_att(file_id, NF90_GLOBAL, "nx", nx)
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF
  !ny
  ierr = nf90_put_att(file_id, NF90_GLOBAL, "ny", ny)
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF
    
  !!!!!!!!!!!!!!!!! Inputting the data !!!!!!!!!!!!!!!!!
  !Need both columns in the data so do DO loops between i=1,2
  !The data for the density, rho
  DO i=1,2
    ierr = nf90_def_dim(file_id, dim_rho(i), sizes_rho(i), dim_id_rho(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Charge Density", NF90_DOUBLE, dim_id_rho, var_id_rho)
  
  !The data for the potential, phi
  DO i=1,2
    ierr = nf90_def_dim(file_id, dim_phi(i),sizes_phi(i), dim_id_phi(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Electric Potential", NF90_DOUBLE, dim_id_phi, var_id_phi)  

  !Electric field in the x-direction Ex
   DO i=1,2
    ierr = nf90_def_dim(file_id, dim_Ex(i),sizes_ex(i), dim_id_ex(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Ex_field", NF90_DOUBLE, dim_id_ex, var_id_ex)
  
  !Electric field in the y-direction Ey
   DO i=1,2
    ierr = nf90_def_dim(file_id, dim_Ey(i),sizes_ey(i), dim_id_ey(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Ey_field", NF90_DOUBLE, dim_id_ey, var_id_ey)  
  
  !The data for electron position
   DO i=1,2
    ierr = nf90_def_dim(file_id, dim_r(i),sizes_r(i), dim_id_r(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Position", NF90_DOUBLE, dim_id_r, var_id_r)  
  
  !The data for electron velocity
   DO i=1,2
    ierr = nf90_def_dim(file_id, dim_v(i),sizes_v(i), dim_id_v(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Velocity", NF90_DOUBLE, dim_id_v, var_id_v)  
  
  !The data for electron acceleration
   DO i=1,2
    ierr = nf90_def_dim(file_id, dim_a(i),sizes_a(i), dim_id_a(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Acceleration", NF90_DOUBLE, dim_id_a, var_id_a)  
  
  
  !!!!!!!!!!!!!!!!! Metadata !!!!!!!!!!!!!!!!!
  
  ierr = nf90_enddef(file_id)
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF  
  
  !!!!!!!!!!!!!!!!! Writing all the data to file !!!!!!!!!!!!!!!!!
  !The density, rho
  ierr = nf90_put_var(file_id, var_id_rho, rho)
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF
  !The potential, phi
  ierr = nf90_put_var(file_id, var_id_phi, phi)
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF
  !Electric field in the x-direction, Ex
  ierr = nf90_put_var(file_id, var_id_ex, Ex)
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF  
  !Electric field in the y-direction, Ey
  ierr = nf90_put_var(file_id, var_id_ey, Ey)
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF  
  !Electron position
  ierr = nf90_put_var(file_id, var_id_r, ele%r)
  IF (ierr /= nf90_noerr) THEN 
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF  
  !Electron velocity
  ierr = nf90_put_var(file_id, var_id_v,ele%v) 
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF  
  !Electron acceleration
  ierr = nf90_put_var(file_id, var_id_a,ele%a)   
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF
 
  !!!!!!!!!!!!!!!!! Closing the NetCDF file !!!!!!!!!!!!!!!!!
  ierr = nf90_close(file_id)
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF  
    
  !Statement to say everything has run 
  PRINT *, "Success in writing file: ", filename
  
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    PRINT *, "Failure in writing file"
    RETURN
  END IF   
  
  END SUBROUTINE write_file
  
END MODULE write_netcdf
