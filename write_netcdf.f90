!This code was written by Charlotte Rogerson

!This code writes and outputs the NetCDF4 file
!Error messages are included in most steps 
!so any part which doesn't compile should be 
!easier to identify


MODULE write_netcdf

  USE ISO_FORTRAN_ENV
  USE netcdf
  
  IMPLICIT NONE
  
  CONTAINS
  
  SUBROUTINE write_file(filename, ierr, rho, phi, 
  
  !!!!!!!!!!!!!!!!! Inputs into the subroutine !!!!!!!!!!!!!!!!!
  
  
  
  !!!!!!!!!!!!!!!!! Creating the variable id's !!!!!!!!!!!!!!!!!
  INTEGER :: var_id_rho, var_id_phi, var_id_ex, var_id_ey, var_id_a
  INTEGER :: var_id_v, var_id_pos
  
  !Dimension id's
  INTEGER, DIMENSION(2) :: dim_id_rho, dim_id_phi, dim_id_ex, dim_id_ey
  INTEGER, DIMENSION(2) :: dim_id_a, dim_id_v, dim_id_pos
  
  !Loop variable
  INTEGER :: i
  !File id
  INTEGER :: file_id
  
  
  !!!!!!!!!!!!!!!!! Creating the file !!!!!!!!!!!!!!!!!
  ! - caution: will overwrite if it already exists
  ierr = nf90_create(filename, NF90_CLOBBER, file_id)
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF
  
  !!!!!!!!!!!!!!!!!! Variable sizes/shape !!!!!!!!!!!!!!!!!
  size_rho = SHAPE(rho)
  size_phi = SHAPE(phi)
  size_ex = SHAPE(ex)
  size_ey = SHAPE(ey)
  size_a = SHAPE(a)
  size_v = SHAPE(v)
  size_pos = SHAPE(pos)
  
  
  !Error messages if code did not run correctly
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF
  
  
  !!!!!!!!!!!!!!!!! Global attributes from the main code !!!!!!!!!!!!!!!!!
  !Problem
  ierr = nf90_put_att(file_id, NF90_GLOBAL, "problem", 
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF
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
    ierr = nf90_def_dim(file_id, rho(i),size_rho(i), dim_id_rho(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Charge Density", NF90_DOUBLE, dim_id_rho, var_id_rho)
  
  !The data for the potential, phi
  DO i=1,2
    ierr = nf90_def_dim(file_id, phi(i),size_phi(i), dim_id_phi(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Electric Potential", NF90_DOUBLE, dim_id_phi, var_id_phi)  

  !Electric field in the x-direction Ex
   DO i=1,2
    ierr = nf90_def_dim(file_id, ex(i),size_ex(i), dim_id_ex(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Ex_field", NF90_DOUBLE, dim_id_ex, var_id_ex)
  
  !Electric field in the y-direction Ey
   DO i=1,2
    ierr = nf90_def_dim(file_id, ey(i),size_ey(i), dim_id_ey(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Ey_field", NF90_DOUBLE, dim_id_ey, var_id_ey)  
  
  !The data for electron position
   DO i=1,2
    ierr = nf90_def_dim(file_id, pos(i),size_pos(i), dim_id_pos(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Position", NF90_DOUBLE, dim_id_pos, var_id_pos)  
  
  !The data for electron velocity
   DO i=1,2
    ierr = nf90_def_dim(file_id, vel(i),size_v(i), dim_id_v(i))
    IF (ierr /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END DO
  ierr = nf90_def_var(file_id, "Velocity", NF90_DOUBLE, dim_id_v, var_id_v)  
  
  !The data for electron acceleration
   DO i=1,2
    ierr = nf90_def_dim(file_id, accel(i),size_a(i), dim_id_a(i))
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
  END I
  !Electric field in the x-direction, Ex
  ierr = nf90_put_var(file_id, var_id_ex, !CHANGE THE LAST ONE HERE!!!)
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF  
  !Electric field in the y-direction, Ey
  ierr = nf90_put_var(file_id, var_id_ey,   
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF  
  !Electron position
  ierr = nf90_put_var(file_id, var_id_pos,   
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF  
  !Electron velocity
  ierr = nf90_put_var(file_id, var_id_v,   
  IF (ierr /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(ierr))
    RETURN
  END IF  
  !Electron acceleration
  ierr = nf90_put_var(file_id, var_id_a,   
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
  
  END SUBROUTINE write_file
  
END MODULE write_netcdf
