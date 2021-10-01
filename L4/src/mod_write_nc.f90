MODULE MOD_WRITE_NC

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !! Very simple for writing a 2D+T array into a Netcdf file
    !! Write  a 2D+T array to a Netcdf file
    !!
    !! Laurent Brodeau, 2011    laurent@misu.su.se
    !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    USE netcdf

    IMPLICIT none

    REAL(4), PARAMETER, PUBLIC :: rmval = -9999. ! missing value

    PRIVATE

    PUBLIC :: WRITE_NC, WRITE_NC_2D

    CHARACTER(len=200), PARAMETER :: &
         &    cv_lon = 'x',  &
         &    cv_lat = 'y',  &
         &    cv_t = 'time'

    CHARACTER(len=512), PARAMETER :: &
         &  cabout = 'Created by StVenant shallow water model. Laurent Brodeau, 2013'

CONTAINS

    SUBROUTINE WRITE_NC(vlon, vlat, vtime, XFLD, cf_out, cv_out, lmask)

        !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !!
        !! Write  a 2D+T array to a Netcdf file
        !!
        !! INPUT:
        !!        * vlon   => longitude or "x" 1D array (REAL)
        !!        * vlat   => latitude  or "y" 1D array (REAL)
        !!        * vtime  => time 1D array  (REAL)
        !!        * XFLD   => 3D (2D+T) array to write  (REAL)
        !!        * cf_out => name of the file to be writen (CHARACTER)
        !!        * cv_out => name of the variable for field XFLD (CHARACTER)
        !!
        !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        REAL, DIMENSION(:), INTENT(in) :: vlon, vlat, vtime
        REAL(4), DIMENSION(:, :, :), INTENT(in) :: XFLD
        CHARACTER(len=*), INTENT(in) :: cf_out, cv_out
        LOGICAL, OPTIONAL, INTENT(in) :: lmask

        INTEGER :: &
             &   jt, ni, nj, nt, &
             &   ierr, id_out, id_lon, id_lat, id_t, id_x, id_y, id_time, id_fld

        LOGICAL :: lm = .FALSE.

        IF (present(lmask)) THEN
            IF (lmask) lm = .TRUE.
        END IF

        PRINT *, ''
        PRINT *, 'Going to write variable ', trim(cv_out), ' into file ', trim(cf_out), ' !'

        !! Get the dimension of XFLD:
        ni = size(XFLD, 1)
        nj = size(XFLD, 2)
        nt = size(XFLD, 3)

        IF ((size(vlon, 1) /= ni) .OR. (size(vlat, 1) /= nj) .OR. (size(vtime, 1) /= nt)) THEN
            PRINT *, 'ERROR (mod_write_nc.f90): Problem of shape!!!'
            PRINT *, '   * size of vx:', size(vlon, 1)
            PRINT *, '   * size of vy:', size(vlat, 1)
            PRINT *, '   * nb of records:', size(vtime, 1)
            PRINT *, '  Shape of ', trim(cv_out), ' : ', ni, nj, nt
            STOP
        END IF

        !! Create file with ID id_out
        ierr = NF90_CREATE(cf_out, NF90_CLOBBER, id_out); call what_error(ierr)

        !! DIMENSIONS
        !! ~~~~~~~~~~

        !! Create longitude DIMENSION with ID id_x:
        ierr = NF90_DEF_DIM(id_out, trim(cv_lon), ni, id_x); call what_error(ierr)

        !! Create latitude DIMENSION with ID id_y:
        ierr = NF90_DEF_DIM(id_out, trim(cv_lat), nj, id_y); call what_error(ierr)

        !! Create time DIMENSION with ID id_t:
        ierr = NF90_DEF_DIM(id_out, trim(cv_t), nf90_unlimited, id_t); call what_error(ierr)

        !! VARIABLES
        !! ~~~~~~~~~~

        !! Create longitude VARIABLE with ID id_lon:
        ierr = NF90_DEF_VAR(id_out, trim(cv_lon), nf90_double, id_x, id_lon); call what_error(ierr)

        !! Create latitude VARIABLE with ID id_lat:
        ierr = NF90_DEF_VAR(id_out, trim(cv_lat), nf90_double, id_y, id_lat); call what_error(ierr)

        !! Create time VARIABLE with ID id_time:
        ierr = NF90_DEF_VAR(id_out, trim(cv_t), nf90_double, id_t, id_time); call what_error(ierr)

        !! Create OUTPUT FIELD VARIABLE (3D -> ni,nj,nt):
        ierr = NF90_DEF_VAR(id_out, trim(cv_out), nf90_float, (/id_x, id_y, id_t/), id_fld)
        call what_error(ierr)

        !! ATRIBUTES (bla bla, units...), for each variable and global file

        !! => NF90_PUT_ATT(...)

        !IF (lm) THEN
        !   ierr = NF90_PUT_ATT(id_out, id_fld, '_FillValue', rmval); call what_error(ierr)
        !END IF

        !! .......................
        !! .......................
        !! .......................

        ierr = NF90_PUT_ATT(id_out, NF90_GLOBAL, 'About', trim(cabout))
        call what_error(ierr)

        !! DEFINITION IS OVER
        !! ~~~~~~~~~~~~~~~~~~
        ierr = NF90_ENDDEF(id_out); call what_error(ierr)

        !! TIME TO WRITE THE VARIABLES
        !! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        !! Writing longitude:
        ierr = NF90_PUT_VAR(id_out, id_lon, vlon); call what_error(ierr)

        !! Writing latitude:
        ierr = NF90_PUT_VAR(id_out, id_lat, vlat); call what_error(ierr)

        !! Writing time:
        ierr = NF90_PUT_VAR(id_out, id_time, vtime); call what_error(ierr)

        DO jt = 1, nt

            ierr = NF90_PUT_VAR(id_out, id_fld, XFLD(:, :, jt), start=(/1, 1, jt/), &
                 &              count=(/ni, nj, 1/))
            CALL what_error(ierr)

        END DO

        ierr = NF90_CLOSE(id_out); call what_error(ierr)

        PRINT *, 'Done.'; PRINT *, ''

    END SUBROUTINE WRITE_NC

    SUBROUTINE WRITE_NC_2D(vlon, vlat, XFLD, cf_out, cv_out)

        !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !!
        !! Write a 2D array to a Netcdf file without time record...
        !!
        !! INPUT:
        !!        * vlon   => longitude or "x" 1D array (REAL)
        !!        * vlat   => latitude  or "y" 1D array (REAL)
        !!        * XFLD   => 3D (2D+T) array to write  (REAL)
        !!        * cf_out => name of the file to be writen (CHARACTER)
        !!        * cv_out => name of the variable for field XFLD (CHARACTER)
        !!
        !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        REAL, DIMENSION(:), INTENT(in) :: vlon, vlat
        REAL(4), DIMENSION(:, :), INTENT(in) :: XFLD
        CHARACTER(len=*), INTENT(in) :: cf_out, cv_out

        INTEGER :: &
             &   ni, nj, ierr, id_out, id_lon, id_lat, id_x, id_y, id_fld

        !! Get the dimension of XFLD:
        ni = size(XFLD, 1)
        nj = size(XFLD, 2)

        IF ((size(vlon, 1) /= ni) .OR. (size(vlat, 1) /= nj)) THEN
            PRINT *, 'ERROR (mod_write_nc.f90): Problem of shape!!!'
            PRINT *, '   * size of vx:', size(vlon, 1)
            PRINT *, '   * size of vy:', size(vlat, 1)
            PRINT *, '  Shape of ', trim(cv_out), ' : ', ni, nj
            STOP
        END IF

        !! Create file with ID id_out
        ierr = NF90_CREATE(cf_out, NF90_CLOBBER, id_out); call what_error(ierr)

        !! DIMENSIONS
        ierr = NF90_DEF_DIM(id_out, trim(cv_lon), ni, id_x); call what_error(ierr)
        ierr = NF90_DEF_DIM(id_out, trim(cv_lat), nj, id_y); call what_error(ierr)

        !! VARIABLES
        ierr = NF90_DEF_VAR(id_out, trim(cv_lon), nf90_double, id_x, id_lon); call what_error(ierr)
        ierr = NF90_DEF_VAR(id_out, trim(cv_lat), nf90_double, id_y, id_lat); call what_error(ierr)
        ierr = NF90_DEF_VAR(id_out, trim(cv_out), nf90_float, (/id_x, id_y/), id_fld); call what_error(ierr)

        ierr = NF90_PUT_ATT(id_out, NF90_GLOBAL, 'About', trim(cabout))
        call what_error(ierr)

        ierr = NF90_ENDDEF(id_out); call what_error(ierr)

        !! Writing...
        ierr = NF90_PUT_VAR(id_out, id_lon, vlon); call what_error(ierr)
        ierr = NF90_PUT_VAR(id_out, id_lat, vlat); call what_error(ierr)
        ierr = NF90_PUT_VAR(id_out, id_fld, XFLD); call what_error(ierr)
        ierr = NF90_CLOSE(id_out); call what_error(ierr)

    END SUBROUTINE WRITE_NC_2D

    SUBROUTINE what_error(ierror)

        IMPLICIT none

        INTEGER, INTENT(in) :: ierror

        IF (ierror /= 0) THEN
            PRINT *, 'There was an error, # = ', ierror
            STOP
        END IF

    END SUBROUTINE what_error

END MODULE MOD_WRITE_NC
