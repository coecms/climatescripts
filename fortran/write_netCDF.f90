PROGRAM write_netCDF
! Program to write a netCDF file using Fortran 90
!! L. Fita, CCRC - ARC CoE CSS, UNSW, Sydney, Australia
!
!!! Comilation
! $NETCDFhome = variable with the location of netCDF libraries
!! gfortran 
! gfortran write_netCDF.f90 -L${NETCDFhome}/lib -lnetcdf -lnetcdff -lm -I${NETCDFhome}/include -o write_netCDF
!
!!!!!!! !!!!!! !!!!! !!!! !!! !! !
  USE netcdf

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  INTEGER                                                :: ix, iy, it
  CHARACTER(LEN=500)                                     :: ofile
  INTEGER                                                :: oid, rcode, funit, ios
  INTEGER                                                :: dimx, dimy, dimt, dimxid, dimyid, dimtid
  REAL                                                   :: first_lon, first_lat, first_time, lon_inc, lat_inc, time_inc
  INTEGER                                                :: vardid, varnumid
  CHARACTER(LEN=50)                                      :: section, varname, varn
  CHARACTER(LEN=200)                                     :: varstdname, varlonname, varunits, gridmap_name, timeunits, timecalendar
  CHARACTER(LEN=250)                                     :: projname, messg, text
  INTEGER, DIMENSION(3)                                  :: vardims, idstart, idcount
  REAL, ALLOCATABLE, DIMENSION(:,:,:)                    :: varvals
  INTEGER, DIMENSION(1)                                  :: dimvdims
  INTEGER, ALLOCATABLE, DIMENSION(:)                     :: dimxvals, dimyvals
  REAL, ALLOCATABLE, DIMENSION(:)                        :: dimtvals
  REAL                                                   :: fillvalue
  LOGICAL                                                :: is_used
  CHARACTER(LEN=200),DIMENSION(3)                        :: basic_attributes

  NAMELIST /datainf/ dimx, dimy, dimt, varname, varstdname, varlonname, varunits, ofile, fillvalue
  NAMELIST /projectioninf/ projname, first_lon, first_lat, first_time, lon_inc, lat_inc, time_inc, timeunits, timecalendar

!!!!!!! Functions & Subroutines
! def_variable: Subroutine to define a variable
! gridmap_name: Function to provide the standard CF-compilant projection name

!!!!!!! Variables
! dim[x/y/t]: dimension of the variable (assumed a 3D time, lat, lon variable)
! varname: name of the variable
! varstdname: standard name of the variable
! varlonname: long name of the variable
! varunits: units of the variable
! ofile: output filename
! fillvalue: value for missing values
! projname: name of the projection (only longitude-latitude)
! first_[lon/lat]: first longitude/latitude/time values (remember coordinates must increase!)
! [lon/lat/time]_inc: longitude, latitude and time increments
! timeunits: units of the time variable
! timecalendar: calendar of the time

! Reading namelist
!!
  DO funit=10,100
    INQUIRE(unit=funit, opened=is_used)
    IF (.not. is_used) EXIT
  END DO

  OPEN(funit,file='namelist.write_netcdf',status='old',form='formatted',iostat=ios)
  IF ( ios /= 0 ) STOP "ERROR opening namelist.write_netCDF"
  READ(funit,datainf)
  READ(funit,projectioninf)
  CLOSE(funit)

  dimxid = 1
  dimyid = 2
  dimtid = 3

  vardims = (/ dimxid, dimyid, dimtid /)
  varnumid = 0

! Running program
!!!!!!! !!!!!! !!!!! !!!! !!! !! !

  section="write_netCDF"
  
  rcode = nf90_create(TRIM(ofile), NF90_CLOBBER, oid)
  !CALL error_nc(section, rcode)

!!
! Deffinition section
!!  
  rcode = nf90_def_dim(oid, 'lon', dimx, dimxid)
  !CALL error_nc(section, rcode)
  rcode = nf90_def_dim(oid, 'lat', dimy, dimyid)
  !CALL error_nc(section, rcode)
  rcode = nf90_def_dim(oid, 'time', NF90_UNLIMITED, dimtid)
  !CALL error_nc(section, rcode)

! Coordinate variables
!!
  varn = 'lon'
  dimvdims = (/ dimxid /)
  basic_attributes(1) = "lon"
  basic_attributes(2) = "longitude"
  basic_attributes(3) = "degrees_east"

  varnumid = varnumid + 1
  CALL def_variable(oid,varn, varnumid, 1, dimvdims, NF90_DOUBLE, basic_attributes)

  varn = 'lat'
  dimvdims = (/ dimyid /)
  basic_attributes(1) = "lat"
  basic_attributes(2) = "latitude"
  basic_attributes(3) = "degrees_north"

  varnumid = varnumid + 1
  CALL def_variable(oid,varn, varnumid, 1, dimvdims, NF90_DOUBLE, basic_attributes)

  varn = 'time'
  dimvdims = (/ dimtid /)
  basic_attributes(1) = "time"
  basic_attributes(2) = "time"
  basic_attributes(3) = timeunits

  varnumid = varnumid + 1
  CALL def_variable(oid,varn, varnumid, 1, dimvdims, NF90_DOUBLE, basic_attributes)

  rcode = nf90_put_att(oid, varnumid, "calendar", TRIM(timecalendar))

! Variable definition
!!
  varn = varname
  basic_attributes(1) = varstdname
  basic_attributes(2) = varlonname
  basic_attributes(3) = varunits

  varnumid = varnumid + 1
  CALL def_variable(oid,varn, varnumid, 3, vardims, NF90_DOUBLE, basic_attributes)
  rcode = nf90_put_att(oid, varnumid, "grid_mapping", TRIM(projname))
  rcode = nf90_put_att(oid, varnumid, "coordinates", "lon lat")

! Variable projection definition
!!
  varn = projname

  varnumid = varnumid + 1
  rcode = nf90_def_var(ncid=oid, name=TRIM(varn), xtype=NF90_CHAR, varid=varnumid)
  rcode = nf90_put_att(oid, varnumid, "grid_mapping_name", TRIM(gridmap_name(projname)))

! Global attributes
!!
  rcode = nf90_put_att(oid, NF90_GLOBAL, "author", "L. Fita")
  text= 'Climate Change Research Center - ARC Center of Excellence for Climate System Science'
  rcode = nf90_put_att(oid, NF90_GLOBAL, "institution", text)

!!
! Clossing definition
!!
  rcode = nf90_enddef(oid)

!!
! Filling file
!!
  varnumid = 0
  IF (ALLOCATED(dimxvals)) DEALLOCATE(dimxvals)
  ALLOCATE(dimxvals(dimx))

  DO ix= 1, dimx
    dimxvals(ix) = first_lon + lon_inc*(ix-1)
  END DO

  varnumid = varnumid + 1
  rcode = nf90_put_var(oid, varnumid, dimxvals)
  !CALL error_nc(section, rcode)
  
  IF (ALLOCATED(dimyvals)) DEALLOCATE(dimyvals)
  ALLOCATE(dimyvals(dimy))

  DO iy= 1, dimy
    dimyvals(iy) = first_lat + lat_inc*(iy-1)
  END DO

  varnumid = varnumid + 1
  rcode = nf90_put_var(oid, varnumid, dimyvals, start = (/ 1 /), count = (/ dimy /))
  !CALL error_nc(section, rcode)

  IF (ALLOCATED(dimtvals)) DEALLOCATE(dimtvals)
  ALLOCATE(dimtvals(dimt))

  DO it= 1, dimt
    dimtvals(it) = first_time + time_inc*(it-1)
  END DO

  varnumid = varnumid + 1
  rcode = nf90_put_var(oid, varnumid, dimtvals, start = (/ 1 /), count = (/ dimt /))
  !CALL error_nc(section, rcode)

  IF (ALLOCATED(varvals)) DEALLOCATE(varvals)
  ALLOCATE(varvals(dimx, dimy, dimt))

  DO it=1,dimt
    DO iy=1,dimy
      DO ix=1,dimx
        varvals(ix,iy,it)=(it-1)*100.+SQRT((dimx/2.-ix)**2.+(dimy/2.-iy)**2.)
      END DO
    END DO
  END DO

  varnumid = varnumid + 1
  rcode = nf90_inquire_variable(oid, varnumid, dimids=idstart)

  idstart=1
  idcount=(/ dimx, dimy, dimt /)

  rcode = nf90_put_var(oid, varnumid, varvals, start=idstart, count=idcount)
  !CALL error_nc(section, rcode)

! Closing
!!
  rcode = nf90_close(oid)
  !CALL error_nc(section, rcode)

END PROGRAM write_netCDF

SUBROUTINE def_variable(ncf,varname, varnid, varNdims, vardims, vartype, basicattrs)
! Subroutine to define a variable

  USE netcdf

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  CHARACTER(LEN=50), INTENT(IN)                          :: varname
  INTEGER, INTENT(IN)                                    :: ncf,varnid, varNdims, vartype
  INTEGER, DIMENSION(varNdims), INTENT(IN)               :: vardims
  CHARACTER(LEN=200), DIMENSION(3), INTENT(IN)           :: basicattrs
  INTEGER                                                :: rc, varid

!!!!!!! Variables
! ncf: netCDF unit number
! varname: name of the variable
! varNdims: number of dimensions of the variable
! vardims: vector with the dimensions id
! vartpye: type of variable
! basicattrs: basic attributes (standard_name, long_name, units)
  varid = varnid

  rc = nf90_def_var(ncf, TRIM(varname), vartype, vardims, varid)

  rc = nf90_put_att(ncf, varnid, "standard_name", TRIM(basicattrs(1)))
  rc = nf90_put_att(ncf, varnid, "long_name", TRIM(basicattrs(2)))
  rc = nf90_put_att(ncf, varnid, "units", TRIM(basicattrs(3)))

END SUBROUTINE def_variable

CHARACTER(LEN=200) FUNCTION gridmap_name(proj)
! Function to provide the standard CF-compilant projection name
!   Following CF-conventions: http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.html#appendix-grid-mappings
  CHARACTER(LEN=200), INTENT(IN)                         :: proj

!!!!!!! Variables
! proj: projection name

  SELECT CASE (proj)
    CASE('albers_conical_equal_area')
      gridmap_name = 'Albers_Equal_Area'
    CASE('azimuthal_equidistant')
      gridmap_name = 'Azimuthal_Equidistant'
    CASE('lambert_azimuthal_equal_area')
      gridmap_name = 'Lambert_azimuthal_equal_area'
    CASE('lambert_conformal_conic')
      gridmap_name = 'Lambert_Conformal'
    CASE('lambert_cylindrical_equal_area')
      gridmap_name = 'Lambert_Cylindrical_Equal_Area'
    CASE('latitude_longitude')
      gridmap_name = 'Latitude_Longitude'
    CASE('mercator')
      gridmap_name = 'Mercator'
    CASE('orthographic')
      gridmap_name = 'Orthographic'
    CASE('polar_stereographic')
      gridmap_name = 'Polar_stereographic'
    CASE('rotated_latitude_longitude')
      gridmap_name = 'Rotated_pole'
    CASE('stereographic')
      gridmap_name = 'Stereographic'
    CASE('transverse_mercator')
      gridmap_name = 'Transverse_Mercator'
    CASE('vertical_perspective')
      gridmap_name = 'Vertical_perspective'
    CASE DEFAULT
      PRINT *, 'ERROR -- error -- ERROR -- error'
      PRINT *, '  gridmap_name: projection name "' // TRIM(proj) // '" is not ready !!!'
      STOP

  END SELECT

END FUNCTION gridmap_name
