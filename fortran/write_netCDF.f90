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

  CHARACTER(LEN=500)                                     :: ofile
  INTEGER                                                :: oid, rcode
  INTEGER                                                :: ix, iy, it
  INTEGER                                                :: dimx, dimy, dimt, dimxid, dimyid, dimtid
  INTEGER                                                :: vardid, varnumid
  CHARACTER(LEN=50)                                      :: section, varname, varn
  CHARACTER(LEN=250)                                     :: messg, text
  INTEGER, DIMENSION(3)                                  :: vardims, idstart, idcount
  REAL, ALLOCATABLE, DIMENSION(:,:,:)                    :: varvals
  INTEGER, DIMENSION(1)                                  :: dimvdims
  INTEGER, ALLOCATABLE, DIMENSION(:)                     :: dimxvals, dimyvals
  REAL, ALLOCATABLE, DIMENSION(:)                        :: dimtvals
  REAL                                                   :: fillvalue

!!!!!!! 

  fillvalue = 1.e20

  dimx = 125
  dimy = 100
  dimt = 10
  dimxid = 1
  dimyid = 2
  dimtid = 3

  varname = 'testvar'
  vardims = (/ dimxid, dimyid, dimtid /)
  varnumid = 0

  ofile = "netCDF_test.nc"

!!!!!!! !!!!!! !!!!! !!!! !!! !! !

  section="write_netCDF"
  
  rcode = nf90_create(TRIM(ofile), NF90_CLOBBER, oid)
  !CALL error_nc(section, rcode)

!!
! Deffinition section
!!  
  rcode = nf90_def_dim(oid, 'x', dimx, dimxid)
  !CALL error_nc(section, rcode)
  rcode = nf90_def_dim(oid, 'y', dimy, dimyid)
  !CALL error_nc(section, rcode)
  rcode = nf90_def_dim(oid, 'time', dimt, dimtid)
  !CALL error_nc(section, rcode)

! Coordinate variables
!!
  varn = 'x'
  vardid = 1
  dimvdims = (/ dimxid /)

  varnumid = varnumid + 1
  rcode = nf90_def_var(oid, TRIM(varn), NF90_DOUBLE, dimvdims, varnumid)
  !CALL error_nc(section, rcode)

  rcode = nf90_put_att(oid, varnumid, "units", "-")
  rcode = nf90_put_att(oid, varnumid, "standard_name", "x")
  rcode = nf90_put_att(oid, varnumid, "long_name", "x-axis")

  varn = 'y'
  vardid = 2
  dimvdims = (/ dimyid /)

  varnumid = varnumid + 1
  rcode = nf90_def_var(oid, TRIM(varn), NF90_DOUBLE, dimvdims, varnumid)
  !CALL error_nc(section, rcode)

  rcode = nf90_put_att(oid, varnumid, "units", "-")
  rcode = nf90_put_att(oid, varnumid, "standard_name", "y")
  rcode = nf90_put_att(oid, varnumid, "long_name", "y-axis")

  varn = 'time'
  vardid = 3
  dimvdims = (/ dimtid /)

  varnumid = varnumid + 1
  rcode = nf90_def_var(oid, TRIM(varn), NF90_DOUBLE, dimvdims, varnumid)
  !CALL error_nc(section, rcode)

  rcode = nf90_put_att(oid, varnumid, "units", "hours since 1949-12-01 00:00:00")
  rcode = nf90_put_att(oid, varnumid, "standard_name", "time")
  rcode = nf90_put_att(oid, varnumid, "long_name", "time")

! Variable definition
!!
  varn = varname

  varnumid = varnumid + 1
  rcode = nf90_def_var(oid, TRIM(varn), NF90_FLOAT, vardims, varnumid)
  !CALL error_nc(section, rcode)

  rcode = nf90_put_att(oid, varnumid, "units", "--")
  rcode = nf90_put_att(oid, varnumid, "coordinates", "x y")
  rcode = nf90_put_att(oid, varnumid, "standard_name", varname)
  rcode = nf90_put_att(oid, varnumid, "long_name", "netCDF test to write a variable using Fortran")
  rcode = nf90_put_att(oid, varnumid, "_FillValue", fillvalue)

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
    dimxvals(ix) = ix
  END DO

  varnumid = varnumid + 1
  rcode = nf90_put_var(oid, varnumid, dimxvals)
  !CALL error_nc(section, rcode)
  
  IF (ALLOCATED(dimyvals)) DEALLOCATE(dimyvals)
  ALLOCATE(dimyvals(dimy))

  DO iy= 1, dimy
    dimyvals(iy) = iy
  END DO

  varnumid = varnumid + 1
  rcode = nf90_put_var(oid, varnumid, dimyvals, start = (/ 1 /), count = (/ dimy /))
  !CALL error_nc(section, rcode)

  IF (ALLOCATED(dimtvals)) DEALLOCATE(dimtvals)
  ALLOCATE(dimtvals(dimt))

  DO it= 1, dimt
    dimtvals(it) = it*1.
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
