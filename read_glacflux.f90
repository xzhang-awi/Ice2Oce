program test
  ! PG: We need to load the netcdf module:
  use netcdf
  implicit none
  INTEGER*8 :: IUGLAC1, IUGLAC2, IUGLAC3, IUGLAC4
  REAL*8, allocatable :: GLACFLUXG(:,:), GLACFLUX2(:,:), GLACFLUX3(:,:)
  CHARACTER (len=16) :: glacfname
  CHARACTER (len=16) :: glacfname2
  CHARACTER (len=16) :: glacfnamenc
  CHARACTER (len=16) :: fw_hosing_varname
  INTEGER :: ncid, varid, x, y  ! PG: Allocate stuff for the netcdf test
  INTEGER, parameter :: ifid_glacflux=552
  INTEGER, parameter :: ifid_glacflux2=553
  INTEGER, parameter :: ifid_glacflux3=554
  ALLOCATE(GLACFLUXG(254, 220))
  write (*,*) "*** EXT TEST ***"
  glacfname="FBFLOMPIOM.ext8"
  open(ifid_glacflux, file=trim(glacfname), form="unformatted", status="old")
  write (*,*) "File opened!"
  read(ifid_glacflux) IUGLAC1, IUGLAC2, IUGLAC3, IUGLAC4
  write (*,*) "Header finished..."
  write (*,*) "IUGLAC1 is ", IUGLAC1, "IUGLAC2 is", IUGLAC2, "IUGLAC3 is ", IUGLAC3, "IUGLAC4 is", IUGLAC4
  write (*,*) "Reading GLACFLUXG..."
  read(ifid_glacflux) GLACFLUXG
  write (*,*) "GLACFLUXG finished"
  write (*,*) SHAPE(GLACFLUXG)
  write (*,*) GLACFLUXG
  
  write(*,*) "################################################################"
  write(*,*) GLACFLUXG(66,126)
  write(*,*) GLACFLUXG(126,66)
  
  close(ifid_glacflux)
  write (*,*) " *** END! ***"


  ! PG: Test 2 with ASCII version
  write (*,*) "*** ASCII TEST ***"
  ALLOCATE(GLACFLUX2(254, 220))
  glacfname2="FBFLOMPIOM.txt"
  open(ifid_glacflux2, file=trim(glacfname2), form="unformatted", status="old")
  write (*,*) "File Opened!"
  read(ifid_glacflux2) GLACFLUX2
  write (*,*) "Reading GLACFLUX2..."
  !write(*,*) GLACFLUX2
  write (*,*) "GLACFLUX2 finished!"
  close(ifid_glacflux2)
  write (*,*) "*** END! ***"

  ! PG: Test 3 as netcdf version
  write (*,*) "*** NETCDF TEST ***"
  ALLOCATE(GLACFLUX3(254, 220)) ! PG: This will be the array where glacflux goes
  fw_hosing_varname="area"
  glacfnamenc="FBFLOMPIOM.nc"   ! PG: The filename
  ! Open the file. NF90_NOWRITE says "read only, dude"
  call check(nf90_open(glacfnamenc, NF90_NOWRITE, ncid)) 
  ! get the varid
  call check(nf90_inq_varid(ncid, fw_hosing_varname, varid))
  ! ...and read the variable
  call check(nf90_get_var(ncid, varid, GLACFLUX3))
  write(*,*) "Here is the netcdf read GLACFLUX3:"
  write(*,*) GLACFLUX3(:,:)
  ! Close the file again:
  call check(nf90_close(ncid))
  write(*,*) "*** END! ***"


CONTAINS
  subroutine check(status)
      integer, intent ( in) :: status
      if(status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped"
      end if
  end subroutine check  
end program test
