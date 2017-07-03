program test
  implicit none
  INTEGER*8 :: IUGLAC1, IUGLAC2, IUGLAC3, IUGLAC4
  REAL*8, allocatable :: GLACFLUXG(:,:), GLACFLUX2(:,:)
  CHARACTER (len=16) :: glacfname
  CHARACTER (len=16) :: glacfname2
  INTEGER, parameter :: ifid_glacflux=552
  INTEGER, parameter :: ifid_glacflux2=553
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
  ! write (*,*) GLACFLUXG
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


end program test
