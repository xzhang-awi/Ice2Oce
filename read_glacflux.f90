program test
  implicit none
  INTEGER*8 :: IUGLAC1, IUGLAC2, IUGLAC3, IUGLAC4
  REAL*8, allocatable :: GLACFLUXG(:,:)
  CHARACTER (len=16) :: glacfname
  INTEGER, parameter :: ifid_glacflux=552
  
  ALLOCATE(GLACFLUXG(254, 220))
  write (*,*) "Hello, World!" 
  write (*,*) "I will read the file FBFLOMPIOM.ext8"
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
end program test
