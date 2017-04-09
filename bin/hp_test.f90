program test

  use pix_tools
  use fitstools, only: getsize_fits, input_map
  use nrtype
  use nrutil
  use ran_state

  IMPLICIT NONE
  DOUBLE PRECISION theta,phi
  integer nside,ipix

  nside=2048
  theta=0.4543
  phi=1.75574

  call ang2pix_ring(nside,theta,phi,ipix)
  write(*,*) 'pixel=', ipix

end program
