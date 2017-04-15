program nr_test

  use healpix_types
  use pix_tools
  use alm_tools
  use fitstools, only: getsize_fits, input_map
  use, intrinsic :: iso_c_binding

  implicit none
  include 'fftw3.f03'

  integer :: npixtot, nmaps, nside, npix
  double precision, dimension(:, :), allocatable :: map
  integer :: err_map, err_dw8, err_alm

  integer :: lmax = 512
  double precision, dimension(2) :: z
  double precision, dimension(:, :), allocatable :: dw8
  complex, dimension(:, :, :), allocatable :: alm

  npixtot = getsize_fits("../../planck_data/COM_CMB_IQU-smica_1024_R2.02_full.fits", nmaps=nmaps, nside=nside)
  npix = nside2npix(nside)

  write(unit=*, fmt=*) 'npixtot', npixtot
  write(unit=*, fmt=*) 'npix', npix

  nmaps = 1

  allocate(map(0:npix-1,1:nmaps), stat=err_map)
  if (err_map /= 0) print *, "map(0:npix-1,1:nmaps): Allocation request denied"

  allocate(dw8(1:2*nside,1:3), stat=err_dw8)
  if (err_dw8 /= 0) print *, "dw8(1:2*nside,1:3): Allocation request denied"

  allocate(alm(1:nmaps, 0:lmax, 0:lmax), stat=err_alm)
  if (err_alm /= 0) print *, "alm(1:3, 0:lmax, 0:lmax): Allocation request denied"


  dw8 = 1.d0

  z = dsin(10.d0 * DEG2RAD)

  call input_map("../../planck_data/COM_CMB_IQU-smica_1024_R2.02_full.fits", map, npix, nmaps, fmissval=0.d0)

  call map2alm(nside, lmax, lmax, map, alm)

  write(unit=*, fmt=*) 'nmaps', nmaps
  write(unit=*, fmt=*) 'nside', nside
  write(unit=*, fmt=*) 'npix-1', npix-1
  write(unit=*, fmt=*) 'npixtot', npixtot
  write(unit=*, fmt=*) '12*nsmax**2-1', nside2npix(nside)-1


  if (allocated(alm)) deallocate(alm, stat=err_alm)
  if (err_alm /= 0) print *, "alm(1:3, 0:lmax, 0:lmax): Deallocation request denied"

  if (allocated(dw8)) deallocate(dw8, stat=err_dw8)
  if (err_dw8 /= 0) print *, "dw8(1:2*nside,1:3): Deallocation request denied"

  if (allocated(map)) deallocate(map, stat=err_map)
  if (err_map /= 0) print *, "map: Deallocation request denied"

end program nr_test
