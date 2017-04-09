program plm_test

  use healpix_types
  use alm_tools, only: plm_gen
  use pix_tools
  use, intrinsic :: iso_c_binding
  use fourier, only: p_lm_gen

  implicit none
  double precision, dimension(:, :), allocatable :: polynom
  integer :: legendre_num = 4
  integer :: err_polynom
  double precision :: theta
  integer :: pix_num = 8

  include 'fftw3.f03'

  integer(I4B) :: nside, lmax, mmax, n_plm, npix
  real(DP), dimension(:,:), allocatable :: plm

  nside=256
  lmax=512
  mmax=lmax
  npix=nside2npix(nside)

  n_plm=nside*(mmax+1)*(2*lmax-mmax+2)
  allocate(plm(0:n_plm-1,1:3))

  call plm_gen(nside, lmax, mmax, plm)


  allocate(polynom(0:legendre_num, 0:legendre_num), stat=err_polynom)
  if (err_polynom /= 0) print *, "polynom: Allocation request denied"

  theta = 2.d0 * PI * (2 - 1) / pix_num

  call p_lm_gen(theta,legendre_num,polynom)

  write(*, *) polynom(:, 2)

  if (allocated(polynom)) deallocate(polynom, stat=err_polynom)
  if (err_polynom /= 0) print *, "polynom: Deallocation request denied"

end program plm_test
