program map_test

  use fourier, only: p_lm_gen, direct_fourier
  use map_coef, only: a_lm_gasdev

  implicit none
  integer :: legendre_num = 3
  complex, dimension(0:3, 0:3) :: a_lms
  integer :: n_pix = 8
  double precision, dimension(1:9, 1:5) :: map
  integer :: i, j

  open(1, file="../../out.dat")

  call a_lm_gasdev(legendre_num, a_lms, 151512526, 0.d0, 1.d0)

  call direct_fourier(n_pix, map, legendre_num, a_lms)

  do i = 1, n_pix+1, 1
    do j = 1, n_pix/2+1, 1
      write(1, *) map(i, j), i, j
    end do
  end do

end program map_test
