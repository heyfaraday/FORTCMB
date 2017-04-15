program map_test2

  use healpix_types

  use map_coef, only: a_lm_gasdev
  use distance, only: s2
  use minkowski, only: area
  use fourier, only: p_lm_gen, direct_fourier, inverse_fourier, direct_point_fourier
  use minkowski, only: area

  use nr, only: gasdev
  use ran_state, only: ran_seed

  implicit none

  ! Precision test
  ! integer(kind=i8b) :: num_i8b = max_i8b
  ! real(kind=dp) :: num_dp = max_dp
  ! complex(kind=dpc) :: num_dpc = (max_dp, max_dp)

  ! write(*, *) 'i8b', num_i8b
  ! write(*, *) 'dp', num_dp
  ! write(*, *) 'dpc', num_dpc

  !
  integer(kind=i8b) :: legendre_num = 1_i8b
  integer(kind=i4b) :: n_pix = 8192_i4b
  real(kind=dp), dimension(:, :), allocatable :: map
  complex(kind=dpc), dimension(:, :), allocatable :: a_lm

  integer(kind=i8b) :: err_map = 0_i8b
  integer(kind=i8b) :: err_a_lm = 0_i8b
  !
  integer(kind=i8b) :: i, j
  real(kind=dp) :: norm, fun

  !
  ! real(kind=sp) :: rand_num_dummy
  ! real(kind=dp) :: rand_num
  !

  !
  allocate(map(1:n_pix+1, 1:n_pix/2+1), stat=err_map)
  if (err_map /= 0) print *, "map: Allocation request denied"

  allocate(a_lm(0:legendre_num, 0:legendre_num), stat=err_a_lm)
  if (err_a_lm /= 0) print *, "a_lm: Allocation request denied"
  !

  !
  map = 0.0_dp
  a_lm = 0.0_dp
  !

  ! Random precision test
  ! call ran_seed(125616)
  ! call gasdev(rand_num_dummy)
  ! rand_num = rand_num_dummy
  ! write(*, *) 'rand_num_dummy', rand_num_dummy
  ! write(*, *) 'rand_num', rand_num
  !

  !
  call a_lm_gasdev(legendre_num, a_lm, 12451525_i4b, 0.0_dp, 1.0_dp)
  write(*, *) '0, 0', a_lm(0, 0)
  write(*, *) '0, 1', a_lm(0, 1)
  write(*, *) '1, 1', a_lm(1, 1)
  write(*, *) '0, 2', a_lm(0, 2)
  write(*, *) '1, 2', a_lm(1, 2)
  write(*, *) '2, 2', a_lm(2, 2)
  write(*, *) '0, 3', a_lm(0, 3)
  write(*, *) '1, 3', a_lm(1, 3)
  write(*, *) '2, 3', a_lm(2, 3)
  write(*, *) '3, 3', a_lm(3, 3)
  !

  !
  call direct_fourier(n_pix, map, legendre_num, a_lm)
  write(*, *) 'j=1', map(:, 1) * map(:, 1)
  write(*, *) 'j=2', map(:, 2) * map(:, 2)
  write(*, *) 'j=3', map(:, 3) * map(:, 3)
  write(*, *) 'j=4', map(:, 4) * map(:, 4)
  write(*, *) 'j=5', map(:, 5) * map(:, 5)

  write(*, *) 'point', direct_point_fourier(2, 2, n_pix, legendre_num, a_lm)
  write(*, *) 'map', map(2, 2)
  !
  call inverse_fourier(n_pix, map, legendre_num, a_lm)
  write(*, *) '0, 0', a_lm(0, 0)
  write(*, *) '0, 1', a_lm(0, 1)
  write(*, *) '1, 1', a_lm(1, 1)
  write(*, *) '0, 2', a_lm(0, 2)
  write(*, *) '1, 2', a_lm(1, 2)
  write(*, *) '2, 2', a_lm(2, 2)
  write(*, *) '0, 3', a_lm(0, 3)
  write(*, *) '1, 3', a_lm(1, 3)
  write(*, *) '2, 3', a_lm(2, 3)
  write(*, *) '3, 3', a_lm(3, 3)
  !

  !
  write(*, *) 'area', area(n_pix, map, 0.0_dp)
  !

  fun = 0.0_dp
  norm = 0.0_dp

  write(*, *) '-------------'

  !
  do i = 1, n_pix, 1
    do j = 1, n_pix/2 + 1 , 1
      fun = fun + map(i, j) * map(i, j) * dsin(2.0_dp * PI * (j - 1) / n_pix)
      norm = norm + dsin(2.0_dp * PI * (j - 1) / n_pix)
    end do
  end do

  write(*, *) fun, norm, fun * 4.0_dp * PI
  !

  !
  if (allocated(map)) deallocate(map, stat=err_map)
  if (err_map /= 0) print *, "map: Deallocation request denied"

  if (allocated(a_lm)) deallocate(a_lm, stat=err_a_lm)
  if (err_a_lm /= 0) print *, "a_lm: Deallocation request denied"
  !

end program map_test2
