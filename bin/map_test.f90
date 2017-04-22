program map_test

  use healpix_types
  use fourier, only: p_lm_gen, direct_fourier, inverse_fourier
  use map_coef, only: a_lm_gasdev
  use minkowski, only: area, length
  use distance, only: s2

  implicit none
  integer(kind=i8b) :: legendre_num = 5 !!! legendre_num = n_pix / 2
  integer(kind=i4b) :: n_pix = 512_i8b
  real(kind=dp), dimension(:, :), allocatable :: map
  complex(kind=dpc), dimension(:, :), allocatable :: a_lms
  integer(kind=i8b) :: i, j
  integer(kind=i8b) :: err_map, err_alms

  allocate(map(1:n_pix+1, 1:n_pix/2+1), stat=err_map)
  if (err_map /= 0) print *, "map: Allocation request denied"

  allocate(a_lms(0:legendre_num, 0:legendre_num), stat=err_alms)
  if (err_alms /= 0) print *, "a_lms: Allocation request denied"

  open(1, file="../../out.dat")

  call a_lm_gasdev(legendre_num, a_lms, 151514512, 0.d0, 1.d0)

  call direct_fourier(n_pix, map, legendre_num, a_lms)

  write(*, *) '0, 0', a_lms(0, 0)
  write(*, *) '0, 1', a_lms(0, 0)
  write(*, *) '1, 1', a_lms(1, 1)
  write(*, *) '0, 2', a_lms(0, 2)
  write(*, *) '1, 2', a_lms(1, 2)
  write(*, *) '2, 2', a_lms(2, 2)

  ! write(*, *) 'map(:, 2)', map(:, 2)
  ! write(*, *) 'map(:, 3)', map(:, 3)
  ! write(*, *) 'map(:, 4)', map(:, 4)

  call inverse_fourier(n_pix, map, legendre_num, a_lms)

  write(*, *) '0, 0', a_lms(0, 0)
  write(*, *) '0, 1', a_lms(0, 1)
  write(*, *) '1, 1', a_lms(1, 1)
  write(*, *) '0, 2', a_lms(0, 2)
  write(*, *) '1, 2', a_lms(1, 2)
  write(*, *) '2, 2', a_lms(2, 2)

  do i = 1, n_pix+1, 1
    do j = 1, n_pix/2+1, 1
      write(1, *) map(i, j), i, j
    end do
  end do

  write(*, *) 'area', area(n_pix, map, 0.d0)
  write(*, *) 'length', length(n_pix, map, 0.d0)


  if (allocated(a_lms)) deallocate(a_lms, stat=err_alms)
  if (err_alms /= 0) print *, "a_lms: Deallocation request denied"

  if (allocated(map)) deallocate(map, stat=err_map)
  if (err_map /= 0) print *, "map: Deallocation request denied"

end program map_test
