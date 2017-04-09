program nr_test

  use pix_tools
  use, intrinsic :: iso_c_binding
  use nr, only: gasdev
  use ran_state, only: ran_seed

  implicit none
  include 'fftw3.f03'
  real :: rand = 1.d0
  integer :: first, second, third

  double precision, dimension(0:3) :: a_coef = (/1.d0, 2.d0, 3.d0, 4.d0/)
  double precision, dimension(0:3) :: b_coef = (/3.d0, 2.d0, 1.d0, 5.d0/)
  complex, dimension(0:3) :: a_b

  type(C_PTR) :: plan
  complex(C_DOUBLE_COMPLEX), dimension(0:3) :: in, out

  plan = fftw_plan_dft_1d(4, in, out, FFTW_FORWARD, FFTW_ESTIMATE)

  in(0) = (0.d0, 2.d0)
  in(1) = (5.d0, 1.d0)
  in(2) = (0.d0, 0.d0)
  in(3) = (0.d0, 0.d0)

  write(unit=*, fmt=*) in

  call fftw_execute_dft(plan, in, out)

  write(unit=*, fmt=*) out

  call fftw_destroy_plan(plan)

  plan = fftw_plan_dft_1d(4, in, out, FFTW_BACKWARD, FFTW_ESTIMATE) ! FFTW_MEASURE

  call fftw_execute_dft(plan, out, in)

  write(unit=*, fmt=*) in / 4.d0

  call fftw_destroy_plan(plan)

  first = 1203124125
  second = 1234145
  third = first + second

  call ran_seed(third)
  call gasdev(rand)
  write(unit=*, fmt=*) rand
  call gasdev(rand)
  write(unit=*, fmt=*) rand

  write(unit=*, fmt=*)

  call ran_seed(first)
  call gasdev(rand)
  write(unit=*, fmt=*) rand

  call ran_seed(second)
  call gasdev(rand)
  write(unit=*, fmt=*) rand

  a_b = cmplx(a_coef, b_coef)

  write(unit=*, fmt=*) 'a_coef'
  write(unit=*, fmt=*) a_coef
  write(unit=*, fmt=*) 'b_coef'
  write(unit=*, fmt=*) b_coef
  write(unit=*, fmt=*) 'a_b'
  write(unit=*, fmt=*) a_b
  write(unit=*, fmt=*) 'real(a_b)'
  write(unit=*, fmt=*) real(a_b)
  write(unit=*, fmt=*) 'aimag(a_b)'
  write(unit=*, fmt=*) aimag(a_b)


end program nr_test
